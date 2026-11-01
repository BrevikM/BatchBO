################################################################################
# effects.R -- latent main effects and two-factor response surfaces
################################################################################


#================== Shared Helpers ==================

#Validate one or more factor labels against the declared design space.
.effect_factors <- function(x, factors = NULL, required = NULL) {
  if (is.null(factors)) {
    factors <- x@design@label
  }

  if (
    !is.character(factors) ||
    !length(factors) ||
    anyNA(factors) ||
    any(!nzchar(trimws(factors))) ||
    anyDuplicated(factors)
  ) {
    stop(
      '"factors" must contain unique, non-empty factor labels.',
      call. = FALSE
    )
  }

  factors <- trimws(factors)
  unknown <- setdiff(factors, x@design@label)
  if (length(unknown)) {
    stop(
      "Unknown factor label(s): ",
      paste(unknown, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (!is.null(required) && length(factors) != required) {
    stop(
      '"factors" must contain exactly ',
      required,
      " factor labels.",
      call. = FALSE
    )
  }

  factors
}


#Validate one grid size or one grid size per varying factor.
.effect_grid_size <- function(grid, length_required) {
  if (
    !is.numeric(grid) ||
    is.logical(grid) ||
    !length(grid) ||
    anyNA(grid) ||
    any(!is.finite(grid)) ||
    any(grid != floor(grid)) ||
    any(grid < 2L) ||
    !length(grid) %in% c(1L, length_required)
  ) {
    stop(
      '"grid" must contain one integer of at least two, or one per factor.',
      call. = FALSE
    )
  }

  rep(as.integer(grid), length.out = length_required)
}


#Transform physical factor values to their natural plotting coordinates. This
#preserves the spacing used by the GP for transformed factors while axis labels
#remain in the original experimental units.
.effect_plot_coordinates <- function(x, factor, values) {
  index <- match(factor, x@design@label)
  transform <- x@design@transform[index]

  if (identical(transform, "log10")) {
    log10(values)
  } else {
    values
  }
}


#Draw an axis in physical units, using logarithmic spacing when declared by the
#design space.
.effect_axis <- function(
    x,
    factor,
    values,
    side,
    las = graphics::par("las")
) {
  index <- match(factor, x@design@label)
  transform <- x@design@transform[index]

  if (identical(transform, "log10")) {
    ticks <- grDevices::axisTicks(
      log10(range(values, finite = TRUE)),
      log = TRUE,
      nint = 5L
    )
    graphics::axis(
      side,
      at = log10(ticks),
      labels = format(
        ticks,
        trim = TRUE,
        scientific = FALSE,
        big.mark = ","
      ),
      las = las
    )
  } else {
    graphics::axis(side, las = las)
  }
}


#Range of the finite values, padded by a fraction of its width (or of its size
#when the values are all equal). Shared by the package's plots.
.padded_range <- function(values, fraction = 0.06) {
  limits <- range(values, finite = TRUE)
  padding <- diff(limits) * fraction
  if (!is.finite(padding) || padding == 0) {
    padding <- max(abs(limits), 1) * fraction
  }
  limits + c(-padding, padding)
}


#Average the latent posterior mean at every row of `varying_grid` over a fixed
#Sobol background. Reusing the background across grid locations removes Monte
#Carlo jitter from comparisons along an effect curve or surface. Averaging is
#performed on the standardized latent scale and decoded once afterwards, so a
#log-response analysis returns exp(E[f]) rather than E[exp(f)].
.marginal_latent_mean <- function(
    x,
    fit,
    varying_indices,
    varying_grid,
    background,
    maximum_evaluations = 20000L
) {
  varying_grid <- as.matrix(varying_grid)
  background <- as.matrix(background)
  D <- ncol(fit$X)
  M <- nrow(varying_grid)

  if (length(varying_indices) == D) {
    coded <- matrix(0.5, nrow = M, ncol = D)
    coded[, varying_indices] <- varying_grid
    return(.decode_fit_response(x, fit, gp_predict(fit, coded)$mean))
  }

  N <- nrow(background)
  cells_per_chunk <- max(
    1L,
    floor(as.integer(maximum_evaluations) / N)
  )
  output <- numeric(M)

  starts <- seq.int(1L, M, by = cells_per_chunk)
  for (first in starts) {
    cells <- seq.int(first, min(M, first + cells_per_chunk - 1L))
    coded <- background[
      rep(seq_len(N), times = length(cells)),
      ,
      drop = FALSE
    ]

    for (j in seq_along(varying_indices)) {
      coded[, varying_indices[j]] <- rep(
        varying_grid[cells, j],
        each = N
      )
    }

    prediction <- gp_predict(fit, coded)$mean
    output[cells] <- colMeans(
      matrix(
        prediction,
        nrow = N,
        ncol = length(cells)
      )
    )
  }

  .decode_fit_response(x, fit, output)
}


#Return the physical and coded condition used to define a conditional slice.
.response_slice_condition <- function(
    x,
    at,
    recommendation_name = NULL
) {
  D <- length(x@design@label)

  if (is.character(at) && length(at) == 1L && !is.na(at)) {
    at <- match.arg(at, c("centre", "recommendation"))

    if (identical(at, "centre")) {
      coded <- matrix(
        0.5,
        nrow = 1L,
        ncol = D,
        dimnames = list(NULL, x@design@label)
      )
      physical <- .decode(x@design, coded)
    } else {
      details <- recommendation(
        x,
        name = recommendation_name,
        strict = TRUE,
        details = TRUE
      )
      coded <- details$coded
      physical <- details$condition[
        ,
        x@design@label,
        drop = FALSE
      ]
    }
  } else {
    physical <- as.data.frame(at, stringsAsFactors = FALSE)
    if (nrow(physical) != 1L) {
      stop(
        '"at" must be "centre", "recommendation", or one physical condition.',
        call. = FALSE
      )
    }

    physical <- .validate_numeric_columns(
      physical,
      x@design@label,
      '"at"'
    )
    .validate_log_domain(physical, x@design, '"at"')
    coded <- .encode(x@design, physical)
  }

  if (any(coded < -1e-8 | coded > 1 + 1e-8)) {
    stop(
      'The response-surface slice condition must lie inside the design space.',
      call. = FALSE
    )
  }

  list(
    physical = physical[, x@design@label, drop = FALSE],
    coded = coded[, x@design@label, drop = FALSE]
  )
}


#================== Stored Main Effects ==================

#runMainEffects() calculates one-factor marginal response profiles E[f | x_i]
#from a stored surrogate's posterior mean. At every plotted factor value,
#predictions are averaged over the same randomized Sobol background for all
#remaining factors, i.e. independent uniform variation on the coded design
#space, as in runSobol(). Averaging is on the transformed latent scale and
#decoded once, so for a log response the curve is exp(E[f]) rather than
#E[exp(f)]. mainEffects() returns the long table (or the stored entry with
#details = TRUE); plotMainEffects() draws one panel per factor.

#' @export
runMainEffects <- function(
    x,
    name = "main_effects",
    surrogate_name = NULL,
    factors = NULL,
    grid = 101L,
    N = 1024L,
    seed = NULL,
    overwrite = FALSE
) {
  if (!methods::is(x, "BOCampaign")) {
    stop('"x" must be a BOCampaign.', call. = FALSE)
  }

  name <- .validate_name(name)
  .validate_overwrite_name(name, names(x@analyses), overwrite, what = "analysis")
  surrogate_name <- .resolve_surrogate_name(x, surrogate_name)
  factors <- .effect_factors(x, factors)
  grid <- .effect_grid_size(grid, length(factors))
  N <- .validate_count(N, "N", minimum = 1L)
  seed <- .validate_seed(seed)

  fit <- x@surrogates[[surrogate_name]]$fit
  D <- ncol(fit$X)
  background <- .sbl(
    n = N,
    D = D,
    lower = rep(0, D),
    upper = rep(1, D),
    seed = seed
  )

  tables <- lapply(seq_along(factors), function(i) {
    factor <- factors[i]
    factor_index <- match(factor, x@design@label)
    coded_value <- seq(0, 1, length.out = grid[i])
    coded_grid <- matrix(coded_value, ncol = 1L)

    latent_mean <- .marginal_latent_mean(
      x = x,
      fit = fit,
      varying_indices = factor_index,
      varying_grid = coded_grid,
      background = background
    )

    physical_grid <- matrix(
      0.5,
      nrow = length(coded_value),
      ncol = D,
      dimnames = list(NULL, x@design@label)
    )
    physical_grid[, factor_index] <- coded_value
    physical_value <- .decode(
      x@design,
      physical_grid
    )[[factor]]

    data.frame(
      factor = factor,
      grid_index = seq_along(coded_value),
      coded_value = coded_value,
      factor_value = physical_value,
      latent_mean = latent_mean,
      stringsAsFactors = FALSE
    )
  })

  table <- do.call(rbind, tables)
  rownames(table) <- NULL
  settings <- list(
    factors = factors,
    grid = stats::setNames(grid, factors),
    N = N,
    seed = seed,
    input_distribution = "independent uniform in the coded design space",
    response_definition = "batch-free latent posterior mean on the response scale"
  )

  x@analyses[[name]] <- list(
    table = table,
    surrogate = surrogate_name,
    settings = settings,
    provenance = .provenance_analysis(
      x = x,
      surrogate_name = surrogate_name,
      kind = "main_effects",
      settings = settings
    ),
    call = match.call()
  )

  methods::validObject(x)
  x
}


#' @export
mainEffects <- function(x, name = NULL, strict = FALSE, details = FALSE) {
  entry <- .analysis_entry(x = x, name = name, kind = "main_effects", strict = strict)
  if (.validate_flag(details, "details")) entry else entry$table
}


#' @export
mainEffectNames <- function(x) {
  .analysis_names(x, "main_effects")
}


#================== Main-Effect Plot ==================

#One panel per factor, with rug marks at the completed sample settings when
#show_design = TRUE. Returns the plotted subset invisibly.

#' @export
plotMainEffects <- function(
    x,
    name = NULL,
    factors = NULL,
    strict = FALSE,
    columns = 2L,
    shared_y = TRUE,
    show_design = TRUE,
    main = NULL,
    panel_titles = FALSE,
    las = 0L,
    ylab = NULL,
    color = "#12355B",
    lwd = 2.5
) {
  table <- mainEffects(x, name = name, strict = strict)
  available <- unique(table$factor)
  factors <- if (is.null(factors)) {
    available
  } else {
    .effect_factors(x, factors)
  }

  missing_factors <- setdiff(factors, available)
  if (length(missing_factors)) {
    stop(
      "The stored analysis does not contain factor(s): ",
      paste(missing_factors, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  columns <- min(.validate_count(columns, "columns"), length(factors))
  rows <- ceiling(length(factors) / columns)
  las <- .validate_las(las)

  if (is.null(ylab)) {
    ylab <- x@response
  }

  selected <- table[table$factor %in% factors, , drop = FALSE]
  shared_limits <- .padded_range(selected$latent_mean)

  #Set the panel layout, and restore exactly these settings on exit.
  panel_margins <- graphics::par("mar")
  if (!panel_titles) {
    panel_margins[3L] <- min(panel_margins[3L], 0.8)
  }
  old_parameters <- graphics::par(
    mfrow = c(rows, columns),
    mar = panel_margins,
    oma = c(0, 0, if (is.null(main)) 0 else 2.2, 0)
  )
  on.exit(graphics::par(old_parameters), add = TRUE)

  for (factor in factors) {
    data <- table[table$factor == factor, , drop = FALSE]
    plot_x <- .effect_plot_coordinates(x, factor, data$factor_value)

    graphics::plot(
      plot_x,
      data$latent_mean,
      type = "n",
      xlab = factor,
      ylab = ylab,
      main = if (panel_titles) factor else NULL,
      ylim = if (shared_y) shared_limits else .padded_range(data$latent_mean),
      xaxt = "n",
      bty = "l"
    )
    .effect_axis(x, factor, data$factor_value, side = 1L, las = las)
    graphics::grid(nx = NA, ny = NULL, col = "#E6E9EC", lty = 1L)
    graphics::lines(plot_x, data$latent_mean, col = color, lwd = lwd)

    if (show_design && nrow(x@observations)) {
      sample_rows <- x@observations$Type == "Sample"
      graphics::rug(
        .effect_plot_coordinates(
          x,
          factor,
          x@observations[[factor]][sample_rows]
        ),
        col = grDevices::adjustcolor("#4A4A4A", alpha.f = 0.45),
        ticksize = 0.025
      )
    }
  }

  if (!is.null(main)) {
    graphics::mtext(main, side = 3L, outer = TRUE, line = 0.4, font = 2L)
  }

  invisible(selected)
}


#================== Stored Response Surfaces ==================

#runResponseSurfaces() evaluates the posterior-mean surface for every pair of
#the selected factors; two factors give one surface. A "marginal" surface
#averages over one shared uniform Sobol background for the unshown factors,
#on the latent scale as in runMainEffects(). A "slice" fixes the unshown factors
#at the coded centre, a stored recommendation, or a supplied physical
#condition. All values describe the batch-free latent posterior mean, not
#future batch-affected observations. responseSurfaces() returns a named list of
#surface tables (or the stored entry with details = TRUE).

#' @export
runResponseSurfaces <- function(
    x,
    name = "response_surfaces",
    surrogate_name = NULL,
    factors = NULL,
    type = c("marginal", "slice"),
    grid = 41L,
    N = 512L,
    at = "centre",
    recommendation_name = NULL,
    seed = NULL,
    overwrite = FALSE
) {
  if (!methods::is(x, "BOCampaign")) {
    stop('"x" must be a BOCampaign.', call. = FALSE)
  }

  name <- .validate_name(name)
  .validate_overwrite_name(name, names(x@analyses), overwrite, what = "analysis")
  surrogate_name <- .resolve_surrogate_name(x, surrogate_name)
  factors <- .effect_factors(x, factors)
  if (length(factors) < 2L) {
    stop("At least two factor labels are required.", call. = FALSE)
  }
  type <- match.arg(type)
  grid <- .effect_grid_size(grid, 2L)
  N <- .validate_count(N, "N", minimum = 1L)
  seed <- .validate_seed(seed)

  fit <- x@surrogates[[surrogate_name]]$fit
  labels <- x@design@label
  D <- length(labels)
  marginal <- identical(type, "marginal")

  #Marginal surfaces share one background; slices share one fixed condition.
  if (marginal) {
    background <- .sbl(n = N, D = D, lower = rep(0, D), upper = rep(1, D), seed = seed)
    slice_condition <- NULL
  } else {
    slice_condition <- .response_slice_condition(x, at, recommendation_name)
  }

  coded_grid <- as.matrix(expand.grid(
    seq(0, 1, length.out = grid[1L]),
    seq(0, 1, length.out = grid[2L]),
    KEEP.OUT.ATTRS = FALSE
  ))

  pairs <- utils::combn(factors, 2L, simplify = FALSE)
  surfaces <- lapply(pairs, function(pair) {
    indices <- match(pair, labels)

    latent_mean <- if (marginal) {
      .marginal_latent_mean(x, fit, indices, coded_grid, background)
    } else {
      evaluation <- matrix(slice_condition$coded, nrow = nrow(coded_grid), ncol = D, byrow = TRUE)
      evaluation[, indices] <- coded_grid
      .decode_fit_response(x, fit, gp_predict(fit, evaluation)$mean)
    }

    #Physical values of the two displayed factors on the grid.
    decode_grid <- matrix(0.5, nrow = nrow(coded_grid), ncol = D, dimnames = list(NULL, labels))
    decode_grid[, indices] <- coded_grid
    table <- .decode(x@design, decode_grid)[, pair, drop = FALSE]
    table$latent_mean <- latent_mean
    table
  })
  names(surfaces) <- vapply(pairs, paste, character(1), collapse = "__")

  settings <- list(
    factors = factors,
    pairs = pairs,
    type = type,
    grid = grid,
    N = if (marginal) N else NULL,
    seed = seed,
    at = at,
    recommendation_name = recommendation_name,
    input_distribution = if (marginal) {
      "independent uniform in the coded design space"
    } else {
      "unshown factors fixed at the stored slice condition"
    },
    response_definition = "batch-free latent posterior mean on the response scale"
  )

  x@analyses[[name]] <- list(
    surfaces = surfaces,
    surrogate = surrogate_name,
    settings = settings,
    slice_condition = slice_condition,
    provenance = .provenance_analysis(
      x = x,
      surrogate_name = surrogate_name,
      kind = "response_surfaces",
      settings = settings
    ),
    call = match.call()
  )

  methods::validObject(x)
  x
}


#' @export
responseSurfaces <- function(x, name = NULL, strict = FALSE, details = FALSE) {
  entry <- .analysis_entry(x = x, name = name, kind = "response_surfaces", strict = strict)
  if (.validate_flag(details, "details")) entry else entry$surfaces
}


#' @export
responseSurfacesNames <- function(x) {
  .analysis_names(x, "response_surfaces")
}


#================== Response-Surface Plot ==================

#Draw one surface. With color_key, filled.contour() adds a response key (and
#manages its own layout); otherwise image() draws a plain panel. Completed
#samples and the recommendation are overlaid when supplied.
.plot_response_surface <- function(
    x,
    table,
    observed,
    recommended,
    colors,
    contours,
    contour_color,
    zlim,
    color_key,
    point_legend,
    main
) {
  factors <- names(table)[1:2]
  x_values <- unique(table[[factors[1L]]])
  y_values <- unique(table[[factors[2L]]])
  plot_x <- .effect_plot_coordinates(x, factors[1L], x_values)
  plot_y <- .effect_plot_coordinates(x, factors[2L], y_values)
  z <- matrix(table$latent_mean, nrow = length(x_values), ncol = length(y_values))

  #A flat surface needs a slightly widened colour range.
  if (is.null(zlim)) zlim <- range(z, finite = TRUE)
  if (zlim[1L] == zlim[2L]) {
    zlim <- zlim + c(-1, 1) * max(abs(zlim[1L]), 1) * sqrt(.Machine$double.eps)
  }

  coordinates <- function(conditions) {
    list(
      x = .effect_plot_coordinates(x, factors[1L], conditions[[factors[1L]]]),
      y = .effect_plot_coordinates(x, factors[2L], conditions[[factors[2L]]])
    )
  }

  draw_details <- function() {
    .effect_axis(x, factors[1L], x_values, side = 1L)
    .effect_axis(x, factors[2L], y_values, side = 2L)

    if (diff(range(z, finite = TRUE)) > 0) {
      graphics::contour(plot_x, plot_y, z, add = TRUE, nlevels = contours, col = contour_color, lwd = 0.8)
    }
    if (!is.null(observed)) {
      graphics::points(coordinates(observed), pch = 21L, cex = 0.75, col = "#303030",
                       bg = grDevices::adjustcolor("white", alpha.f = 0.72))
    }
    if (!is.null(recommended)) {
      graphics::points(coordinates(recommended), pch = 23L, cex = 1.25, col = "#303030", bg = "#F2A900")
    }

    shown <- c(!is.null(observed), !is.null(recommended))
    if (point_legend && any(shown)) {
      graphics::legend(
        "topright",
        legend = c("Completed sample", "Recommendation")[shown],
        pch = c(21L, 23L)[shown],
        pt.bg = c("white", "#F2A900")[shown],
        col = "#303030",
        bty = "n",
        cex = 0.82
      )
    }
  }

  palette <- grDevices::colorRampPalette(colors)
  if (color_key) {
    graphics::filled.contour(
      x = plot_x,
      y = plot_y,
      z = z,
      zlim = zlim,
      color.palette = palette,
      nlevels = length(colors),
      xlab = factors[1L],
      ylab = factors[2L],
      main = main,
      key.title = graphics::title(main = x@response, cex.main = 0.8),
      plot.axes = draw_details()
    )
  } else {
    graphics::image(
      x = plot_x,
      y = plot_y,
      z = z,
      zlim = zlim,
      col = palette(length(colors)),
      axes = FALSE,
      xlab = factors[1L],
      ylab = factors[2L],
      main = main,
      useRaster = TRUE
    )
    draw_details()
    graphics::box()
  }
}


#plotResponseSurfaces() draws a single surface with its colour key (color_key)
#and point legend (point_legend), or several surfaces as panels on a shared
#colour range (shared_scale) or the supplied zlim. Returns the plotted tables
#invisibly.

#' @export
plotResponseSurfaces <- function(
    x,
    name = NULL,
    strict = FALSE,
    columns = 4L,
    shared_scale = TRUE,
    panel_titles = FALSE,
    show_observations = FALSE,
    recommendation_name = NULL,
    main = NULL,
    colors = grDevices::hcl.colors(64L, "Blues 3"),
    contours = 10L,
    contour_color = "#303030",
    color_key = TRUE,
    point_legend = TRUE,
    zlim = NULL
) {
  surfaces <- responseSurfaces(x, name = name, strict = strict)
  columns <- min(.validate_count(columns, "columns"), length(surfaces))

  sample_rows <- x@observations$Type == "Sample"
  observed <- if (show_observations && any(sample_rows)) x@observations[sample_rows, , drop = FALSE]
  recommended <- if (!is.null(recommendation_name)) {
    recommendation(x, name = recommendation_name, strict = TRUE)
  }
  if (is.null(zlim) && shared_scale) {
    zlim <- range(unlist(lapply(surfaces, `[[`, "latent_mean")), finite = TRUE)
  }

  draw <- function(table, color_key, point_legend, main) {
    .plot_response_surface(
      x, table, observed, recommended, colors, contours, contour_color,
      zlim, color_key, point_legend, main
    )
  }

  if (length(surfaces) == 1L) {
    draw(surfaces[[1L]], color_key, point_legend, main)
    return(invisible(surfaces))
  }

  #Several surfaces: a panel layout, restored on exit.
  rows <- ceiling(length(surfaces) / columns)
  old_parameters <- graphics::par(
    mfrow = c(rows, columns),
    mar = c(3.4, 3.6, if (panel_titles) 2.4 else 0.8, 0.7),
    mgp = c(2.0, 0.65, 0),
    cex = 0.78,
    oma = c(0, 0, if (is.null(main)) 0 else 2.2, 0)
  )
  on.exit(graphics::par(old_parameters), add = TRUE)

  for (table in surfaces) {
    draw(table, FALSE, FALSE, if (panel_titles) paste(names(table)[1:2], collapse = " and "))
  }
  for (panel in seq_len(rows * columns - length(surfaces))) {
    graphics::plot.new()
  }
  if (!is.null(main)) {
    graphics::mtext(main, side = 3L, outer = TRUE, line = 0.4, font = 2L)
  }

  invisible(surfaces)
}
