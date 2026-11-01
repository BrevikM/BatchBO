################################################################################
# sobol.R -- global sensitivity analysis of fitted GP surfaces
################################################################################


#================== Stored Sobol Analysis ==================

#runSobol() estimates first-order, total-order and, optionally, pairwise
#second-order Sobol indices of a stored surrogate's posterior-mean surface, over
#a uniform distribution on each coded factor. sensobol builds the integration
#design (shifted Sobol' points) and computes the indices; it also checks its own
#arguments (N, conf, first, total, matrices and interval_type). Bootstrap
#intervals (off by default) cover numerical integration only, not uncertainty in
#the fitted GP. With scale = "response", the mean is decoded to the reported
#response scale first; this changes the indices only for a log response.
#
#Reference: Puy, A., Lo Piano, S., Saltelli, A., and Levin, S. A. (2022).
#sensobol: An R Package to Compute Variance-Based Sensitivity Indices.
#Journal of Statistical Software, 102(5), 1-37. doi:10.18637/jss.v102.i05

#' @export
runSobol <- function(
    x,
    name = "sobol",
    surrogate_name = NULL,
    N = 4096L,
    bootstrap = 0L,
    conf = 0.95,
    second_order = FALSE,
    scale = c("model", "response"),
    seed = NULL,
    overwrite = FALSE,
    first = "saltelli",
    total = "jansen",
    interval_type = "percent",
    matrices = c("A", "B", "AB")
) {
  if (!methods::is(x, "BOCampaign")) {
    stop('"x" must be a BOCampaign.', call. = FALSE)
  }

  name <- .validate_name(name)
  .validate_overwrite_name(name, names(x@analyses), overwrite, what = "analysis")
  surrogate_name <- .resolve_surrogate_name(x, surrogate_name)
  bootstrap <- .validate_count(bootstrap, "bootstrap", minimum = 0L)
  second_order <- .validate_flag(second_order, "second_order")
  seed <- .validate_seed(seed)
  scale <- match.arg(scale)

  fit <- x@surrogates[[surrogate_name]]$fit
  labels <- x@design@label
  D <- length(labels)
  parameter_names <- paste0("x", seq_len(D))
  pairs <- if (second_order && D > 1L) utils::combn(D, 2L, simplify = FALSE) else list()
  sobol_order <- if (length(pairs)) "second" else "first"

  evaluate_surface <- function(X) {
    value <- gp_predict(fit, X)$mean
    if (identical(scale, "response")) value <- .decode_fit_response(x, fit, value)
    value
  }

  #Integration design from sensobol: Cranley-Patterson shifted Sobol' points,
  #stacked as the blocks that sobol_indices() expects.
  design <- sensobol::sobol_matrices(
    matrices = matrices,
    N = N,
    params = parameter_names,
    order = sobol_order,
    type = "QRN",
    scrambling = "shift",
    seed = seed
  )

  #Evaluate the surface block by block (N rows at a time) to bound memory use.
  blocks <- split(seq_len(nrow(design)), (seq_len(nrow(design)) - 1L) %/% N)
  Y <- unlist(
    lapply(blocks, function(rows) evaluate_surface(design[rows, , drop = FALSE])),
    use.names = FALSE
  )

  result <- .with_seed(
    .seed_offset(seed, 1L),
    sensobol::sobol_indices(
      matrices = matrices,
      Y = Y,
      N = N,
      params = parameter_names,
      first = first,
      total = total,
      order = sobol_order,
      boot = bootstrap > 0L,
      R = if (bootstrap > 0L) bootstrap else NULL,
      conf = conf,
      type = interval_type
    )
  )
  results <- as.data.frame(result$results)

  #One column of sensobol's results for one kind of index, in the given order.
  #Interval limits exist only when bootstrapped, and are NA otherwise.
  index <- function(kind, names, column) {
    rows <- results[results$sensitivity == kind, , drop = FALSE]
    value <- rows[[column]][match(names, rows$parameters)]
    if (is.null(value)) rep(NA_real_, length(names)) else value
  }
  pair_names <- vapply(pairs, function(pair) paste(parameter_names[pair], collapse = "."), character(1))

  table <- data.frame(
    factor = labels,
    first_order = index("Si", parameter_names, "original"),
    first_lower = index("Si", parameter_names, "low.ci"),
    first_upper = index("Si", parameter_names, "high.ci"),
    total_order = index("Ti", parameter_names, "original"),
    total_lower = index("Ti", parameter_names, "low.ci"),
    total_upper = index("Ti", parameter_names, "high.ci"),
    stringsAsFactors = FALSE
  )

  interaction_table <- data.frame(
    factor_1 = labels[vapply(pairs, `[`, integer(1), 1L)],
    factor_2 = labels[vapply(pairs, `[`, integer(1), 2L)],
    second_order = index("Sij", pair_names, "original"),
    second_lower = index("Sij", pair_names, "low.ci"),
    second_upper = index("Sij", pair_names, "high.ci"),
    stringsAsFactors = FALSE
  )

  if (anyNA(table$first_order) || anyNA(table$total_order) || anyNA(interaction_table$second_order)) {
    stop("sensobol returned missing indices.", call. = FALSE)
  }

  settings <- list(
    N = N,
    bootstrap = bootstrap,
    conf = conf,
    second_order = second_order,
    scale = scale,
    seed = seed,
    sampling = "sensobol QRN, Cranley-Patterson shift",
    matrices = matrices,
    first = first,
    total = total,
    interval_type = if (bootstrap > 0L) interval_type else NULL,
    sensobol_version = as.character(utils::packageVersion("sensobol"))
  )

  x@analyses[[name]] <- list(
    table = table,
    interactions = interaction_table,
    surrogate = surrogate_name,
    settings = settings,
    provenance = .provenance_analysis(
      x = x,
      surrogate_name = surrogate_name,
      kind = "sobol",
      settings = settings
    ),
    call = match.call()
  )

  methods::validObject(x)
  x
}


#================== Accessors ==================

#sobolIndices() returns the factor table, sobolInteractions() the pair table,
#and with details = TRUE either returns the complete stored entry.

#' @export
sobolIndices <- function(x, name = NULL, strict = FALSE, details = FALSE) {
  entry <- .analysis_entry(x = x, name = name, kind = "sobol", strict = strict)
  if (.validate_flag(details, "details")) entry else entry$table
}


#' @export
sobolNames <- function(x) {
  .analysis_names(x, "sobol")
}


#' @export
sobolInteractions <- function(x, name = NULL, strict = FALSE, details = FALSE) {
  entry <- .analysis_entry(x = x, name = name, kind = "sobol", strict = strict)
  if (!nrow(entry$interactions)) {
    stop(
      "Second-order indices were not calculated for this analysis; ",
      "rerun runSobol() with second_order = TRUE.",
      call. = FALSE
    )
  }
  if (.validate_flag(details, "details")) entry else entry$interactions
}


#================== Sobol-Index Plot ==================

#plotSobol() draws first- and total-order indices per factor, either stacked
#(first order, then the interaction part S_Ti - S_i on top) or side by side,
#with factors in decreasing order of total-order index.
#Further arguments (main, xlab, ylab, ylim, ...) go to plot(); line widths and
#other settings follow par(). Both plotters return the plotted table invisibly.

#' @export
plotSobol <- function(
    x,
    name = NULL,
    strict = FALSE,
    stacked = TRUE,
    colors = c(first_order = "#E76F61", total_order = "#00AEB5"),
    border = "#333333",
    show_intervals = TRUE,
    legend = TRUE,
    las = 2L,
    ...
) {
  las <- .validate_las(las)

  table <- sobolIndices(x, name = name, strict = strict)
  table <- table[order(table$total_order, decreasing = TRUE), , drop = FALSE]
  labels <- gsub("_", " ", table$factor, fixed = TRUE)
  at <- seq_len(nrow(table))

  if (stacked) {
    #The bar reaches the total index, and the first-order part is clipped to
    #[0, total] for display only; Monte Carlo estimates can fall slightly outside.
    top <- pmax(table$total_order, 0)
    first <- pmin(pmax(table$first_order, 0), top)
    left <- list(at = at, bottom = 0, top = first)
    right <- list(at = at, bottom = first, top = top)
    width <- 0.66
    offset <- 0.10
  } else {
    left <- list(at = at - 0.20, bottom = 0, top = table$first_order)
    right <- list(at = at + 0.20, bottom = 0, top = table$total_order)
    width <- 0.34
    offset <- 0.20
  }

  values <- c(0, left$top, right$top)
  if (show_intervals) {
    values <- c(values, table$first_lower, table$first_upper, table$total_lower, table$total_upper)
  }
  #Bars start at zero unless an estimate or interval is negative.
  limits <- .padded_range(values, 0.07)
  if (min(values, na.rm = TRUE) >= 0) limits[1L] <- 0

  #Widen the bottom margin for long rotated factor names.
  axis_cex <- 0.86
  margins <- graphics::par("mai")
  if (las %in% c(2L, 3L)) {
    margins[1L] <- max(margins[1L], max(graphics::strwidth(labels, "inches", cex = axis_cex)) + 0.75)
  }
  old_margins <- graphics::par(mai = margins)
  on.exit(graphics::par(old_margins), add = TRUE)

  do.call(graphics::plot, utils::modifyList(
    list(
      x = NA_real_, y = NA_real_,
      xlim = c(0.4, length(at) + 0.6), ylim = limits,
      xaxt = "n", xlab = "", ylab = "Sobol index", bty = "l"
    ),
    list(...)
  ))
  graphics::axis(1L, at = at, labels = labels, las = las, cex.axis = axis_cex)
  graphics::abline(h = graphics::axTicks(2L), col = "#E8E8E8")
  graphics::abline(h = 0, col = "#A0A0A0", lty = 3L)

  graphics::rect(left$at - width / 2, left$bottom, left$at + width / 2, left$top,
                 col = colors[["first_order"]], border = border)
  graphics::rect(right$at - width / 2, right$bottom, right$at + width / 2, right$top,
                 col = colors[["total_order"]], border = border)

  #Error bars with caps; rows without intervals (NA) are skipped.
  if (show_intervals) {
    graphics::arrows(at - offset, table$first_lower, at - offset, table$first_upper,
                     angle = 90, code = 3, length = 0.04, col = border)
    graphics::arrows(at + offset, table$total_lower, at + offset, table$total_upper,
                     angle = 90, code = 3, length = 0.04, col = border)
  }

  if (legend) {
    graphics::legend(
      "topright",
      legend = c(
        expression("First order " * (S[i])),
        if (stacked) expression("Interactions " * (S[T[i]] - S[i])) else expression("Total order " * (S[T[i]]))
      ),
      fill = unname(colors[c("first_order", "total_order")]),
      border = border,
      bty = "n",
      cex = 0.9
    )
  }

  invisible(table)
}


#plotSobolInteractions() draws the pairwise second-order indices as horizontal
#bars, or as points with dots = TRUE, with intervals when they were bootstrapped,
#the largest index at the top.

#' @export
plotSobolInteractions <- function(
    x,
    name = NULL,
    strict = FALSE,
    dots = FALSE,
    color = "#00AEB5",
    interval_color = "#333333",
    pch = 19L,
    show_intervals = TRUE,
    las = 2L,
    ...
) {
  las <- .validate_las(las)

  table <- sobolInteractions(x, name = name, strict = strict)
  table <- table[order(table$second_order), , drop = FALSE]
  labels <- paste(
    gsub("_", " ", table$factor_1, fixed = TRUE),
    "x",
    gsub("_", " ", table$factor_2, fixed = TRUE)
  )
  at <- seq_len(nrow(table))

  values <- c(0, table$second_order)
  if (show_intervals) values <- c(values, table$second_lower, table$second_upper)
  limits <- .padded_range(values, 0.08)

  #Widen the left margin for long horizontal pair labels.
  axis_cex <- 0.76
  margins <- graphics::par("mai")
  if (las %in% c(1L, 2L)) {
    margins[2L] <- max(margins[2L], max(graphics::strwidth(labels, "inches", cex = axis_cex)) + 0.45)
  }
  old_margins <- graphics::par(mai = margins)
  on.exit(graphics::par(old_margins), add = TRUE)

  do.call(graphics::plot, utils::modifyList(
    list(
      x = NA_real_, y = NA_real_,
      xlim = limits, ylim = c(0.5, length(at) + 0.5),
      yaxt = "n", xlab = expression("Second-order Sobol index " * (S[ij])), ylab = "", bty = "l"
    ),
    list(...)
  ))
  graphics::axis(2L, at = at, labels = labels, las = las, cex.axis = axis_cex)
  graphics::abline(v = 0, col = "#A0A0A0", lty = 3L)
  graphics::abline(h = at, col = "#EEEEEE")

  if (!dots) {
    graphics::rect(0, at - 0.33, table$second_order, at + 0.33, col = color, border = interval_color)
    if (show_intervals) {
      graphics::arrows(table$second_lower, at, table$second_upper, at,
                       angle = 90, code = 3, length = 0.04, col = interval_color)
    }
  } else {
    if (show_intervals) {
      graphics::segments(table$second_lower, at, table$second_upper, at, col = interval_color)
    }
    graphics::points(table$second_order, at, pch = pch, col = color)
  }

  invisible(table)
}
