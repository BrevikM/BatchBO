################################################################################
# cumulative-best.R -- observed and batch-adjusted best-so-far trajectories
################################################################################


#================== Batch-Effect Adjustment ==================

#Posterior means (BLUPs) of the standardized batch effects:
#E[b | y] = sigma_b^2 Z^T V^{-1}(y - beta0*1).
.batch_effect_means <- function(g) {
  levels <- unique(as.character(g$batch))
  effect <- vapply(
    levels,
    function(level) {
      g$hp$sigma_b2 * sum(g$a[as.character(g$batch) == level])
    },
    numeric(1)
  )
  stats::setNames(effect, levels)
}


#================== Stored Cumulative Best ==================

#runCumulativeBest() follows the best tested condition as batches accumulate. At
#stage i only runs from batches 1, ..., i are considered, in the order the
#batches were added. "observed" is the running best of the measured responses;
#"batch_adjusted" first subtracts each run's estimated batch effect, using the
#batch effects of the chosen surrogate at every stage. Both are maxima of noisy
#values and therefore tend to overstate the true best slightly.
#
#The trajectory uses the surrogate's own data, so it covers the batches that
#surrogate was fitted to, and it is current exactly while that surrogate is.
#Sample runs always count; reference runs count only inside the domain, and only
#with include_references = TRUE. cumulativeBest() returns the stored table (or
#the stored entry with details = TRUE).

#' @export
runCumulativeBest <- function(
    x,
    name = "cumulative_best",
    surrogate_name = NULL,
    include_references = TRUE,
    overwrite = FALSE
) {
  if (!methods::is(x, "BOCampaign")) {
    stop('"x" must be a BOCampaign.', call. = FALSE)
  }

  name <- .validate_name(name)
  .validate_overwrite_name(name, names(x@analyses), overwrite, what = "analysis")
  surrogate_name <- .resolve_surrogate_name(x, surrogate_name)
  include_references <- .validate_flag(include_references, "include_references")

  entry <- x@surrogates[[surrogate_name]]
  g <- entry$fit
  type <- entry$provenance$state$type
  batch <- as.character(g$batch)
  batch_levels <- unique(batch)
  batch_number <- match(batch, batch_levels)
  eligible <- .eligible_runs(g$X, type, include_references)

  #Batch-adjusted response: each run minus the posterior mean of its batch effect.
  adjusted <- g$y - unname(.batch_effect_means(g)[batch])

  rows <- lapply(seq_along(batch_levels), function(i) {
    eligible_i <- which(eligible & batch_number <= i)
    observed_row <- eligible_i[which.max(g$y[eligible_i])]
    adjusted_row <- eligible_i[which.max(adjusted[eligible_i])]

    output <- data.frame(
      batch = batch_levels[i],
      batch_index = i,
      completed_runs = sum(batch_number <= i),
      eligible_runs = length(eligible_i),
      observed_best = .decode_fit_response(x, g, g$y[observed_row]),
      observed_best_source_batch = batch[observed_row],
      observed_best_type = type[observed_row],
      batch_adjusted_best = .decode_fit_response(x, g, adjusted[adjusted_row]),
      batch_adjusted_source_batch = batch[adjusted_row],
      batch_adjusted_type = type[adjusted_row],
      stringsAsFactors = FALSE
    )

    observed_condition <- .decode(x@design, g$X[observed_row, , drop = FALSE])
    adjusted_condition <- .decode(x@design, g$X[adjusted_row, , drop = FALSE])
    for (factor in x@design@label) {
      output[[paste0("observed_", factor)]] <- observed_condition[[factor]]
      output[[paste0("batch_adjusted_", factor)]] <- adjusted_condition[[factor]]
    }

    output
  })

  table <- do.call(rbind, rows)
  rownames(table) <- NULL
  settings <- list(include_references = include_references)

  x@analyses[[name]] <- list(
    table = table,
    surrogate = surrogate_name,
    settings = settings,
    provenance = .provenance_analysis(
      x = x,
      surrogate_name = surrogate_name,
      kind = "cumulative_best",
      settings = settings
    ),
    call = match.call()
  )

  methods::validObject(x)
  x
}


#' @export
cumulativeBest <- function(x, name = NULL, strict = FALSE, details = FALSE) {
  entry <- .analysis_entry(x = x, name = name, kind = "cumulative_best", strict = strict)
  if (.validate_flag(details, "details")) entry else entry$table
}


#' @export
cumulativeBestNames <- function(x) {
  .analysis_names(x, "cumulative_best")
}


#================== Cumulative-Best Plot ==================

#plotCumulativeBest() draws a stored trajectory against the completed batches,
#and returns the plotted table invisibly.

#' @export
plotCumulativeBest <- function(
    x,
    name = NULL,
    strict = FALSE,
    curves = c("raw", "batch_adjusted"),
    legend = TRUE,
    las = 1L,
    main = NULL,
    ylab = NULL,
    colors = c(raw = "#B8C4CE", batch_adjusted = "#6F7B86"),
    lty = c(raw = 3L, batch_adjusted = 2L),
    lwd = 2.5,
    pch = c(raw = 16L, batch_adjusted = 17L)
) {
  legend <- .validate_flag(legend, "legend")
  las <- .validate_las(las)

  result <- cumulativeBest(x, name = name, strict = strict)

  column_by_curve <- c(
    raw = "observed_best",
    batch_adjusted = "batch_adjusted_best"
  )
  label_by_curve <- c(
    raw = "Raw observed best",
    batch_adjusted = "Batch adjusted"
  )

  if (
    !is.character(curves) ||
    !length(curves) ||
    anyNA(curves) ||
    any(!curves %in% names(column_by_curve))
  ) {
    stop(
      '"curves" must contain one or more of: ',
      paste(names(column_by_curve), collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  curves <- unique(curves)

  colors <- colors[curves]
  lty <- lty[curves]
  pch <- pch[curves]

  if (
    length(colors) != length(curves) ||
    anyNA(colors) ||
    length(lty) != length(curves) ||
    anyNA(lty) ||
    length(pch) != length(curves) ||
    anyNA(pch)
  ) {
    stop(
      '"colors", "lty", and "pch" must provide values for every requested curve.',
      call. = FALSE
    )
  }

  if (is.null(ylab)) {
    ylab <- paste0("Best ", x@response, " so far")
  }

  plot_margins <- graphics::par("mai")
  if (is.null(main)) {
    plot_margins[3L] <- min(plot_margins[3L], 0.20)
  }
  if (las %in% c(2L, 3L)) {
    axis_cex <- 0.8
    label_extent <- max(
      graphics::strwidth(
        as.character(result$batch),
        units = "inches",
        cex = axis_cex
      ),
      na.rm = TRUE
    )
    plot_margins[1L] <- max(plot_margins[1L], label_extent + 0.65)
  }
  #Restore only the margins changed here, so multi-panel layouts keep working.
  old_margins <- graphics::par(mai = plot_margins)
  on.exit(graphics::par(old_margins), add = TRUE)

  trajectory <- do.call(cbind, lapply(
    unname(column_by_curve[curves]),
    function(column) result[[column]]
  ))
  colnames(trajectory) <- curves

  limits <- .padded_range(trajectory)

  graphics::plot(
    x = result$batch_index,
    y = trajectory[, 1L],
    type = "n",
    xaxt = "n",
    xlab = "Completed batch",
    ylab = ylab,
    main = main,
    ylim = limits,
    bty = "l"
  )

  graphics::axis(
    side = 1L,
    at = result$batch_index,
    labels = result$batch,
    las = las,
    cex.axis = 0.8
  )
  graphics::grid(nx = NA, ny = NULL, col = "#E6E9EC", lty = 1L)

  graphics::matlines(
    x = result$batch_index,
    y = trajectory,
    type = "o",
    lty = unname(lty),
    lwd = lwd,
    pch = unname(pch),
    col = unname(colors)
  )

  if (legend) {
    graphics::legend(
      "bottomright",
      legend = unname(label_by_curve[curves]),
      col = unname(colors),
      lty = unname(lty),
      lwd = lwd,
      pch = unname(pch),
      bty = "n",
      cex = 0.85
    )
  }

  invisible(result)
}
