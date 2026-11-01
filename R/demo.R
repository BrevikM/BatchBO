################################################################################
# Demo.R -- self-contained synthetic laboratory for package examples
################################################################################


#================== Internal Hartmann-6 Surface ==================

.demo_hartmann6_spec <- function(
    interpretation = c("mass_spectrometry", "unit_cube")
) {
  interpretation <- match.arg(interpretation)

  if (identical(interpretation, "mass_spectrometry")) {
    return(list(
      interpretation = interpretation,
      labels = c(
        "Capillary_voltage",
        "Nebulizer_pressure",
        "Nozzle_voltage",
        "Drying_gas_flow",
        "Drying_gas_temperature",
        "Sheath_gas_flow"
      ),
      lower = c(2500, 20, 0, 8, 200, 8),
      upper = c(4500, 50, 2000, 12, 350, 12),
      units = c("V", "psi", "V", "L/min", "degC", "L/min"),
      response = "Ion_response"
    ))
  }

  list(
    interpretation = interpretation,
    labels = paste0("x", seq_len(6L)),
    lower = rep(0, 6L),
    upper = rep(1, 6L),
    units = rep("", 6L),
    response = "Response"
  )
}


.demo_hartmann6_to_unit <- function(data, spec) {
  physical <- as.matrix(data[, spec$labels, drop = FALSE])
  coded <- sweep(physical, 2L, spec$lower, "-")
  coded <- sweep(coded, 2L, spec$upper - spec$lower, "/")
  colnames(coded) <- spec$labels
  coded
}


.demo_hartmann6_from_unit <- function(X, spec) {
  X <- as.matrix(X)
  physical <- sweep(X, 2L, spec$upper - spec$lower, "*")
  physical <- sweep(physical, 2L, spec$lower, "+")
  colnames(physical) <- spec$labels
  as.data.frame(physical, stringsAsFactors = FALSE)
}


.demo_hartmann6_reward <- function(X) {
  X <- as.matrix(X)
  alpha <- c(1.0, 1.2, 3.0, 3.2)

  A <- matrix(
    c(
      10, 3, 17, 3.5, 1.7, 8,
      0.05, 10, 17, 0.1, 8, 14,
      3, 3.5, 1.7, 10, 17, 8,
      17, 8, 0.05, 10, 0.1, 14
    ),
    nrow = 4L,
    byrow = TRUE
  )

  P <- 1e-4 * matrix(
    c(
      1312, 1696, 5569, 124, 8283, 5886,
      2329, 4135, 8307, 3736, 1004, 9991,
      2348, 1451, 3522, 2883, 3047, 6650,
      4047, 8828, 8732, 5743, 1091, 381
    ),
    nrow = 4L,
    byrow = TRUE
  )

  components <- vapply(seq_len(4L), function(i) {
    centre <- matrix(P[i, ], nrow(X), ncol(X), byrow = TRUE)
    weights <- matrix(A[i, ], nrow(X), ncol(X), byrow = TRUE)
    alpha[i] * exp(-rowSums(weights * (X - centre)^2))
  }, numeric(nrow(X)))

  rowSums(matrix(components, nrow = nrow(X)))
}


.demo_hartmann6_response <- function(data, spec) {
  1000 + 500 * .demo_hartmann6_reward(
    .demo_hartmann6_to_unit(data, spec)
  )
}


#The demonstration functions take the surface returned by demoSurface().
.check_demo_surface <- function(x) {
  if (!inherits(x, "DemoSurface")) {
    stop('"x" must be a demoSurface().', call. = FALSE)
  }
  invisible(x)
}


#The response column a simulated batch is expected to populate. An explicit
#name always wins, so a campaign may name its response whatever it likes.
.demo_response_name <- function(x, response = NULL) {
  if (is.null(response)) x$response else .validate_name(response, "response")
}


.validate_demo_sd <- function(value, name) {
  if (
    !is.numeric(value) ||
    is.logical(value) ||
    length(value) != 1L ||
    is.na(value) ||
    !is.finite(value) ||
    value < 0
  ) {
    stop('"', name, '" must be one finite non-negative number.', call. = FALSE)
  }

  as.numeric(value)
}


#================== Public Demonstration Interface ==================

#demoSurface() returns the synthetic laboratory on its own: the Hartmann-6
#surface together with the between-batch and residual standard deviations used
#to contaminate it. Pass it to simulateBatch(), demoResponse() and demoTruth()
#to generate measurements for a campaign declared with design_space() and
#campaign(). The interpretation changes only factor names and physical ranges.
#' @export
demoSurface <- function(
    batch_sd = 80,
    residual_sd = 40,
    interpretation = c("mass_spectrometry", "unit_cube")
) {
  batch_sd <- .validate_demo_sd(batch_sd, "batch_sd")
  residual_sd <- .validate_demo_sd(residual_sd, "residual_sd")
  spec <- .demo_hartmann6_spec(interpretation)

  factors <- lapply(seq_along(spec$labels), function(i) {
    def_factor(spec$labels[i], spec$lower[i], spec$upper[i])
  })

  structure(
    list(
      surface = "hartmann6",
      interpretation = spec$interpretation,
      units = stats::setNames(spec$units, spec$labels),
      batch_sd = batch_sd,
      residual_sd = residual_sd,
      design = do.call(design_space, factors),
      response = spec$response
    ),
    class = "DemoSurface"
  )
}


#' @method print DemoSurface
#' @export
print.DemoSurface <- function(x, ...) {
  cat("Hidden demonstration response surface (Hartmann-6)\n")
  cat("  Interpretation:  ", x$interpretation, "\n", sep = "")
  cat("  Batch SD:        ", x$batch_sd, "\n", sep = "")
  cat("  Residual SD:     ", x$residual_sd, "\n", sep = "")
  cat("  Default response:", x$response, "\n", sep = " ")
  cat("  Factors:\n")
  print(
    data.frame(
      factor = x$design@label,
      lower = x$design@lower,
      upper = x$design@upper,
      units = unname(x$units),
      check.names = FALSE
    ),
    row.names = FALSE
  )
  invisible(x)
}


#simulateBatch() simulates one completed laboratory batch: one batch offset per
#call plus independent residual error, in the row order of the run sheet.
#' @export
simulateBatch <- function(x, run_sheet, seed = NULL, response = NULL) {
  .check_demo_surface(x)
  spec <- .demo_hartmann6_spec(x$interpretation)
  design <- x$design
  response <- .demo_response_name(x, response)
  seed <- .validate_seed(seed)
  run_sheet <- as.data.frame(run_sheet, stringsAsFactors = FALSE)

  if (!nrow(run_sheet)) {
    stop('"run_sheet" must contain at least one condition.', call. = FALSE)
  }

  run_sheet <- .validate_numeric_columns(
    run_sheet,
    design@label,
    '"run_sheet"'
  )
  .validate_log_domain(run_sheet, design, '"run_sheet"')

  if (
    response %in% names(run_sheet) &&
    any(!is.na(run_sheet[[response]]))
  ) {
    stop(
      '"run_sheet" already contains one or more response values.',
      call. = FALSE
    )
  }

  simulated <- .with_seed(seed, {
    batch_offset <- stats::rnorm(1L, sd = x$batch_sd)
    .demo_hartmann6_response(run_sheet, spec) +
      batch_offset +
      stats::rnorm(nrow(run_sheet), sd = x$residual_sd)
  })

  as.numeric(simulated)
}


#demoResponse() evaluates the hidden response without simulated noise. Intended
#only for checking the synthetic vignette; a real campaign has no such oracle.
#' @export
demoResponse <- function(x, newdata) {
  .check_demo_surface(x)
  spec <- .demo_hartmann6_spec(x$interpretation)
  design <- x$design
  newdata <- as.data.frame(newdata, stringsAsFactors = FALSE)
  newdata <- .validate_numeric_columns(
    newdata,
    design@label,
    '"newdata"'
  )
  .validate_log_domain(newdata, design, '"newdata"')
  .demo_hartmann6_response(newdata, spec)
}


#demoTruth() returns the known optimum and its batch-free response.
#' @export
demoTruth <- function(x, response = NULL) {
  .check_demo_surface(x)
  spec <- .demo_hartmann6_spec(x$interpretation)
  design <- x$design
  response <- .demo_response_name(x, response)

  coded <- matrix(
    c(0.20169, 0.150011, 0.476874, 0.275332, 0.311652, 0.6573),
    nrow = 1L,
    dimnames = list(NULL, design@label)
  )

  condition <- .demo_hartmann6_from_unit(coded, spec)
  condition[[response]] <- .demo_hartmann6_response(condition, spec)
  condition
}


#demoCompare() compares one condition, such as a recommendation, with the known
#optimum: both conditions, their batch-free responses, and the regret.
#' @export
demoCompare <- function(x, condition, response = NULL) {
  .check_demo_surface(x)
  response <- .demo_response_name(x, response)
  condition <- as.data.frame(condition, stringsAsFactors = FALSE)

  if (nrow(condition) != 1L) {
    stop('"condition" must contain exactly one condition.', call. = FALSE)
  }

  truth <- demoTruth(x, response = response)
  labels <- x$design@label
  value <- demoResponse(x, condition)

  data.frame(
    condition = c(unlist(condition[1L, labels]), value, truth[[response]] - value),
    optimum = c(unlist(truth[1L, labels]), truth[[response]], 0),
    row.names = c(labels, response, "regret")
  )
}
