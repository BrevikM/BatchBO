################################################################################
# Recommendation.R -- optimize and store the fitted latent response surface
################################################################################


#================== Posterior-Mean Optimization ==================

#Numerically maximize the batch-free posterior mean over the coded unit cube.
.optimize_posterior_mean <- function(
    g,
    n_screen = 20000L,
    restarts = 32L,
    maxit = 500L,
    seed = NULL
) {
  D <- ncol(g$X)
  controls <- .validate_acquisition_controls(
    n_screen = n_screen,
    restarts = restarts,
    maxit = maxit,
    lower = rep(0, D),
    upper = rep(1, D),
    D = D,
    max_screen_factor = 1L,
    seed = seed
  )

  objective <- function(X) {
    gp_predict(g, as.matrix(X))$mean
  }

  screen <- .sbl(
    n = controls$n_screen,
    D = D,
    lower = controls$lower,
    upper = controls$upper,
    seed = controls$seed
  )

  .optimize_acquisition(
    objective = objective,
    screen = screen,
    lower = controls$lower,
    upper = controls$upper,
    restarts = controls$restarts,
    maxit = controls$maxit,
    max_screen_factor = 1L,
    seed = controls$seed
  )
}


#================== Stored Recommendation ==================

#recommend() numerically optimizes the stored GP's batch-free posterior mean over
#the declared design space. Unlike an acquisition function, it contains no
#exploration term: it answers which condition the model currently believes is
#best. A Sobol screen is followed by bounded L-BFGS-B refinement from the best
#screened starts. recommendation() retrieves a stored result.

#' @export
recommend <- function(
    x,
    name = "recommendation",
    surrogate_name = NULL,
    n_screen = 20000L,
    restarts = 32L,
    maxit = 500L,
    seed = NULL,
    level = 0.95,
    overwrite = FALSE
) {
  if (!methods::is(x, "BOCampaign")) {
    stop('"x" must be a BOCampaign.', call. = FALSE)
  }

  name <- .validate_name(name)
  overwrite <- .validate_overwrite_name(
    name,
    names(x@analyses),
    overwrite,
    what = "analysis"
  )
  surrogate_name <- .resolve_surrogate_name(x, surrogate_name)
  seed <- .validate_seed(seed)

  fit <- surrogate(x, surrogate_name, strict = FALSE)

  optimization <- .with_seed(
    seed,
    .optimize_posterior_mean(
      g = fit,
      n_screen = n_screen,
      restarts = restarts,
      maxit = maxit,
      seed = seed
    )
  )

  coded_optimum <- matrix(
    optimization$x,
    nrow = 1L,
    dimnames = list(NULL, x@design@label)
  )
  physical_optimum <- .decode(x@design, coded_optimum)
  predicted_optimum <- predictSurrogate(
    x,
    physical_optimum,
    surrogate_name = surrogate_name,
    level = level,
    strict = FALSE
  )

  settings <- list(
    n_screen = as.integer(n_screen),
    restarts = as.integer(restarts),
    maxit = as.integer(maxit),
    seed = seed,
    level = level
  )

  result <- list(
    condition = predicted_optimum,
    coded = coded_optimum,
    surrogate = surrogate_name,
    settings = settings,
    diagnostics = list(
      standardized_posterior_mean = optimization$value,
      convergence = optimization$convergence,
      message = .optimizer_message(
        optimization$convergence,
        optimization$message
      ),
      screen_max = optimization$screen_max,
      screen_attempts = optimization$screen_attempts,
      screen_size = optimization$screen_size
    ),
    provenance = .provenance_analysis(
      x = x,
      surrogate_name = surrogate_name,
      kind = "recommendation",
      settings = settings
    ),
    call = match.call()
  )

  x@analyses[[name]] <- result
  methods::validObject(x)
  x
}


#' @export
recommendation <- function(
    x,
    name = NULL,
    strict = FALSE,
    details = FALSE
) {
  entry <- .analysis_entry(
    x = x,
    name = name,
    kind = "recommendation",
    strict = strict
  )

  details <- .validate_flag(details, "details")
  if (details) entry else entry$condition
}


#' @export
recommendationNames <- function(x) {
  .analysis_names(x, "recommendation")
}
