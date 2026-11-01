################################################################################
# Prediction.R -- practitioner-facing predictions from a stored surrogate
################################################################################


#================== Surrogate Prediction ==================

#predictSurrogate() predicts the batch-free latent response at physical
#conditions. Its interval describes uncertainty in the fitted latent surface,
#not the variation of a future batch-affected measurement.

#' @export
predictSurrogate <- function(
    x,
    newdata,
    surrogate_name = NULL,
    level = 0.95,
    strict = FALSE
) {
  if (!methods::is(x, "BOCampaign")) {
    stop('"x" must be a BOCampaign.', call. = FALSE)
  }

  surrogate_name <- .resolve_surrogate_name(x, surrogate_name)
  strict <- .validate_flag(strict, "strict")

  if (
    !is.numeric(level) ||
    is.logical(level) ||
    length(level) != 1L ||
    is.na(level) ||
    !is.finite(level) ||
    level <= 0 ||
    level >= 1
  ) {
    stop('"level" must be a single number strictly between zero and one.', call. = FALSE)
  }

  newdata <- as.data.frame(newdata, stringsAsFactors = FALSE)
  if (!nrow(newdata)) {
    stop('"newdata" must contain at least one condition.', call. = FALSE)
  }

  newdata <- .validate_numeric_columns(
    newdata,
    x@design@label,
    '"newdata"'
  )
  .validate_log_domain(newdata, x@design, '"newdata"')

  fit <- surrogate(
    x,
    name = surrogate_name,
    strict = strict
  )

  coded_newdata <- .encode(x@design, newdata)
  prediction <- gp_predict(fit, coded_newdata)

  standardized_mean <- prediction$mean
  model_mean <- fit$y_mean +
    fit$sign * standardized_mean * fit$y_sd
  model_sd <- sqrt(pmax(prediction$var, 0)) * fit$y_sd

  critical_value <- stats::qnorm((1 + level) / 2)
  model_lower <- model_mean - critical_value * model_sd
  model_upper <- model_mean + critical_value * model_sd

  central_prediction <- switch(
    x@response_transform,
    identity = model_mean,
    log = exp(model_mean)
  )

  response_mean <- switch(
    x@response_transform,
    identity = model_mean,
    log = exp(model_mean + 0.5 * model_sd^2)
  )

  response_lower <- switch(
    x@response_transform,
    identity = model_lower,
    log = exp(model_lower)
  )

  response_upper <- switch(
    x@response_transform,
    identity = model_upper,
    log = exp(model_upper)
  )

  output <- cbind(
    newdata[, x@design@label, drop = FALSE],
    data.frame(
      prediction = central_prediction,
      response_mean = response_mean,
      prediction_lower = response_lower,
      prediction_upper = response_upper,
      model_mean = model_mean,
      model_sd = model_sd,
      stringsAsFactors = FALSE
    )
  )

  attr(output, "level") <- level
  attr(output, "response_transform") <- x@response_transform
  attr(output, "prediction_definition") <- if (
    identical(x@response_transform, "log")
  ) {
    "prediction is the median obtained by inverse-transforming the latent posterior mean"
  } else {
    "prediction is the latent posterior mean on the reported response scale"
  }

  output
}
