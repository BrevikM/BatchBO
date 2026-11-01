\name{predictSurrogate}
\alias{predictSurrogate}
\title{Predict the response surface at given conditions}
\description{
Predict the batch-free response, with a credible interval, at conditions given
in physical units, from a surrogate stored in the campaign.
}
\usage{
predictSurrogate(x, newdata, surrogate_name = NULL, level = 0.95,
                 strict = FALSE)
}
\arguments{
  \item{x}{A \code{\linkS4class{BOCampaign}}.}
  \item{newdata}{A data frame with one column per factor, in physical units.
    Other columns are ignored. Conditions outside the bounds are allowed (the
    prediction then extrapolates); log-scaled factors must be positive.}
  \item{surrogate_name}{Name of the stored surrogate to use. \code{NULL}
    (default) uses the most recently stored current surrogate.}
  \item{level}{Probability of the central credible interval (default 0.95).}
  \item{strict}{Logical; stop if the surrogate is historical (see
    \code{\link{summary.BOCampaign}}).}
}
\details{
The prediction is for \eqn{\beta_0 + f(x)}{beta0 + f(x)}, the response surface without batch
offset or residual noise (see \code{\link{runSurrogate}}): the response
expected for a typical batch. The interval therefore describes uncertainty
about the surface, not the spread of a future single measurement, which would
also include batch-to-batch and run-to-run variation.

The interval is computed on the model scale (after any log transform) as the
posterior mean \eqn{\pm}{+/-} a normal quantile times the posterior standard
deviation, and then transformed back. For an identity response,
\code{prediction} and \code{response_mean} are the same. For a log response,
\code{prediction} is \eqn{\exp(m)}{exp(m)}, the posterior median on the
response scale, while \code{response_mean} is
\eqn{\exp(m + s^2/2)}{exp(m + s^2/2)}, the posterior mean, where \eqn{m} and
\eqn{s} are the posterior mean and standard deviation on the log scale.
}
\value{
A data frame with the factor columns of \code{newdata} followed by
\item{prediction}{Central prediction in response units.}
\item{response_mean}{Posterior mean in response units.}
\item{prediction_lower, prediction_upper}{Limits of the credible interval in
  response units.}
\item{model_mean, model_sd}{Posterior mean and standard deviation on the model
  scale (the response, or its logarithm).}
The attributes \code{level}, \code{response_transform} and
\code{prediction_definition} record how the columns were obtained.
}
\seealso{
\code{\link{runSurrogate}}, \code{\link{recommend}} for the condition with the
best prediction.
}
\examples{
surface <- demoSurface()
study <- campaign(surface$design, response = "Ion_response")
study <- proposeSpaceFilling(study, batch = c("B1", "B2"), q = 8, seed = 1)
for (b in c("B1", "B2")) {
  sheet <- proposal(study, b)
  sheet$Ion_response <- simulateBatch(surface, sheet)
  study <- addBatch(study, sheet)
}
study <- runSurrogate(study, "B2_fit", seed = 1)

## Predict at the measured Sample conditions and compare with the truth
runs <- observations(study)
runs <- runs[runs$Type == "Sample", ]
pred <- predictSurrogate(study, runs)
plot(demoResponse(surface, runs), pred$prediction,
     xlab = "True surface", ylab = "Prediction")
abline(0, 1, lty = 2)
}
