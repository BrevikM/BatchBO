\name{summary.GP}
\alias{summary.GP}
\alias{print.GP}
\alias{print.summary.GP}
\title{Summarize a fitted surrogate}
\description{
Print a short overview of a fitted surrogate, or summarize its variance
components, length scales and optimization diagnostics.
}
\usage{
\method{print}{GP}(x, ...)
\method{summary}{GP}(object, scaled = TRUE, ...)
\method{print}{summary.GP}(x, digits = 4, ...)
}
\arguments{
  \item{x}{A fitted surrogate (class \code{"GP"}) from
    \code{\link{surrogate}}, or a \code{summary.GP} object.}
  \item{object}{A fitted surrogate.}
  \item{scaled}{Logical; report variances on the standardized response scale
    (default \code{TRUE}), or in response units (squared), on the log scale for
    a log response.}
  \item{digits}{Number of significant digits printed.}
  \item{...}{Ignored.}
}
\details{
\strong{Variance components.} The \emph{Response Surface} variance
\eqn{\sigma_f^2}{sigma_f^2} is the variation explained by the factors, the
\emph{Between-Batch} variance \eqn{\sigma_b^2}{sigma_b^2} the variation of the batch
offsets, and the \emph{Residual} variance \eqn{\sigma^2}{sigma^2} the run-to-run noise
within a batch (see \code{\link{runSurrogate}}). On the standardized scale the
Sample responses have unit variance, so the components can be compared with 1
and with each other. The Between-Batch row is \code{NA} for a batch-naive fit.

\strong{Length scales.} One per factor, on the coded scale where each factor
spans 1. A short length scale means the response changes quickly along that
factor; a length scale well above 1 means the response is nearly flat across
the whole range.

\code{at_bound} flags estimates that ended on an optimization bound, and the
optimization section reports how many starts converged. Per-start results are
in the \code{starts} element.
}
\value{
\code{summary()} returns an object of class \code{summary.GP}: a list with
\code{model} (model type, kernel, response, direction, numbers of runs,
batches and factors, negative log marginal likelihood and convergence),
\code{variance_components} (a data frame with \code{component},
\code{variance}, \code{sd} and \code{at_bound}), \code{length_scales} (a data
frame with \code{factor}, \code{length_scale} and \code{at_bound}),
\code{optimization} (numbers of starts and converged starts, and the optimizer
message), \code{starts} (the initial and final negative log likelihood,
convergence code and message of every start) and \code{scaled}. The print
methods return their argument invisibly.
}
\seealso{
\code{\link{runSurrogate}}, \code{\link{surrogate}}.
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

fit <- surrogate(study)
fit
summary(fit)
summary(fit, scaled = FALSE)$variance_components
}
