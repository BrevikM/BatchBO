\name{recommendation}
\alias{recommendation}
\alias{recommend}
\alias{recommendationNames}
\title{Recommend the best condition}
\description{
\code{recommend()} finds the condition within the design space where the
fitted surrogate predicts the best response, and stores it in the campaign.
\code{recommendation()} retrieves it.
}
\usage{
recommend(x, name = "recommendation", surrogate_name = NULL,
          n_screen = 20000L, restarts = 32L, maxit = 500L, seed = NULL,
          level = 0.95, overwrite = FALSE)
recommendation(x, name = NULL, strict = FALSE, details = FALSE)
recommendationNames(x)
}
\arguments{
  \item{x}{A \code{\linkS4class{BOCampaign}}.}
  \item{name}{Name under which to store the recommendation (default
    \code{"recommendation"}). For \code{recommendation()}, the name to
    retrieve; \code{NULL} (default) gives the most recently stored current
    recommendation.}
  \item{surrogate_name}{Name of the stored surrogate to use. \code{NULL}
    (default) uses the most recently stored current surrogate.}
  \item{n_screen}{Number of points in the Sobol' screen (default 20000).}
  \item{restarts}{Number of best screened points refined with L-BFGS-B
    (default 32).}
  \item{maxit}{Maximum number of L-BFGS-B iterations per start (default 500).}
  \item{seed}{\code{NULL} (default) or an integer seed for the Sobol' screen.
    The global random-number stream is left untouched.}
  \item{level}{Probability of the central credible interval (default 0.95).}
  \item{overwrite}{Logical; allow replacing a stored analysis with the same
    name. Default \code{FALSE}.}
  \item{strict}{Logical; stop if the requested recommendation is historical
    (see \code{\link{summary.BOCampaign}}).}
  \item{details}{Logical; return the complete stored entry instead of the
    recommended condition.}
}
\details{
The recommendation maximizes (or, for a minimized response, minimizes) the
posterior mean of the batch-free surface \eqn{\beta_0 + f(x)}{beta0 + f(x)} over the
declared bounds. Unlike the acquisition functions of
\code{\link{proposeBatch}}, it rewards no exploration: it is the condition the
model currently believes is best, and is meant as the final answer once the
experimental budget is spent. Its prediction refers to a typical batch (see
\code{\link{predictSurrogate}}).

The posterior mean is evaluated on a randomized Sobol' set of
\code{n_screen} points and refined with bounded L-BFGS-B from the best
\code{restarts} of them.

Recommendations share their names with the other stored analyses
(\code{\link{runSobol}}, \code{\link{runMainEffects}},
\code{\link{runResponseSurfaces}}, \code{\link{runCumulativeBest}}), so each
name can be used only once among them. A recommendation is current as long as
its surrogate is.
}
\value{
\code{recommend()} returns the updated campaign. \code{recommendation()}
returns a one-row data frame with the recommended condition in physical units
and the prediction columns of \code{\link{predictSurrogate}}; with
\code{details = TRUE}, the stored entry: the \code{condition}, the condition
on the coded scale (\code{coded}), the \code{surrogate} name,
\code{settings}, optimizer \code{diagnostics} (standardized posterior mean,
convergence code and message, and screen maximum, attempts and size),
\code{provenance} and the \code{call}. \code{recommendationNames()} returns
the names of all stored recommendations.
}
\seealso{
\code{\link{predictSurrogate}}, \code{\link{bestObserved}} for the best
measured run, \code{\link{plotResponseSurfaces}} to show the recommendation
on the fitted surface.
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

study <- recommend(study, seed = 1)
best <- recommendation(study)
best

## Compare with the true surface and its optimum
demoResponse(surface, best)
demoTruth(surface)
recommendation(study, details = TRUE)$diagnostics
}
