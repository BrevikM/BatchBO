\name{runCumulativeBest}
\alias{runCumulativeBest}
\alias{cumulativeBest}
\alias{cumulativeBestNames}
\alias{plotCumulativeBest}
\title{Best response found as the campaign progressed}
\description{
Follow the best tested condition as the batches accumulate, both as measured
and after removing each batch's estimated offset. \code{runCumulativeBest()}
computes and stores the trajectory, \code{cumulativeBest()} retrieves it and
\code{plotCumulativeBest()} draws it.
}
\usage{
runCumulativeBest(x, name = "cumulative_best", surrogate_name = NULL,
                  include_references = TRUE, overwrite = FALSE)
cumulativeBest(x, name = NULL, strict = FALSE, details = FALSE)
cumulativeBestNames(x)
plotCumulativeBest(x, name = NULL, strict = FALSE,
                   curves = c("raw", "batch_adjusted"), legend = TRUE,
                   las = 1L, main = NULL, ylab = NULL,
                   colors = c(raw = "#B8C4CE", batch_adjusted = "#6F7B86"),
                   lty = c(raw = 3L, batch_adjusted = 2L), lwd = 2.5,
                   pch = c(raw = 16L, batch_adjusted = 17L))
}
\arguments{
  \item{x}{A \code{\linkS4class{BOCampaign}}.}
  \item{name}{Name under which to store the trajectory (default
    \code{"cumulative_best"}). For the accessor and plot, the name to
    retrieve; \code{NULL} (default) gives the most recently stored current
    trajectory.}
  \item{surrogate_name}{Name of the stored surrogate whose data and batch
    effects are used. \code{NULL} (default) uses the most recently stored
    current surrogate.}
  \item{include_references}{Logical; let reference runs inside the bounds
    count as tested conditions (default \code{TRUE}). Sample runs always
    count; references outside the bounds never do.}
  \item{overwrite}{Logical; allow replacing a stored analysis with the same
    name. Default \code{FALSE}.}
  \item{strict}{Logical; stop if the requested trajectory is historical (see
    \code{\link{summary.BOCampaign}}).}
  \item{details}{Logical; return the complete stored entry instead of the
    table.}
  \item{curves}{Which curves to draw: one or both of \code{"raw"} and
    \code{"batch_adjusted"}.}
  \item{legend}{Logical; draw a legend.}
  \item{las}{Orientation of the batch labels on the horizontal axis, as in
    \code{\link[graphics]{par}}; use 2 for long labels.}
  \item{main}{Plot title; \code{NULL} for none.}
  \item{ylab}{Vertical-axis label; \code{NULL} gives
    \code{"Best <response> so far"}.}
  \item{colors, lty, pch}{Colour, line type and point symbol of each curve,
    as vectors named \code{raw} and \code{batch_adjusted}.}
  \item{lwd}{Line width.}
}
\details{
At stage \eqn{i}, the runs of the first \eqn{i} batches (in the order they were
added) are considered, and the best one is reported in two ways:
\describe{
  \item{Raw}{The best measured response so far.}
  \item{Batch adjusted}{The best response after subtracting each run's
    estimated batch offset, the posterior mean (best linear unbiased
    predictor) \eqn{E[u_b \mid y] = \sigma_b^2 \sum_{i \in b} a_i}{E[u_b | y] = sigma_b^2 sum_(i in b) a_i},
    with \eqn{a = V^{-1}(y - \beta_0)}{a = V^-1 (y - beta0)}, from the
    surrogate (see \code{\link{runSurrogate}}).}
}
A batch that happens to read high can make the raw curve jump without any
real improvement; the batch-adjusted curve puts all batches on a common level,
so its rise reflects better conditions. Both curves are maxima of noisy values
and so tend to overstate the true best slightly; the batch-adjusted curve of a
batch-naive surrogate is essentially the raw curve.

The batch offsets of the chosen surrogate are used at every stage, so the
trajectory covers the batches that surrogate was fitted to, and it is current
exactly while that surrogate is.

For a minimized response, \emph{best} means smallest. Values are reported in
response units.
}
\value{
\code{runCumulativeBest()} returns the updated campaign.
\code{cumulativeBest()} returns a data frame with one row per batch and
columns \code{batch}, \code{batch_index}, \code{completed_runs},
\code{eligible_runs} (runs counted so far), \code{observed_best},
\code{observed_best_source_batch}, \code{observed_best_type},
\code{batch_adjusted_best}, \code{batch_adjusted_source_batch},
\code{batch_adjusted_type}, and the conditions of both best runs
(\code{observed_<factor>} and \code{batch_adjusted_<factor>}); with
\code{details = TRUE}, the stored entry (\code{table}, \code{surrogate},
\code{settings}, \code{provenance} and \code{call}).
\code{cumulativeBestNames()} returns the names of all stored trajectories.
\code{plotCumulativeBest()} draws the plot and returns the table invisibly.
}
\seealso{
\code{\link{bestObserved}} for the best measured run,
\code{\link{recommend}} for the best predicted condition.
}
\examples{
surface <- demoSurface()
study <- campaign(surface$design, response = "Ion_response")
study <- proposeSpaceFilling(study, batch = paste0("B", 1:4), q = 6, seed = 1)
for (b in paste0("B", 1:4)) {
  sheet <- proposal(study, b)
  sheet$Ion_response <- simulateBatch(surface, sheet)
  study <- addBatch(study, sheet)
}
study <- runSurrogate(study, "B4_fit", seed = 1)

study <- runCumulativeBest(study)
cumulativeBest(study)[, c("batch", "observed_best", "batch_adjusted_best")]
plotCumulativeBest(study)
plotCumulativeBest(study, curves = "batch_adjusted", legend = FALSE)
}
