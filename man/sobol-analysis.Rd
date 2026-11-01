\name{sobol-analysis}
\alias{sobol-analysis}
\alias{runSobol}
\alias{sobolIndices}
\alias{sobolInteractions}
\alias{sobolNames}
\alias{plotSobol}
\alias{plotSobolInteractions}
\title{Sobol' sensitivity analysis of the fitted surface}
\description{
Estimate how much of the variation of the fitted response surface each factor
explains, alone and through interactions. \code{runSobol()} computes first-,
total- and optionally second-order Sobol' indices of a stored surrogate and
stores them in the campaign; the accessors retrieve them and the plot
functions draw them.
}
\usage{
runSobol(x, name = "sobol", surrogate_name = NULL, N = 4096L, bootstrap = 0L,
         conf = 0.95, second_order = FALSE, scale = c("model", "response"),
         seed = NULL, overwrite = FALSE, first = "saltelli", total = "jansen",
         interval_type = "percent", matrices = c("A", "B", "AB"))
sobolIndices(x, name = NULL, strict = FALSE, details = FALSE)
sobolInteractions(x, name = NULL, strict = FALSE, details = FALSE)
sobolNames(x)
plotSobol(x, name = NULL, strict = FALSE, stacked = TRUE,
          colors = c(first_order = "#E76F61", total_order = "#00AEB5"),
          border = "#333333", show_intervals = TRUE, legend = TRUE, las = 2L,
          ...)
plotSobolInteractions(x, name = NULL, strict = FALSE, dots = FALSE,
                      color = "#00AEB5", interval_color = "#333333",
                      pch = 19L, show_intervals = TRUE, las = 2L, ...)
}
\arguments{
  \item{x}{A \code{\linkS4class{BOCampaign}}.}
  \item{name}{Name under which to store the analysis (default
    \code{"sobol"}). For the accessors and plots, the name to retrieve;
    \code{NULL} (default) gives the most recently stored current analysis.}
  \item{surrogate_name}{Name of the stored surrogate to analyse. \code{NULL}
    (default) uses the most recently stored current surrogate.}
  \item{N}{Base sample size of the integration design (default 4096). The
    surface is evaluated \eqn{N(D + 2)} times for first- and total-order
    indices, and more with \code{second_order = TRUE}.}
  \item{bootstrap}{Number of bootstrap replicates for confidence intervals;
    0 (default) gives point estimates only.}
  \item{conf}{Confidence level of the bootstrap intervals (default 0.95).}
  \item{second_order}{Logical; also estimate the second-order index of every
    pair of factors. Default \code{FALSE}.}
  \item{scale}{\code{"model"} (default) to analyse the surface on the scale it
    was fitted (after any log transform), or \code{"response"} to analyse it in
    response units.}
  \item{seed}{\code{NULL} (default) or an integer seed for the integration
    design and the bootstrap. The global random-number stream is left
    untouched.}
  \item{overwrite}{Logical; allow replacing a stored analysis with the same
    name. Default \code{FALSE}.}
  \item{first, total}{Estimators of the first- and total-order indices, passed
    to \code{\link[sensobol]{sobol_indices}}. Defaults \code{"saltelli"} and
    \code{"jansen"}.}
  \item{interval_type}{Type of bootstrap interval, passed to
    \code{\link[sensobol]{sobol_indices}} (default \code{"percent"}).}
  \item{matrices}{Layout of the integration design, passed to
    \code{\link[sensobol]{sobol_matrices}} (default
    \code{c("A", "B", "AB")}); it must suit the chosen estimators.}
  \item{strict}{Logical; stop if the requested analysis is historical (see
    \code{\link{summary.BOCampaign}}).}
  \item{details}{Logical; return the complete stored entry instead of the
    table.}
  \item{stacked}{Logical; in \code{plotSobol()}, draw one bar per factor
    with the first-order index at the bottom and the interaction part
    \eqn{S_{T_i} - S_i}{S_Ti - S_i} on top, so the full bar is the total-order index
    (default). \code{FALSE} draws the two indices side by side.}
  \item{colors}{Colours of the first- and total-order bars, named
    \code{first_order} and \code{total_order}.}
  \item{border}{Colour of bar borders and intervals.}
  \item{show_intervals}{Logical; draw bootstrap intervals when they were
    computed.}
  \item{legend}{Logical; draw a legend.}
  \item{las}{Orientation of the axis labels, as in
    \code{\link[graphics]{par}}.}
  \item{dots}{Logical; in \code{plotSobolInteractions()}, draw points instead
    of bars. Default \code{FALSE}.}
  \item{color}{Colour of the second-order bars or points.}
  \item{interval_color}{Colour of bar borders and intervals in
    \code{plotSobolInteractions()}.}
  \item{pch}{Point symbol when \code{dots = TRUE}.}
  \item{...}{Further arguments passed to \code{\link[graphics]{plot}}, such
    as \code{main}, \code{xlab}, \code{ylab}, \code{xlim} or \code{ylim}.}
}
\details{
\strong{Indices.} The analysis concerns the posterior mean of the batch-free
surface \eqn{\beta_0 + f(x)}{beta0 + f(x)} (see \code{\link{runSurrogate}}), with every factor
varying independently and uniformly over its range on the coded scale (so
log-scaled factors vary uniformly in \eqn{\log_{10}}{log10}). The first-order index
\eqn{S_i}{S_i} is the share of the surface's variance explained by factor \eqn{i}
alone; the total-order index \eqn{S_{T_i}}{S_Ti} is the share in which factor
\eqn{i} takes part in any way, alone or in interaction with others. Their
difference \eqn{S_{T_i} - S_i}{S_Ti - S_i} measures how much factor \eqn{i} acts through
interactions, and a factor with \eqn{S_{T_i}}{S_Ti} near zero could be fixed
anywhere in its range with little effect on the response. The second-order
index \eqn{S_{ij}}{S_ij} is the share due to the interaction of factors \eqn{i} and
\eqn{j} alone.

\strong{Computation.} The integration design (Sobol' points with a random
Cranley--Patterson shift) and the estimators come from \pkg{sensobol}. The
default estimators are those of Saltelli et al. (2010) for \eqn{S_i}{S_i} and
Jansen (1999) for \eqn{S_{T_i}}{S_Ti}. Estimates can be slightly negative or exceed
one through numerical error, particularly for small indices and small
\code{N}. Bootstrap intervals describe this numerical error only, not the
uncertainty of the fitted surrogate itself, which is usually much larger.

With \code{scale = "response"}, the posterior mean is decoded to response
units before the indices are computed. This gives the same indices for an
untransformed response, but can differ for a log response, where effects that
add on the log scale multiply on the response scale.

\strong{Plots.} \code{plotSobol()} orders the factors by decreasing
total-order index, and \code{plotSobolInteractions()} draws the pairs as
horizontal bars with the largest index at the top.

Analyses share their names with the other stored analyses (see
\code{\link{recommend}}). An analysis is current as long as its surrogate is.
}
\value{
\code{runSobol()} returns the updated campaign.
\code{sobolIndices()} returns a data frame with one row per factor and
columns \code{factor}, \code{first_order}, \code{first_lower},
\code{first_upper}, \code{total_order}, \code{total_lower} and
\code{total_upper} (interval limits are \code{NA} without bootstrap).
\code{sobolInteractions()} returns a data frame with one row per pair and
columns \code{factor_1}, \code{factor_2}, \code{second_order},
\code{second_lower} and \code{second_upper}; it stops if the analysis has no
second-order indices. With \code{details = TRUE}, both return the stored
entry (\code{table}, \code{interactions}, \code{surrogate}, \code{settings},
\code{provenance} and \code{call}). \code{sobolNames()} returns the names of
all stored Sobol' analyses. The plot functions draw the plot and return the
plotted table invisibly.
}
\references{
Jansen, M. J. W. (1999). Analysis of variance designs for model output.
\emph{Computer Physics Communications}, \bold{117}, 35--43.

Puy, A., Lo Piano, S., Saltelli, A. and Levin, S. A. (2022). sensobol: An R
package to compute variance-based sensitivity indices. \emph{Journal of
Statistical Software}, \bold{102}(5), 1--37. \doi{10.18637/jss.v102.i05}

Saltelli, A., Annoni, P., Azzini, I., Campolongo, F., Ratto, M. and Tarantola,
S. (2010). Variance based sensitivity analysis of model output. Design and
estimator for the total sensitivity index. \emph{Computer Physics
Communications}, \bold{181}(2), 259--270.
}
\seealso{
\code{\link{runMainEffects}} for the shape of each factor's effect,
\code{\link{runResponseSurfaces}} for two-factor views.
}
\examples{
surface <- demoSurface()
study <- campaign(surface$design, response = "Ion_response")
study <- proposeSpaceFilling(study, batch = paste0("B", 1:3), q = 10, seed = 1)
for (b in paste0("B", 1:3)) {
  sheet <- proposal(study, b)
  sheet$Ion_response <- simulateBatch(surface, sheet)
  study <- addBatch(study, sheet)
}
study <- runSurrogate(study, "B3_fit", seed = 1)

\donttest{
study <- runSobol(study, N = 1024, second_order = TRUE, seed = 1)
sobolIndices(study)
plotSobol(study)
plotSobol(study, stacked = FALSE, main = "First and total order")

head(sobolInteractions(study))
plotSobolInteractions(study)
plotSobolInteractions(study, dots = TRUE)

## Bootstrap intervals for the numerical error
study <- runSobol(study, name = "sobol_boot", N = 1024, bootstrap = 100,
                  seed = 1)
plotSobol(study, "sobol_boot")
}
}
