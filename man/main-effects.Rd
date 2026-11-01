\name{runMainEffects}
\alias{runMainEffects}
\alias{mainEffects}
\alias{mainEffectNames}
\alias{plotMainEffects}
\title{Main effects of the factors on the fitted surface}
\description{
Show how the fitted response changes along each factor, averaged over the
other factors. \code{runMainEffects()} computes the curves from a stored
surrogate and stores them in the campaign; \code{mainEffects()} retrieves them
and \code{plotMainEffects()} draws one panel per factor.
}
\usage{
runMainEffects(x, name = "main_effects", surrogate_name = NULL,
               factors = NULL, grid = 101L, N = 1024L, seed = NULL,
               overwrite = FALSE)
mainEffects(x, name = NULL, strict = FALSE, details = FALSE)
mainEffectNames(x)
plotMainEffects(x, name = NULL, factors = NULL, strict = FALSE, columns = 2L,
                shared_y = TRUE, show_design = TRUE, main = NULL,
                panel_titles = FALSE, las = 0L, ylab = NULL,
                color = "#12355B", lwd = 2.5)
}
\arguments{
  \item{x}{A \code{\linkS4class{BOCampaign}}.}
  \item{name}{Name under which to store the analysis (default
    \code{"main_effects"}). For the accessor and plot, the name to retrieve;
    \code{NULL} (default) gives the most recently stored current analysis.}
  \item{surrogate_name}{Name of the stored surrogate to analyse. \code{NULL}
    (default) uses the most recently stored current surrogate.}
  \item{factors}{Factor labels to compute or plot. \code{NULL} (default) uses
    all factors (all stored factors for the plot).}
  \item{grid}{Number of evenly spaced points along each factor (on the coded
    scale), either one value or one per factor (default 101).}
  \item{N}{Number of points in the background design over which the other
    factors are averaged (default 1024).}
  \item{seed}{\code{NULL} (default) or an integer seed for the background
    design. The global random-number stream is left untouched.}
  \item{overwrite}{Logical; allow replacing a stored analysis with the same
    name. Default \code{FALSE}.}
  \item{strict}{Logical; stop if the requested analysis is historical (see
    \code{\link{summary.BOCampaign}}).}
  \item{details}{Logical; return the complete stored entry instead of the
    table.}
  \item{columns}{Number of panel columns (default 2).}
  \item{shared_y}{Logical; use the same vertical range in all panels (default
    \code{TRUE}), which makes the size of the effects comparable.}
  \item{show_design}{Logical; mark the settings of the completed Sample runs
    along each horizontal axis (default \code{TRUE}).}
  \item{main}{Overall title; \code{NULL} for none.}
  \item{panel_titles}{Logical; give each panel the factor name as title.}
  \item{las}{Orientation of the axis labels, as in
    \code{\link[graphics]{par}}.}
  \item{ylab}{Vertical-axis label; \code{NULL} gives the response name.}
  \item{color, lwd}{Colour and width of the curves.}
}
\details{
The main effect of factor \eqn{i} at value \eqn{x_i}{x_i} is the posterior mean of
the batch-free surface \eqn{\beta_0 + f(x)}{beta0 + f(x)} (see \code{\link{runSurrogate}})
averaged over the other factors, which vary independently and uniformly over
their coded ranges, as in \code{\link{runSobol}}. The average is taken over one
randomized Sobol' set of \code{N} points that is shared by every grid point
and every factor, so the curves are smooth and comparable rather than jittered
by Monte Carlo noise.

The averaging is done on the model scale and the result is transformed back
once. For a log response the curve is therefore \eqn{\exp(E[\log y])}{exp(E[log y])},
a geometric rather than arithmetic mean, consistent with the reported
predictions.

A flat curve means the factor has little effect on average, but it may still
act through interactions; compare with the total-order indices of
\code{\link{runSobol}} or with \code{\link{runResponseSurfaces}}. Factors on a
log scale are drawn on a logarithmic axis.
}
\value{
\code{runMainEffects()} returns the updated campaign. \code{mainEffects()}
returns a data frame in long format with columns \code{factor},
\code{grid_index}, \code{coded_value}, \code{factor_value} (in physical units)
and \code{latent_mean} (the main effect in response units); with
\code{details = TRUE}, the stored entry (\code{table}, \code{surrogate},
\code{settings}, \code{provenance} and \code{call}).
\code{mainEffectNames()} returns the names of all stored main-effect analyses.
\code{plotMainEffects()} draws the plot and returns the plotted rows
invisibly.
}
\seealso{
\code{\link{runSobol}}, \code{\link{runResponseSurfaces}},
\code{\link{predictSurrogate}}.
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

study <- runMainEffects(study, grid = 51, N = 512, seed = 1)
head(mainEffects(study))
plotMainEffects(study)
plotMainEffects(study, factors = c("Capillary_voltage", "Nozzle_voltage"),
                columns = 2, shared_y = FALSE, panel_titles = TRUE)
}
