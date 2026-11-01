\name{runResponseSurfaces}
\alias{runResponseSurfaces}
\alias{responseSurfaces}
\alias{responseSurfacesNames}
\alias{plotResponseSurfaces}
\title{Two-factor views of the fitted surface}
\description{
Show the fitted response over pairs of factors. \code{runResponseSurfaces()}
evaluates the surface of a stored surrogate on a grid for every pair of the
chosen factors and stores the result in the campaign;
\code{responseSurfaces()} retrieves it and \code{plotResponseSurfaces()} draws
filled contour plots.
}
\usage{
runResponseSurfaces(x, name = "response_surfaces", surrogate_name = NULL,
                    factors = NULL, type = c("marginal", "slice"), grid = 41L,
                    N = 512L, at = "centre", recommendation_name = NULL,
                    seed = NULL, overwrite = FALSE)
responseSurfaces(x, name = NULL, strict = FALSE, details = FALSE)
responseSurfacesNames(x)
plotResponseSurfaces(x, name = NULL, strict = FALSE, columns = 4L,
                     shared_scale = TRUE, panel_titles = FALSE,
                     show_observations = FALSE, recommendation_name = NULL,
                     main = NULL,
                     colors = grDevices::hcl.colors(64L, "Blues 3"),
                     contours = 10L, contour_color = "#303030",
                     color_key = TRUE, point_legend = TRUE, zlim = NULL)
}
\arguments{
  \item{x}{A \code{\linkS4class{BOCampaign}}.}
  \item{name}{Name under which to store the analysis (default
    \code{"response_surfaces"}). For the accessor and plot, the name to
    retrieve; \code{NULL} (default) gives the most recently stored current
    analysis.}
  \item{surrogate_name}{Name of the stored surrogate to analyse. \code{NULL}
    (default) uses the most recently stored current surrogate.}
  \item{factors}{Two or more factor labels; a surface is computed for every
    pair. \code{NULL} (default) uses all factors.}
  \item{type}{\code{"marginal"} (default) to average over the factors not
    shown, or \code{"slice"} to hold them fixed at \code{at}.}
  \item{grid}{Number of grid points along each displayed factor, one value or
    two (default 41, giving \eqn{41 \times 41}{41 x 41} points).}
  \item{N}{For \code{type = "marginal"}, number of points in the background
    design over which the other factors are averaged (default 512).}
  \item{at}{For \code{type = "slice"}, where to hold the other factors:
    \code{"centre"} (default, the middle of each coded range),
    \code{"recommendation"} (a stored current recommendation), or a one-row
    data frame giving a condition in physical units inside the bounds.}
  \item{recommendation_name}{Name of the stored recommendation to use with
    \code{at = "recommendation"}, or to mark on the plot; \code{NULL} gives the
    most recent current one in \code{runResponseSurfaces()} and no mark in
    the plot.}
  \item{seed}{\code{NULL} (default) or an integer seed for the background
    design. The global random-number stream is left untouched.}
  \item{overwrite}{Logical; allow replacing a stored analysis with the same
    name. Default \code{FALSE}.}
  \item{strict}{Logical; stop if the requested analysis is historical (see
    \code{\link{summary.BOCampaign}}).}
  \item{details}{Logical; return the complete stored entry instead of the
    surfaces.}
  \item{columns}{Number of panel columns when several surfaces are drawn
    (default 4).}
  \item{shared_scale}{Logical; use one colour range for all panels (default
    \code{TRUE}), so that colours mean the same response everywhere.}
  \item{panel_titles}{Logical; give each panel the names of its two factors
    as title.}
  \item{show_observations}{Logical; mark the completed Sample runs (projected
    onto each pair of factors). Default \code{FALSE}.}
  \item{main}{Title; \code{NULL} for none.}
  \item{colors}{Colour ramp of the surface, from low to high response.}
  \item{contours}{Approximate number of contour lines.}
  \item{contour_color}{Colour of the contour lines.}
  \item{color_key}{Logical; draw a colour key. Only for a single surface.}
  \item{point_legend}{Logical; explain the plotted points. Only for a single
    surface.}
  \item{zlim}{Colour range as two values in response units; \code{NULL}
    (default) uses the range of the surfaces.}
}
\details{
Values are the posterior mean of the batch-free surface
\eqn{\beta_0 + f(x)}{beta0 + f(x)} (see \code{\link{runSurrogate}}), in response units.
\describe{
  \item{Marginal}{The average over the factors not shown, which vary
    independently and uniformly over their coded ranges, using one
    randomized Sobol' set of \code{N} points shared by all grid points and
    pairs. This is the two-factor analogue of \code{\link{runMainEffects}}
    and shows the average behaviour of the pair, including their
    interaction. As there, averaging is done on the model scale.}
  \item{Slice}{The surface with the other factors held at one condition, for
    example at the recommendation, which shows what changing two factors
    around that condition would do. A slice can look very different from the
    marginal surface when the factors interact.}
}
Grid points are evenly spaced on the coded scale, and factors on a log scale
are drawn on logarithmic axes. With \eqn{k} factors there are
\eqn{k(k-1)/2} surfaces, each needing \eqn{N} evaluations per grid point for a
marginal surface, so for many factors a smaller \code{grid} or \code{N}, or a
subset of \code{factors}, keeps the computation short.

A single surface is drawn with \code{\link[graphics]{filled.contour}}, which
sets up its own layout. Several surfaces are drawn as panels with
\code{\link[graphics]{image}}, without colour key.
}
\value{
\code{runResponseSurfaces()} returns the updated campaign.
\code{responseSurfaces()} returns a named list with one data frame per pair
(named \code{"<factor1>__<factor2>"}), with the two factor columns in physical
units and \code{latent_mean}, the surface in response units; with
\code{details = TRUE}, the stored entry (\code{surfaces}, \code{surrogate},
\code{settings}, \code{slice_condition}, \code{provenance} and \code{call}).
\code{responseSurfacesNames()} returns the names of all stored response-surface
analyses. \code{plotResponseSurfaces()} draws the plot and returns the
surfaces invisibly.
}
\seealso{
\code{\link{runMainEffects}}, \code{\link{runSobol}} with
\code{second_order = TRUE} to find which pairs interact,
\code{\link{recommend}}.
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
study <- recommend(study, name = "best", seed = 1)

## One marginal surface, with the data and the recommendation
pair <- c("Capillary_voltage", "Nozzle_voltage")
study <- runResponseSurfaces(study, name = "pair", factors = pair,
                             grid = 31, N = 256, seed = 1)
plotResponseSurfaces(study, "pair", show_observations = TRUE,
                     recommendation_name = "best")

\donttest{
## All pairs of three factors, sliced through the recommendation
study <- runResponseSurfaces(study, name = "slices",
                             factors = c(pair, "Drying_gas_temperature"),
                             type = "slice", at = "recommendation",
                             recommendation_name = "best", grid = 31)
plotResponseSurfaces(study, "slices", columns = 3, panel_titles = TRUE)
names(responseSurfaces(study, "slices"))
}
}
