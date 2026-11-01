\name{demonstration}
\encoding{UTF-8}
\alias{demonstration}
\alias{demoSurface}
\alias{print.DemoSurface}
\alias{simulateBatch}
\alias{demoResponse}
\alias{demoTruth}
\alias{demoCompare}
\title{A simulated laboratory for trying out the package}
\description{
A known six-factor response surface with simulated batch offsets and
measurement noise, for examples, teaching and testing optimization settings.
\code{demoSurface()} creates the simulated laboratory, \code{simulateBatch()}
"measures" a run sheet, and \code{demoResponse()} and \code{demoTruth()} give
the true surface and its optimum, which a real campaign never has.
\code{demoCompare()} compares a condition, such as a recommendation, with that
optimum. Declare
the campaign itself with \code{\link{campaign}}, using the design space of
the surface; see Examples.
}
\usage{
demoSurface(batch_sd = 80, residual_sd = 40,
            interpretation = c("mass_spectrometry", "unit_cube"))
\method{print}{DemoSurface}(x, ...)
simulateBatch(x, run_sheet, seed = NULL, response = NULL)
demoResponse(x, newdata)
demoTruth(x, response = NULL)
demoCompare(x, condition, response = NULL)
}
\arguments{
  \item{batch_sd}{Standard deviation of the simulated batch offset (default
    80).}
  \item{residual_sd}{Standard deviation of the simulated run-to-run noise
    (default 40).}
  \item{interpretation}{\code{"mass_spectrometry"} (default) to present the
    factors as the settings of an electrospray ion source, or
    \code{"unit_cube"} for factors \code{x1}, \ldots, \code{x6} on
    \eqn{[0, 1]}. Only names and ranges differ; the surface is the same.}
  \item{x}{A \code{DemoSurface} from \code{demoSurface()}.}
  \item{run_sheet}{A run sheet, as returned by \code{\link{proposal}}, for one
    batch. Its response column, if present, must be empty.}
  \item{seed}{\code{NULL} (default) or an integer seed for the simulated
    offset and noise. The global random-number stream is left untouched.}
  \item{response}{Name of the response, if it differs from the one declared
    by \code{x} (\code{"Ion_response"}, or \code{"Response"} for the unit
    cube).}
  \item{newdata}{A data frame of conditions in physical units.}
  \item{condition}{A one-row data frame with a condition in physical units,
    for example the result of \code{\link{recommendation}}. Other columns are
    ignored.}
  \item{...}{Ignored.}
}
\details{
\strong{Surface.} The true response is \eqn{1000 + 500\,h(z)}{1000 + 500 h(z)},
where \eqn{h} is the six-dimensional Hartmann function (Dixon and Szeg\enc{ö}{o},
1978) and \eqn{z} the settings rescaled to \eqn{[0, 1]^6}{[0, 1]^6}. It has several local
maxima and a global maximum of about 2661, at the condition returned by
\code{demoTruth()}. The mass-spectrometry interpretation uses the following
factors and ranges (all on a linear scale):
\tabular{lll}{
  \code{Capillary_voltage} \tab 2500--4500 \tab V\cr
  \code{Nebulizer_pressure} \tab 20--50 \tab psi\cr
  \code{Nozzle_voltage} \tab 0--2000 \tab V\cr
  \code{Drying_gas_flow} \tab 8--12 \tab L/min\cr
  \code{Drying_gas_temperature} \tab 200--350 \tab \eqn{^\circ}{deg}C\cr
  \code{Sheath_gas_flow} \tab 8--12 \tab L/min
}
These names and ranges are illustrative; the surface does not model a real
instrument.

\strong{Simulated measurements.} Each call of \code{simulateBatch()} draws
one batch offset from a normal distribution with mean 0 and standard deviation
\code{batch_sd}, shared by all runs of the sheet, and independent noise with
standard deviation \code{residual_sd} for every run, and adds them to the true
response. Call it once per batch.
}
\value{
\code{demoSurface()} returns an object of class \code{DemoSurface}: a list
with the \code{design} (a \code{\linkS4class{DesignSpace}}), the default
\code{response} name, \code{batch_sd}, \code{residual_sd},
\code{interpretation} and the \code{units} of the factors; printing it lists
these. \code{simulateBatch()} returns the simulated
responses, in the row order of \code{run_sheet}. \code{demoResponse()} returns
the true response at each row of \code{newdata}. \code{demoTruth()} returns a
one-row data frame with the optimal condition and its true response.
\code{demoCompare()} returns a data frame with columns \code{condition} and
\code{optimum}, and one row per factor followed by the true (batch-free)
response and the \code{regret}, the shortfall of the condition's true response
from the optimum (for the optimum, 0).
}
\references{
Dixon, L. C. W. and Szeg\enc{ö}{o}, G. P. (1978). The global optimization
problem: an introduction. In Dixon, L. C. W. and Szeg\enc{ö}{o}, G. P. (eds.),
\emph{Towards Global Optimisation 2}, pp. 1--15. North-Holland.
}
\seealso{
\code{\link{campaign}}, \code{\link{proposeSpaceFilling}},
\code{\link{addBatch}}; the package vignette for a complete simulated campaign.
}
\examples{
surface <- demoSurface(batch_sd = 80, residual_sd = 40)
surface
demoTruth(surface)

## Measure one batch of a campaign over the same design space
study <- campaign(surface$design, response = "Ion_response")
study <- proposeSpaceFilling(study, batch = "B1", q = 6, seed = 1)
sheet <- proposal(study, "B1")
sheet$Ion_response <- simulateBatch(surface, sheet, seed = 1)
study <- addBatch(study, sheet)

## Measured versus true: the difference is the batch offset plus noise
cbind(measured = sheet$Ion_response,
      true = demoResponse(surface, sheet))

## How close is the best measured condition to the optimum?
best <- observations(study)[bestObserved(study)$row, ]
demoCompare(surface, best)
}
