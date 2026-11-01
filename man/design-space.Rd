\name{design-space}
\alias{design-space}
\alias{def_factor}
\alias{design_space}
\alias{factorNames}
\alias{as.data.frame,DesignSpace-method}
\title{Declare the experimental factors and their ranges}
\description{
\code{def_factor()} declares one continuous experimental factor by its name and
its lower and upper bound in physical units. \code{design_space()} combines one
or more factors into the design space of a campaign. \code{factorNames()}
returns the factor names, and \code{as.data.frame()} lists the factors with
their bounds and coding.
}
\usage{
def_factor(label, lower, upper, transform = c("linear", "log10"))
design_space(...)
factorNames(x)
\S4method{as.data.frame}{DesignSpace}(x, row.names = NULL, optional = FALSE, ...)
}
\arguments{
  \item{label}{Name of the factor. It becomes a column name in run sheets and
    in the completed data, so it must be unique (ignoring capitalization) and
    must not be one of the reserved column names listed in Details.}
  \item{lower, upper}{Lower and upper bound of the factor in physical units,
    with \code{lower < upper}. Both must be positive for
    \code{transform = "log10"}.}
  \item{transform}{How the factor is coded internally: \code{"linear"}
    (default) or \code{"log10"}.}
  \item{...}{For \code{design_space()}, one or more factors created by
    \code{def_factor()}, in the order they should appear. For
    \code{as.data.frame()}, ignored.}
  \item{x}{A \code{DesignSpace} or, for \code{factorNames()}, also a
    \code{BOCampaign}.}
  \item{row.names, optional}{Standard \code{as.data.frame()} arguments.
    \code{row.names} is passed to the result; \code{optional} is ignored.}
}
\details{
Internally every factor is coded to \eqn{[0, 1]}: the lower bound maps to 0
and the upper bound to 1. With \code{transform = "linear"} the coding is linear
in the physical value. With \code{transform = "log10"} it is linear in
\eqn{\log_{10}}{log10} of the value, which suits factors whose effect is naturally
multiplicative, such as concentrations or ion-accumulation targets spanning
orders of magnitude. The transform only affects the internal coding: bounds,
run sheets, completed data and every reported condition stay in physical
units, and plots of a log10-coded factor use a logarithmic axis.

The Gaussian process uses one length scale per factor on the coded scale, so
the coding determines what "a short distance" means for each factor.

The names \code{Batch}, \code{Type} and \code{Run} (run-sheet columns), and
\code{prediction}, \code{response_mean}, \code{prediction_lower},
\code{prediction_upper}, \code{model_mean}, \code{model_sd} and
\code{latent_mean} (output columns) are reserved and cannot be used as factor
labels.
}
\value{
\code{def_factor()} returns a factor declaration for \code{design_space()}.
\code{design_space()} returns a \code{\linkS4class{DesignSpace}}.
\code{factorNames()} returns a character vector of factor names.
\code{as.data.frame()} returns a data frame with columns \code{factor},
\code{lower}, \code{upper} and \code{transform}, one row per factor.
}
\seealso{
\code{\link{campaign}}, which combines a design space with the response and
reference condition.
}
\examples{
space <- design_space(
  def_factor("Temperature", lower = 20, upper = 80),
  def_factor("Concentration", lower = 0.1, upper = 100, transform = "log10"),
  def_factor("pH", lower = 4, upper = 9)
)
space
factorNames(space)
as.data.frame(space)
}
