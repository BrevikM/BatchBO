\name{BatchBO-classes}
\docType{class}
\alias{BatchBO-classes}
\alias{DesignSpace-class}
\alias{BOCampaign-class}
\alias{show,DesignSpace-method}
\alias{show,BOCampaign-method}
\title{The DesignSpace and BOCampaign classes}
\description{
\code{DesignSpace} holds the factors of an experiment and their ranges;
\code{BOCampaign} holds a whole optimization campaign: its specification, the
completed data, and every stored surrogate, proposal and analysis.
}
\section{Objects from the classes}{
Create a \code{DesignSpace} with \code{\link{design_space}} and a
\code{BOCampaign} with \code{\link{campaign}}. Change a campaign only through
the package functions, such as \code{\link{addBatch}} and
\code{\link{runSurrogate}}, which keep its parts consistent; slots should not
be modified directly. A campaign is an ordinary R object, so it is saved and
resumed with \code{\link[base]{saveRDS}} and \code{\link[base]{readRDS}}.
}
\section{Slots of DesignSpace}{
\describe{
  \item{\code{label}}{Factor names.}
  \item{\code{lower}, \code{upper}}{Bounds in physical units.}
  \item{\code{transform}}{Coding of each factor, \code{"linear"} or
    \code{"log10"}.}
}
}
\section{Slots of BOCampaign}{
\describe{
  \item{\code{design}}{The \code{DesignSpace}.}
  \item{\code{response}}{Name of the response column.}
  \item{\code{response_transform}}{\code{"identity"} or \code{"log"}.}
  \item{\code{direction}}{\code{"maximize"} or \code{"minimize"}.}
  \item{\code{controls}}{The reference condition, a one-row data frame in
    physical units, or an empty data frame.}
  \item{\code{control_replicates}}{Number of reference runs per batch.}
  \item{\code{observations}}{The completed runs (see
    \code{\link{observations}}).}
  \item{\code{coded}}{The factor settings of the completed runs on the coded
    scale (see \code{\link{coded}}).}
  \item{\code{standardisation}}{The mean and standard deviation of the
    (transformed) Sample responses and the sign used for minimization (see
    \code{\link{campaign}}). Updated whenever data are added or removed.}
  \item{\code{surrogates}}{Named list of fitted surrogates
    (\code{\link{runSurrogate}}).}
  \item{\code{analyses}}{Named list of recommendations and analyses
    (\code{\link{recommend}}, \code{\link{runSobol}},
    \code{\link{runMainEffects}}, \code{\link{runResponseSurfaces}},
    \code{\link{runCumulativeBest}}).}
  \item{\code{proposals}}{Named list of proposed run sheets
    (\code{\link{proposeSpaceFilling}}, \code{\link{proposeBatch}}).}
}
Every entry of the last three slots carries a provenance record, used to
mark it as current or historical (see \code{\link{summary.BOCampaign}}).
}
\section{Methods}{
\describe{
  \item{show}{\code{signature(object = "DesignSpace")}: print the factors and
    their ranges.}
  \item{show}{\code{signature(object = "BOCampaign")}: print a short overview
    of the campaign and the names of its stored results.}
  \item{as.data.frame}{\code{signature(x = "DesignSpace")}: see
    \code{\link{design_space}}.}
  \item{summary}{\code{\link{summary.BOCampaign}}.}
}
}
\seealso{
\code{\link{design_space}}, \code{\link{campaign}},
\code{\link{campaign-accessors}}.
}
\examples{
space <- design_space(
  def_factor("Temperature", lower = 20, upper = 80),
  def_factor("Concentration", lower = 0.1, upper = 100, transform = "log10")
)
space
is(space, "DesignSpace")

study <- campaign(space, response = "Yield")
study
slotNames(study)
}
\keyword{classes}
