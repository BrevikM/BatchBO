\name{campaign-accessors}
\alias{campaign-accessors}
\alias{design}
\alias{observations}
\alias{coded}
\alias{nruns}
\alias{nbatches}
\alias{surrogateNames}
\alias{bestObserved}
\title{Inspect the contents of a campaign}
\description{
Simple accessors for the design space, the completed data, their coded form,
the numbers of runs and batches, the names of stored surrogates, and the best
measured Sample run.
}
\usage{
design(x)
observations(x)
coded(x)
nruns(x)
nbatches(x)
surrogateNames(x)
bestObserved(x)
}
\arguments{
  \item{x}{A \code{\linkS4class{BOCampaign}}.}
}
\details{
\code{observations()} returns the runs exactly as entered with
\code{\link{addBatch}} (in physical units, with \code{Type} normalized to
\code{"Sample"} or \code{"Ref"}), in the order the batches were added.
\code{coded()} returns the same factor settings coded to the unit cube, which
is the scale the Gaussian process works on (see \code{\link{design_space}}).

\code{bestObserved()} considers Sample runs only and returns the largest
measured response, or the smallest when the campaign minimizes. It is the raw
measured best and so includes that run's batch offset and measurement noise;
\code{\link{runCumulativeBest}} gives a batch-adjusted version, and
\code{\link{recommend}} the condition the model expects to be best.
}
\value{
\code{design()} returns the \code{\linkS4class{DesignSpace}};
\code{observations()} a data frame with one row per run; \code{coded()} a
numeric matrix with one row per run and one column per factor, in
\eqn{[0, 1]} for conditions inside the bounds; \code{nruns()} and
\code{nbatches()} integers; \code{surrogateNames()} the names of all stored
surrogates (current and historical); \code{bestObserved()} a list with the
best \code{value}, its \code{row} in \code{observations()} and its
\code{batch}, or \code{NULL} when no Sample run has been entered.
}
\seealso{
\code{\link{summary.BOCampaign}} for an overview of stored results,
\code{\link{proposalNames}}, \code{\link{sobolNames}} and the other
\code{*Names()} functions for stored analyses.
}
\examples{
surface <- demoSurface()
study <- campaign(surface$design, response = "Ion_response")
study <- proposeSpaceFilling(study, batch = "B1", q = 6, seed = 1)
sheet <- proposal(study, "B1")
sheet$Ion_response <- simulateBatch(surface, sheet, seed = 1)
study <- addBatch(study, sheet)

nruns(study)
nbatches(study)
head(observations(study))
head(coded(study))
bestObserved(study)
}
