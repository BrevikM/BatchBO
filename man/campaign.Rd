\name{campaign}
\alias{campaign}
\title{Create an optimization campaign}
\description{
Create the empty, resumable object that carries an optimization campaign: the
design space, the measured response and whether to maximize or minimize it, and
the fixed reference condition repeated in every batch. Completed batches are
added with \code{\link{addBatch}}; fitted surrogates, proposals and analyses
are stored in the same object as the campaign progresses.
}
\usage{
campaign(design, response, transform = c("identity", "log"),
         direction = c("maximize", "minimize"), controls = NULL,
         replicates = 2)
}
\arguments{
  \item{design}{A \code{\linkS4class{DesignSpace}} from \code{\link{design_space}}.}
  \item{response}{Name of the measured response column, for example
    \code{"Peptides"}. It must differ from the factor names and from
    \code{Batch}, \code{Type} and \code{Run}.}
  \item{transform}{\code{"identity"} (default) to model the response as
    measured, or \code{"log"} to model its natural logarithm. The response must
    then be strictly positive.}
  \item{direction}{\code{"maximize"} (default) or \code{"minimize"}.}
  \item{controls}{\code{NULL} for no reference condition, or a one-row data
    frame giving the fixed reference condition in physical units, with one
    column per factor.}
  \item{replicates}{Number of reference runs added to every proposed batch
    (default 2). Ignored when \code{controls} is \code{NULL}.}
}
\details{
\strong{Reference condition.} A reference condition is strongly recommended
when batch effects are expected. Runs of the same condition within one batch
measure residual noise directly, and the same condition across batches anchors
the batches to a common point; together they let the model separate batch
offsets from residual noise and from the response surface. Every proposed run
sheet contains \code{replicates} reference runs (two or more are needed to
measure residual noise within a batch), and \code{\link{addBatch}} checks that
reference runs use exactly this condition.
The reference may lie outside the declared bounds (for example a vendor default
setting); it then informs the model but cannot itself be proposed or
recommended.

\strong{Response scale.} Before fitting, the (optionally log-transformed)
response is standardized to mean 0 and standard deviation 1 using the Sample
runs only, and its sign is flipped when minimizing, so that larger values are
always better internally. Reference runs enter the model but not these
constants. All reported predictions, recommendations and analyses are converted
back to the units of the response. For a log response, reported central values
are back-transformed posterior means, i.e. posterior medians on the response
scale; \code{\link{predictSurrogate}} also reports the mean.

The design, response, transform, direction and reference condition are fixed
once the campaign is created.
}
\value{
An empty \code{\linkS4class{BOCampaign}}. Printing it gives a short overview.
}
\seealso{
\code{\link{design_space}}, \code{\link{proposeSpaceFilling}} for the initial
batches, \code{\link{addBatch}}, \code{\link{campaign-accessors}},
\code{\link{summary.BOCampaign}}.
}
\examples{
space <- design_space(
  def_factor("Temperature", lower = 20, upper = 80),
  def_factor("Concentration", lower = 0.1, upper = 100, transform = "log10")
)

study <- campaign(
  design     = space,
  response   = "Yield",
  direction  = "maximize",
  controls   = data.frame(Temperature = 50, Concentration = 10),
  replicates = 2
)
study
}
