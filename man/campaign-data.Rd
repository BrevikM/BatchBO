\name{campaign-data}
\alias{campaign-data}
\alias{addBatch}
\alias{removeBatch}
\title{Add or remove completed batches}
\description{
\code{addBatch()} enters measured runs into a campaign, one batch or several at
once. \code{removeBatch()} removes batches by their labels.
}
\usage{
addBatch(x, data, overwrite = FALSE)
removeBatch(x, batches)
}
\arguments{
  \item{x}{A \code{\linkS4class{BOCampaign}}.}
  \item{data}{The completed runs: a data frame, or a list of data frames (for
    example one per batch) with the same columns. Usually the run sheets from
    \code{\link{proposal}} with the measured response filled in.}
  \item{overwrite}{Replace a batch that is already in the campaign with
    different data. Defaults to \code{FALSE}.}
  \item{batches}{Character vector of batch labels to remove.}
}
\details{
\strong{Required columns.} \code{data} must contain exactly these columns:
\describe{
  \item{\code{Batch}}{The batch label. Labels are free text; a proposed batch
    may be renamed here, for example to record the date it was run.}
  \item{\code{Type}}{\code{"Sample"} for a proposed condition, or \code{"Ref"}
    (also accepted: \code{"Reference"}, any capitalization) for a reference
    run.}
  \item{one column per factor}{The settings actually run, in physical units.}
  \item{the response}{The measured values, named as in \code{\link{campaign}}.
    All values must be finite, and positive for a log response.}
  \item{\code{Run}}{Optional run order within the batch: positive whole
    numbers, unique within each batch.}
}
Any other column is refused, so the data cannot silently carry unintended
information into the model. Reference runs must use the campaign's reference
condition (to a relative tolerance of \eqn{10^{-6}}{10^-6}).

\strong{Batches already in the campaign.} When a batch label is already
present, the incoming runs are compared with the stored ones (the same Type,
factor settings and responses, in any row order, to a small numerical
tolerance). An identical batch is skipped with a message, so re-entering the
complete data set, or a file written to CSV and read back, adds only the new
batches. A batch with different data is refused unless
\code{overwrite = TRUE}, in which case it replaces the stored batch in its
original position.

\strong{Checks after adding.} A warning is given for any new batch without
reference runs when the campaign declares a reference condition, since such a
batch cannot be anchored to the others. A message is given when the reference
condition lies outside the declared bounds.

Adding or removing data does not delete stored surrogates, proposals or
analyses; they become \emph{historical} (see \code{\link{summary.BOCampaign}}).
Refit the surrogate with \code{\link{runSurrogate}} before proposing again.
}
\value{
The updated \code{\linkS4class{BOCampaign}}.
}
\seealso{
\code{\link{proposal}} for the run sheets, \code{\link{observations}} for the
stored data, \code{\link{runSurrogate}}.
}
\examples{
surface <- demoSurface(batch_sd = 80, residual_sd = 40)
reference <- data.frame(Capillary_voltage = 3500, Nebulizer_pressure = 35,
                        Nozzle_voltage = 1000, Drying_gas_flow = 10,
                        Drying_gas_temperature = 275, Sheath_gas_flow = 10)
study <- campaign(surface$design, response = "Ion_response",
                  controls = reference, replicates = 2)
study <- proposeSpaceFilling(study, batch = c("B1", "B2"), q = 4, seed = 1)

sheet <- proposal(study, c("B1", "B2"))
sheets <- split(sheet, sheet$Batch)
for (i in seq_along(sheets)) {
  sheets[[i]]$Ion_response <- simulateBatch(surface, sheets[[i]], seed = i)
}

study <- addBatch(study, sheets)   # a list of batches
study <- addBatch(study, sheets)   # identical batches are skipped
nbatches(study)

study <- removeBatch(study, "B2")
nbatches(study)
}
