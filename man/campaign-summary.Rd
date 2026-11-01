\name{summary.BOCampaign}
\alias{summary.BOCampaign}
\alias{print.summary.BOCampaign}
\title{Summarize a campaign and its stored results}
\description{
List the campaign's data and every stored surrogate, analysis and proposal,
marking each as \emph{current} or \emph{historical}.
}
\usage{
\method{summary}{BOCampaign}(object, ...)
\method{print}{summary.BOCampaign}(x, ...)
}
\arguments{
  \item{object}{A \code{\linkS4class{BOCampaign}}.}
  \item{x}{A \code{summary.BOCampaign} object.}
  \item{...}{Ignored.}
}
\details{
Every stored result keeps a snapshot of what it was computed from: the design
space, the response specification and the complete coded data (factor
settings, responses, batch labels and run types). A result is \strong{current}
when its snapshot still matches the campaign, and \strong{historical}
otherwise. In addition,
\itemize{
  \item proposals become historical when the declared reference condition
    changes, since it is part of their run sheets;
  \item proposals and analyses become historical when the surrogate they were
    computed from is no longer stored, or has been refitted under the same
    name;
  \item an analysis (recommendation, Sobol indices, main effects, response
    surfaces, cumulative best) inherits the snapshot of its surrogate, so it is
    current exactly while that surrogate is.
}
Adding a batch therefore makes all existing results historical; this is
expected, and they remain available by name. Accessors called without a name,
such as \code{surrogate(study)} or \code{sobolIndices(study)}, return the most
recent \emph{current} result of that kind. With \code{strict = TRUE}, any
accessor stops with the reason when the requested result is historical.
}
\value{
\code{summary()} returns an object of class \code{summary.BOCampaign}: a list
with elements \code{campaign} (factor names, response, transform, direction,
numbers of runs and batches) and three data frames with one row per stored
result. \code{surrogates} has columns \code{name}, \code{model}
(\code{"batch-aware"} or \code{"naive"}), \code{kernel}; \code{analyses} has
\code{name}, \code{type}, \code{surrogate}; \code{proposals} has
\code{name} (the batch label), \code{method}, \code{surrogate} (empty for
space-filling proposals). All three also have \code{batches} and \code{runs}
(the data the result was computed from), \code{stored_at} and \code{status}
(\code{"current"} or \code{"historical"}). \code{print()} returns
its argument invisibly.
}
\seealso{
\code{\link{campaign-accessors}}, \code{\link{surrogate}}.
}
\examples{
surface <- demoSurface()
study <- campaign(surface$design, response = "Ion_response")
study <- proposeSpaceFilling(study, batch = c("B1", "B2"), q = 5, seed = 1)
summary(study)
}
