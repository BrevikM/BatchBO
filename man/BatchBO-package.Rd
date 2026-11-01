\name{BatchBO-package}
\alias{BatchBO-package}
\alias{BatchBO}
\docType{package}
\title{Batch-aware Bayesian optimization for laboratory experiments}
\description{
\pkg{BatchBO} optimizes costly laboratory experiments that are acquired in
batches. Measurements from the same batch share an unknown offset, for example
from a fresh reagent lot, a recalibrated instrument or a new chromatographic
column. The package models these offsets explicitly with a batch-aware Gaussian
process, so that new conditions are chosen, and the final condition is
recommended, on the reproducible batch-free response surface rather than on
whichever batch happened to run high.
}
\details{
A campaign normally proceeds as follows.
\enumerate{
  \item Declare the factors with \code{\link{def_factor}} and
    \code{\link{design_space}}, and the response, optimization direction and
    a fixed reference condition with \code{\link{campaign}}.
  \item Generate an initial space-filling design with
    \code{\link{proposeSpaceFilling}} and retrieve the run sheets with
    \code{\link{proposal}}.
  \item Run the batches, add the measured responses to the sheets and enter
    them with \code{\link{addBatch}}.
  \item Fit the surrogate with \code{\link{runSurrogate}} and propose the next
    batch with \code{\link{proposeBatch}}. Repeat steps 3 and 4.
  \item Recommend a condition with \code{\link{recommend}}, and interpret the
    fitted surface with \code{\link{runCumulativeBest}},
    \code{\link{runSobol}}, \code{\link{runMainEffects}} and
    \code{\link{runResponseSurfaces}}.
}

\strong{The model.} With factor settings coded to the unit cube, the
measurement of run \eqn{r} in batch \eqn{b(r)} is modelled as
\deqn{y_r = \beta_0 + f(x_r) + u_{b(r)} + \varepsilon_r,}{y_r = beta_0 + f(x_r) + u_b(r) + e_r,}
where \eqn{f} is a zero-mean Gaussian process with an ARD Matern or squared
exponential kernel and variance \eqn{\sigma_f^2}{sigma_f^2}, the batch offsets \eqn{u_b}
are independent \eqn{N(0, \sigma_b^2)}{N(0, sigma_b^2)} and the residuals \eqn{\varepsilon_r}{e_r}
are independent \eqn{N(0, \sigma^2)}{N(0, sigma^2)}. Hyperparameters are estimated by
maximizing the marginal likelihood (ML-II). Acquisition, recommendation and all
analyses use the posterior of the batch-free surface \eqn{\beta_0 + f}{beta0 + f}.

\strong{Reference runs.} Repeating one fixed reference condition in every batch
gives a direct estimate of residual noise within a batch and anchors the batches
to a common point, which is what lets the model separate batch offsets from
residual noise and from the surface.

\strong{Stored results.} A campaign object carries everything needed to resume
it: the observations, every fitted surrogate, every proposed batch and every
analysis, each stored under a name together with a snapshot of the data it was
computed from. \code{\link{summary.BOCampaign}} lists them and marks each as
\emph{current} (computed from the present data) or \emph{historical}. Nothing
is deleted when data change; historical results stay available by name.

\strong{Reproducibility.} Every function that draws random numbers accepts a
\code{seed}. A supplied seed is used locally: the global random-number state is
left unchanged.

The vignette \code{vignette("hartmann_ESI", package =
"BatchBO")} walks through a complete simulated campaign.
}
\references{
Rasmussen, C. E. and Williams, C. K. I. (2006). \emph{Gaussian Processes for
Machine Learning}. MIT Press.
}
\seealso{
\code{\link{campaign}}, \code{\link{BatchBO-classes}},
\code{\link{demoSurface}} for a synthetic laboratory to try the package on.
}
\keyword{package}
