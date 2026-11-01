\name{runSurrogate}
\alias{runSurrogate}
\alias{surrogate}
\title{Fit and retrieve a batch-aware Gaussian-process surrogate}
\description{
\code{runSurrogate()} fits a Gaussian process to the completed runs of a
campaign, with a random batch offset in addition to residual noise, and stores
it in the campaign under a name. \code{surrogate()} retrieves a stored fit.
}
\usage{
runSurrogate(x, name, kernel = c("matern52", "matern32", "rbf"),
             naive = FALSE, n_starts = 8, maxit = 200L, jitter = 1e-08,
             seed = NULL, overwrite = FALSE, ...)
surrogate(x, name = NULL, strict = FALSE)
}
\arguments{
  \item{x}{A \code{\linkS4class{BOCampaign}}.}
  \item{name}{Name under which to store the surrogate (required), for example
    \code{"B3_fit"}. For \code{surrogate()}, the name to retrieve; \code{NULL}
    (default) gives the most recently stored current surrogate.}
  \item{kernel}{Correlation function: \code{"matern52"} (default),
    \code{"matern32"} or \code{"rbf"} (squared exponential).}
  \item{naive}{Logical; fit a batch-naive model without the batch offset, for
    comparison. Default \code{FALSE}.}
  \item{n_starts}{Number of random starting points for the
    hyperparameter optimization (default 8).}
  \item{maxit}{Maximum number of L-BFGS-B iterations per start (default 200).}
  \item{jitter}{Small constant added to the diagonal of the covariance matrix
    for numerical stability (default \code{1e-8}).}
  \item{seed}{\code{NULL} (default) or an integer seed for the starting points.
    The global random-number stream is left untouched.}
  \item{overwrite}{Logical; allow replacing a stored surrogate with the same
    name. Default \code{FALSE}.}
  \item{...}{Further optimization settings; see Details.}
  \item{strict}{Logical; stop if the requested surrogate is historical (see
    \code{\link{summary.BOCampaign}}).}
}
\details{
\strong{Model.} With factor settings \eqn{x} coded to the unit cube (see
\code{\link{design_space}}), the standardized response of run \eqn{i} in batch
\eqn{b} is modelled as
\deqn{y_i = \beta_0 + f(x_i) + u_b + \varepsilon_i,}{y_i = beta0 + f(x_i) + u_b + e_i,}
where \eqn{f} is a zero-mean Gaussian process with variance
\eqn{\sigma_f^2}{sigma_f^2} and correlation function \eqn{R} with one length scale per
factor (automatic relevance determination), \eqn{u_b \sim N(0, \sigma_b^2)}{u_b ~ N(0, sigma_b^2)}
is a batch offset shared by all runs of batch \eqn{b}, and
\eqn{\varepsilon_i \sim N(0, \sigma^2)}{e_i ~ N(0, sigma^2)} is residual noise.
The covariance of the responses is therefore
\deqn{V = \sigma_f^2 R + \sigma_b^2 S + \sigma^2 I,}{V = sigma_f^2 R + sigma_b^2 S + sigma^2 I,}
where \eqn{S_{ij} = 1}{S_ij = 1} when runs \eqn{i} and \eqn{j} share a batch. The
constant \eqn{\beta_0}{beta0} is estimated by generalized least squares. Predictions,
proposals and analyses all concern \eqn{\beta_0 + f(x)}{beta0 + f(x)}, the batch-free
response surface. With \code{naive = TRUE} the batch term is omitted, so batch
offsets are absorbed by \eqn{f} and the residual noise.

Before fitting, the response is log-transformed if requested, standardized
using the Sample runs and sign-flipped when minimizing (see
\code{\link{campaign}}). Reference runs are included in the fit.

\strong{Estimation.} The length scales and the three variances are estimated
by maximizing the marginal likelihood (type-II maximum likelihood) with
L-BFGS-B and analytic gradients, on the log scale and within bounds. The
optimization is started from \code{n_starts} points drawn uniformly (on the
log scale) within the bounds, and the best solution is kept. A warning is
given when the selected solution did not converge or lies on a bound; a length
scale at its upper bound usually means the factor has little detectable
effect, and a variance at its lower bound that the component is not
detectable.

Settings accepted through \code{...}:
\describe{
  \item{\code{lower}, \code{upper}}{Bounds as four values, for the length
    scales (on the coded scale), signal variance, batch variance and residual
    variance (both on the standardized scale). Defaults
    \code{c(0.05, 1e-4, 1e-8, 1e-8)} and \code{c(50, 5, 5, 5)}.}
  \item{\code{hp_start}}{A list with elements \code{ell} (one value per
    factor), \code{sigma_f2}, \code{sigma_b2} and \code{sigma2}, used as the
    first starting point, for example the hyperparameters of a previous fit,
    \code{surrogate(study, "B2_fit")$hp}.}
}

At least three completed runs, and two Sample runs with different responses,
are needed. A batch-aware fit to a single batch is possible but warns, since
the batch and residual variances cannot then be told apart.

The fit is stored with a snapshot of the data it used (see
\code{\link{summary.BOCampaign}}); after further batches are added it becomes
historical, and a new surrogate should be fitted under a new name.
}
\value{
\code{runSurrogate()} returns the updated campaign. \code{surrogate()} returns
the fitted model, an object of class \code{"GP"} (see
\code{\link{summary.GP}}), a list containing the coded factor settings
\code{X}, standardized responses \code{y}, \code{batch} labels, the
hyperparameters \code{hp} (\code{ell}, \code{sigma_f2}, \code{sigma_b2},
\code{sigma2}, \code{kernel} and optimization diagnostics), \code{beta0}, the
Cholesky factor \code{U} of \eqn{V}, the negative log marginal likelihood
\code{nll} and the standardization constants.
}
\references{
Rasmussen, C. E. and Williams, C. K. I. (2006). \emph{Gaussian Processes for
Machine Learning}. MIT Press.
}
\seealso{
\code{\link{summary.GP}} for variance components and length scales,
\code{\link{predictSurrogate}}, \code{\link{proposeBatch}},
\code{\link{recommend}}.
}
\examples{
surface <- demoSurface()
study <- campaign(surface$design, response = "Ion_response")
study <- proposeSpaceFilling(study, batch = c("B1", "B2"), q = 8, seed = 1)
for (b in c("B1", "B2")) {
  sheet <- proposal(study, b)
  sheet$Ion_response <- simulateBatch(surface, sheet)
  study <- addBatch(study, sheet)
}

study <- runSurrogate(study, "B2_fit", seed = 1)
fit <- surrogate(study)
fit
summary(fit)

## A batch-naive fit for comparison
study <- runSurrogate(study, "B2_naive", naive = TRUE, seed = 1)
surrogateNames(study)
}
