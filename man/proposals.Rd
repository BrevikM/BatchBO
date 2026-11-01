\name{proposals}
\alias{proposals}
\alias{proposeSpaceFilling}
\alias{proposeBatch}
\alias{proposal}
\alias{proposalNames}
\title{Propose the next batches of experiments}
\description{
\code{proposeSpaceFilling()} proposes the initial batches from a space-filling
design. \code{proposeBatch()} uses a fitted surrogate and a batch acquisition
function to propose the next batch. Both store the run sheet in the campaign
under its batch label; \code{proposal()} retrieves it, ready to be run, filled
in and passed to \code{\link{addBatch}}.
}
\usage{
proposeSpaceFilling(x, batch = NULL, n_batches = NULL, q,
                    method = c("lhs", "sobol"), seed = NULL,
                    randomize = FALSE, reference_order = NULL,
                    overwrite = FALSE)
proposeBatch(x, batch, surrogate_name = NULL, method = c("qnei", "kb", "ts"),
             q = 5L, seed = NULL, randomize = FALSE, reference_order = NULL,
             overwrite = FALSE, ...)
proposal(x, batch, file = NULL, strict = FALSE, details = FALSE)
proposalNames(x)
}
\arguments{
  \item{x}{A \code{\linkS4class{BOCampaign}}.}
  \item{batch}{Batch label(s). The label is written to the \code{Batch} column
    of the run sheet and is also the name under which the proposal is stored;
    it must not already be used by a completed batch. For
    \code{proposeSpaceFilling()}, one unique label per batch, or \code{NULL}
    for \code{"initial_01"}, \code{"initial_02"}, \ldots. Required for
    \code{proposeBatch()}. For \code{proposal()}, one or more stored labels.}
  \item{n_batches}{Number of initial batches. \code{NULL} (default) uses
    \code{length(batch)}, or 1 when \code{batch} is \code{NULL}.}
  \item{q}{Number of Sample conditions per batch. Reference runs are added on
    top (see \code{\link{campaign}}).}
  \item{method}{For \code{proposeSpaceFilling()}, \code{"lhs"} (default) for a
    Latin hypercube or \code{"sobol"} for a randomized Sobol' sequence. For
    \code{proposeBatch()}, \code{"qnei"} (default), \code{"kb"} or \code{"ts"};
    see Details.}
  \item{surrogate_name}{Name of the stored surrogate to use. \code{NULL}
    (default) uses the most recently stored current surrogate. A historical
    surrogate may be named explicitly, with a warning.}
  \item{seed}{\code{NULL} (default) or an integer seed. Seeds for the design,
    the acquisition and the run order are derived from it, and the global
    random-number stream is left untouched.}
  \item{randomize}{Logical; shuffle the order of the Sample runs and place the
    reference runs at random positions within the chosen policy. Default
    \code{FALSE} keeps the proposed order.}
  \item{reference_order}{Where the reference runs go in the run order:
    \code{"spread"} (one in each of equally sized segments of the run),
    \code{"random"} (anywhere; requires \code{randomize = TRUE}),
    \code{"last"} (after all Sample runs) or \code{"bracket"} (first and last,
    any others spread in between). \code{NULL} (default) gives
    \code{"spread"} when randomized and \code{"last"} otherwise.}
  \item{overwrite}{Logical; allow replacing a stored proposal with the same
    label. Default \code{FALSE}.}
  \item{...}{Named settings passed to the acquisition; see Details.}
  \item{file}{\code{NULL} (default), or the path of a CSV file to which
    \code{proposal()} writes the requested run sheet(s), with the empty
    response column ready to be filled in. An existing file is replaced; the
    directory must exist.}
  \item{strict}{Logical; stop if a requested proposal is historical (see
    \code{\link{summary.BOCampaign}}).}
  \item{details}{Logical; return the complete stored entry instead of the run
    sheet.}
}
\details{
\strong{Run sheets.} A run sheet has columns \code{Run} (run order),
\code{Batch}, \code{Type} (\code{"Sample"} or \code{"Ref"}), one column per
factor in physical units, and an empty response column. \code{proposal()}
hides the response column while it is still empty and stacks the sheets when
several labels are requested. With \code{file}, the same sheets are written
to one CSV file, including the response column, so the completed file can be
read back with \code{\link[utils]{read.csv}} and passed to
\code{\link{addBatch}}. The \code{Batch} label may be changed before the
completed sheet is added.

\strong{Initial batches.} \code{proposeSpaceFilling()} draws one design of
\code{n_batches * q} points over the whole design space (on the coded scale,
so log-scaled factors are spread evenly in \eqn{\log_{10}}{log10}) and allocates its
points to the batches at random, so that no batch is tied to one region of
the space. The Latin hypercube places each point uniformly at random within
its cell; the Sobol' design uses a digitally shifted Sobol' sequence from
\pkg{qrng}.

\strong{Adaptive batches.} \code{proposeBatch()} selects the \code{q}
conditions one at a time, each time maximizing an acquisition function that
accounts for the conditions already chosen for the batch. All three methods
work on the batch-free surface \eqn{f(x)}, so batch offsets do not influence
where to measure next. The incumbent (the level to improve upon) is defined by
every distinct Sample condition measured so far, and by reference conditions
that lie inside the bounds.
\describe{
  \item{\code{"qnei"}}{Noisy expected improvement (Letham et al., 2019),
    estimated from \code{S_n} joint posterior draws at the incumbent
    conditions, the conditions already chosen and the candidate. Measurement
    noise in the incumbent is accounted for, and the same draws are reused
    for every candidate, which makes the estimate smooth enough to
    optimize.}
  \item{\code{"kb"}}{Kriging believer (Ginsbourger et al., 2010): expected
    improvement (Jones et al., 1998), where each chosen condition is treated as
    if measured at its predicted value. This reduces the uncertainty near it
    and steers the next choice elsewhere. The fastest of the three.}
  \item{\code{"ts"}}{Thompson sampling: each batch member maximizes its own
    independent draw of \eqn{f(x)} from the posterior, approximated with
    \code{F_n} random Fourier features (Rahimi and Recht, 2007) and
    conditioned on the data by pathwise updating (Wilson et al., 2020). The
    most exploratory of the three; with few features the draws are rough.}
}
Each acquisition is maximized by evaluating it on a randomized Sobol' set of
\code{n_screen} points and running L-BFGS-B from the best \code{restarts} of
them. For qNEI and KB, if no screened point gives a positive value, the
screen is redrawn up to \code{max_screen_factor} times larger.

Settings accepted through \code{...}:
\describe{
  \item{\code{n_screen}}{Size of the Sobol' screen (default 4096).}
  \item{\code{restarts}}{Number of local optimizations (default 12).}
  \item{\code{maxit}}{Maximum L-BFGS-B iterations per start (default 500).}
  \item{\code{lower}, \code{upper}}{Search bounds on the coded scale
    (vectors with one value per factor in \eqn{[0, 1]}; default the whole
    design space).}
  \item{\code{max_screen_factor}}{qNEI and KB only (default 4).}
  \item{\code{S_n}}{qNEI only: number of posterior draws (default 256).}
  \item{\code{F_n}}{TS only: number of random Fourier features (default
    1024).}
}

\strong{Stored entries.} With \code{details = TRUE}, \code{proposal()} returns
the stored entry: the run sheet (\code{batch}), \code{batch_label}, the
proposed Sample conditions on the coded scale (\code{coded}), the
\code{surrogate} name, \code{method}, \code{settings}, optimizer
\code{diagnostics} (acquisition value, convergence code and message, screen
maximum, size and attempts per batch member, and the smallest
length-scale-weighted distance between batch members), \code{provenance} and
the \code{call}.
}
\value{
\code{proposeSpaceFilling()} and \code{proposeBatch()} return the updated
campaign. \code{proposal()} returns a data frame with the run sheet(s), or,
with \code{details = TRUE}, the stored entry (a list of entries if several
labels are requested). \code{proposalNames()} returns the labels of all stored
proposals.
}
\references{
Ginsbourger, D., Le Riche, R. and Carraro, L. (2010). Kriging is well-suited to
parallelize optimization. In Tenne, Y. and Goh, C.-K. (eds.),
\emph{Computational Intelligence in Expensive Optimization Problems},
pp. 131--162. Springer.

Jones, D. R., Schonlau, M. and Welch, W. J. (1998). Efficient global
optimization of expensive black-box functions. \emph{Journal of Global
Optimization}, \bold{13}, 455--492.

Letham, B., Karrer, B., Ottoni, G. and Bakshy, E. (2019). Constrained Bayesian
optimization with noisy experiments. \emph{Bayesian Analysis}, \bold{14}(2),
495--519.

McKay, M. D., Beckman, R. J. and Conover, W. J. (1979). A comparison of three
methods for selecting values of input variables in the analysis of output from
a computer code. \emph{Technometrics}, \bold{21}(2), 239--245.

Rahimi, A. and Recht, B. (2007). Random features for large-scale kernel
machines. \emph{Advances in Neural Information Processing Systems}, \bold{20}.

Wilson, J. T., Borovitskiy, V., Terenin, A., Mostowsky, P. and Deisenroth, M. P.
(2020). Efficiently sampling functions from Gaussian process posteriors.
\emph{Proceedings of the 37th International Conference on Machine Learning}.
}
\seealso{
\code{\link{runSurrogate}} to fit the surrogate, \code{\link{addBatch}} to add
the completed run sheet, \code{\link{recommend}} for the final
recommendation.
}
\examples{
surface <- demoSurface()
reference <- data.frame(Capillary_voltage = 3500, Nebulizer_pressure = 35,
                        Nozzle_voltage = 1000, Drying_gas_flow = 10,
                        Drying_gas_temperature = 275, Sheath_gas_flow = 10)
study <- campaign(surface$design, response = "Ion_response",
                  controls = reference, replicates = 2)

## Two initial batches of six conditions, references spread through the run
study <- proposeSpaceFilling(study, batch = c("B1", "B2"), q = 6,
                             randomize = TRUE, seed = 1)
proposalNames(study)
proposal(study, "B1")

## Write both sheets to one file for the lab
sheet_file <- file.path(tempdir(), "initial_batches.csv")
proposal(study, c("B1", "B2"), file = sheet_file)
head(read.csv(sheet_file))

## Measure them (here simulated) and add them
for (b in c("B1", "B2")) {
  sheet <- proposal(study, b)
  sheet$Ion_response <- simulateBatch(surface, sheet)
  study <- addBatch(study, sheet)
}

\donttest{
## Fit a surrogate and propose the next batch
study <- runSurrogate(study, "B2_fit", seed = 1)
study <- proposeBatch(study, batch = "B3", method = "kb", q = 6, seed = 1)
proposal(study, "B3")
proposal(study, "B3", details = TRUE)$diagnostics

## Restrict the search to the lower half of the first factor
study <- proposeBatch(study, batch = "B3", method = "qnei", q = 6, seed = 1,
                      upper = c(0.5, 1, 1, 1, 1, 1), overwrite = TRUE)
}
}
