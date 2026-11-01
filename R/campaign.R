################################################################################
# Campaign.R -- construction, growth, standardisation, and surrogate storage
################################################################################


#================== Campaign Construction ==================

#campaign() creates the empty, resumable object that carries an optimization
#campaign. addBatch() appends or replaces completed batches, and removeBatch()
#removes named batches. #The response is transformed and standardized using 
#completed sample rows (reference runs enter the GP likelihood but not these constants).

#' @export
campaign <- function(
    design,
    response,
    transform = c("identity", "log"),
    direction = c("maximize", "minimize"),
    controls = NULL,
    replicates = 2
) {
  if (!methods::is(design, "DesignSpace")) {
    stop('"design" must be a DesignSpace.', call. = FALSE)
  }

  response <- .validate_response_name(response, design@label)

  transform <- match.arg(transform)
  direction <- match.arg(direction)

  if (is.null(controls) || !nrow(as.data.frame(controls))) {
    controls <- data.frame()
    replicates <- integer(0)
  } else {
    controls <- as.data.frame(controls, stringsAsFactors = FALSE)

    if (nrow(controls) != 1) {
      stop(
        '"controls" must contain exactly one fixed reference condition.',
        call. = FALSE
      )
    }

    controls <- .validate_numeric_columns(
      controls,
      design@label,
      '"controls"'
    )
    .validate_log_domain(controls, design, '"controls"')
    replicates <- .validate_replicates(replicates, nrow(controls))
  }

  object <- methods::new( #Initiate an object from a class
    "BOCampaign",
    design = design,
    response = response,
    response_transform = transform,
    direction = direction,
    controls = controls,
    control_replicates = replicates,
    observations = data.frame(),
    coded = matrix( #Define a matrix with coded variables (initially zero data)
      numeric(0), #Define it as a matrix for ease in subsequent matrix algebra
      nrow = 0,
      ncol = length(design@label),
      dimnames = list(NULL, design@label)
    ),
    standardisation = list() #Filled in by .standardise() below
  )
  object@standardisation <- .standardise(object)

  methods::validObject(object)
  object
}


#================== Response Handling ==================

#Transform a response before standardisation.
.transform_response <- function(x, y) { #For some campaign x, and data y
  switch(
    x@response_transform,
    identity = y, #If identity (default), keep y as is
    log = log(y) #If log, do a natural-log transform
  )
}


#Calculate response-standardisation constants from Sample rows only.
.standardise <- function(x) {
  sign <- if (x@direction == "maximize") 1 else -1
  sample_rows <- x@observations$Type == "Sample" 
  n_sample <- sum(sample_rows)

  if (n_sample < 2) {
    return(
      list(
        y_mean = NA,
        y_sd = NA,
        sign = sign,
        n_sample = n_sample
      )
    )
  }

  transformed <- .transform_response(
    x,
    x@observations[[x@response]]
  )

  y_sd <- stats::sd(transformed[sample_rows])

  if (!is.finite(y_sd) || y_sd <= 0) {
    stop(
      "Standardization failed. The data has zero variance.",
      call. = FALSE
    )
  }

  list(
    y_mean = mean(transformed[sample_rows]),
    y_sd = y_sd,
    sign = sign,
    n_sample = n_sample
  )
}


#Return the transformed and standardised response, with larger values always
#representing better conditions for the acquisition functions.
.response_coded <- function(x) {
  n <- nrow(x@observations)
  if (!n) {
    return(numeric(0))
  }

  standardisation <- x@standardisation

  if (
    !is.finite(standardisation$y_mean) ||
    !is.finite(standardisation$y_sd)
  ) {
    return(rep(NA, n))
  }

  transformed <- .transform_response(
    x,
    x@observations[[x@response]]
  )

  standardisation$sign*(transformed - standardisation$y_mean)/standardisation$y_sd
}


#Return a coded response to its original scale
.decode_response <- function(x, z, standardisation = x@standardisation) {
  transformed <- standardisation$y_mean+standardisation$sign*z*standardisation$y_sd

  switch(
    x@response_transform,
    identity = transformed,
    log = exp(transformed)
  )
}


#Decode with a fit's own standardisation, which is the right scale also for a
#historical surrogate. For a log response this exponentiates the latent mean
#without a variance correction.
.decode_fit_response <- function(x, fit, value) {
  .decode_response(
    x,
    value,
    standardisation = list(y_mean = fit$y_mean, y_sd = fit$y_sd, sign = fit$sign)
  )
}


#Rebuild all state derived directly from the observations.
.rebuild <- function(x) {
  x@coded <- if (nrow(x@observations)) {
    .encode(x@design, x@observations)
  } else {
    matrix(
      numeric(0),
      nrow = 0,
      ncol = length(x@design@label),
      dimnames = list(NULL, x@design@label)
    )
  }

  x@standardisation <- .standardise(x)
  x
}


#================== Campaign Growth ==================

#' @export
addBatch <- function(x, data, overwrite = FALSE) {
  if (!methods::is(x, "BOCampaign")) {
    stop('"x" must be a BOCampaign.', call. = FALSE)
  }

  overwrite <- .validate_flag(overwrite, "overwrite")

  #Accept one data frame, or a list of data frames (for example one per batch).
  if (is.list(data) && !is.data.frame(data)) {
    if (!length(data) || !all(vapply(data, is.data.frame, logical(1)))) {
      stop('"data" must be a data frame or a list of data frames.', call. = FALSE)
    }
    if (length(unique(lapply(data, function(d) sort(names(d))))) > 1) {
      stop('All data frames in "data" must have the same columns.', call. = FALSE)
    }
    data <- do.call(rbind, lapply(data, function(d) d[, names(data[[1]]), drop = FALSE]))
  }

  incoming <- as.data.frame(data, stringsAsFactors = FALSE)

  if (!nrow(incoming)) {
    stop('"data" must contain at least one completed run.', call. = FALSE)
  }

  #Only the columns of a completed run sheet are accepted: Batch, Type, the
  #factors and the response are required, and Run is optional.
  required <- c("Batch", "Type", x@design@label, x@response)
  missing_columns <- setdiff(required, names(incoming))
  if (length(missing_columns)) {
    stop(
      '"data" is missing required column(s): ',
      paste(missing_columns, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  unexpected <- setdiff(names(incoming), c("Run", required))
  if (length(unexpected)) {
    stop(
      '"data" contains unexpected column(s): ',
      paste(unexpected, collapse = ", "),
      ". Only Batch, Type, Run (optional), the factors and the response ",
      "are accepted.",
      call. = FALSE
    )
  }

  incoming <- incoming[
    ,
    intersect(c("Run", required), names(incoming)),
    drop = FALSE
  ]

  incoming$Batch <- trimws(as.character(incoming$Batch))

  if (anyNA(incoming$Batch) || any(!nzchar(incoming$Batch))) {
    stop("Batch labels must be non-missing and non-empty.", call. = FALSE)
  }

  type <- tolower(trimws(as.character(incoming$Type)))
  incoming$Type <- ifelse(
    type %in% c("ref", "reference"),
    "Ref",
    ifelse(type == "sample", "Sample", NA_character_)
  )

  if (anyNA(incoming$Type)) {
    stop(
      '"Type" must contain only "Sample", "Ref", or "Reference".',
      call. = FALSE
    )
  }

  incoming <- .validate_numeric_columns(
    incoming,
    c(x@design@label, x@response),
    "Incoming batch"
  )
  .validate_log_domain(incoming, x@design, "Incoming factor values")
  .validate_response_values(
    incoming[[x@response]],
    x@response_transform,
    "The incoming response"
  )

  #Run, when supplied, must be positive whole numbers, unique within each batch.
  if ("Run" %in% names(incoming)) {
    incoming <- .validate_numeric_columns(incoming, "Run", "Incoming batch")
    run_order <- incoming$Run

    if (any(run_order < 1 | run_order != floor(run_order))) {
      stop('"Run" must contain positive whole numbers.', call. = FALSE)
    }
    if (anyDuplicated(data.frame(incoming$Batch, run_order))) {
      stop('"Run" must be unique within each batch.', call. = FALSE)
    }
  }

  #Reference runs must use the campaign's declared reference condition.
  reference_rows <- which(incoming$Type == "Ref")

  if (length(reference_rows)) {
    if (!nrow(x@controls)) {
      stop(
        '"data" contains reference runs, but the campaign declares no reference condition.',
        call. = FALSE
      )
    }

    reference <- as.numeric(x@controls[1, x@design@label])
    supplied <- as.matrix(incoming[reference_rows, x@design@label, drop = FALSE])
    differs <- abs(sweep(supplied, 2, reference, "-")) >
      1e-6 * pmax(1, abs(reference))

    if (any(differs)) {
      stop(
        "Reference runs must use the declared reference condition; they differ in: ",
        paste(x@design@label[apply(differs, 2, any)], collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  }

  #Batches already in the campaign: an identical re-supplied batch is skipped,
  #while a batch with different data is refused unless overwrite = TRUE.
  registered <- intersect(
    unique(incoming$Batch),
    unique(x@observations$Batch)
  )
  is_identical <- vapply(registered, function(label) {
    .same_batch(
      x@observations[x@observations$Batch == label, , drop = FALSE],
      incoming[incoming$Batch == label, , drop = FALSE],
      c(x@design@label, x@response)
    )
  }, logical(1))
  identical_batches <- registered[is_identical]
  duplicate_batches <- registered[!is_identical]

  if (length(duplicate_batches) && !overwrite) {
    stop(
      "Batch label(s) already registered with different data: ",
      paste(duplicate_batches, collapse = ", "),
      '. Set "overwrite = TRUE" to replace them.',
      call. = FALSE
    )
  }

  if (length(identical_batches)) {
    message(
      "Batch(es) already registered with identical data, skipped: ",
      paste(identical_batches, collapse = ", "),
      "."
    )
    incoming <- incoming[!incoming$Batch %in% identical_batches, , drop = FALSE]

    if (!nrow(incoming)) {
      return(x)
    }
  }

  if (nrow(x@observations)) {
    #Existing and incoming data can differ only in whether Run is present.
    all_columns <- union(names(x@observations), names(incoming))

    align_columns <- function(object) {
      for (column in setdiff(all_columns, names(object))) {
        object[[column]] <- NA
      }
      object[, all_columns, drop = FALSE]
    }

    existing <- align_columns(x@observations)
    incoming <- align_columns(incoming)

    if (length(duplicate_batches)) {
      #Keep each existing batch in its original chronological position.
      #Incoming batches not previously observed are appended at the end.
      existing_order <- unique(as.character(existing$Batch))
      incoming_order <- unique(as.character(incoming$Batch))
      batch_order <- c(
        existing_order,
        setdiff(incoming_order, existing_order)
      )

      x@observations <- do.call(rbind, lapply(batch_order, function(label) {
        if (label %in% incoming_order) {
          incoming[incoming$Batch == label, , drop = FALSE]
        } else {
          existing[existing$Batch == label, , drop = FALSE]
        }
      }))
    } else {
      x@observations <- rbind(existing, incoming)
    }
  } else {
    x@observations <- incoming
  }

  rownames(x@observations) <- NULL

  x <- .rebuild(x)
  methods::validObject(x)

  if (nrow(x@controls)) {
    without_reference <- setdiff(
      unique(incoming$Batch),
      unique(incoming$Batch[incoming$Type == "Ref"])
    )
    if (length(without_reference)) {
      warning(
        "The campaign declares a reference condition, but these batches ",
        "contain no reference runs: ",
        paste(without_reference, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  }

  if (any(incoming$Type == "Ref")) {
    reference_coded <- .encode(x@design, x@controls)
    outside <- reference_coded < -1e-8 | reference_coded > 1 + 1e-8

    if (any(outside)) {
      message(
        "Reference condition lies outside the declared bounds for: ",
        paste(x@design@label[outside], collapse = ", "),
        ". This is allowed, but the acquisition cannot propose there."
      )
    }
  }

  x
}


#Do two sets of runs contain the same runs (Type, factor settings and response),
#in any order? Numeric values are compared to a small relative tolerance, so a
#batch written to CSV and read back still counts as identical.
.same_batch <- function(a, b, columns) {
  if (nrow(a) != nrow(b)) {
    return(FALSE)
  }

  sorted <- function(d) {
    d <- d[, c("Type", columns), drop = FALSE]
    d[do.call(base::order, unname(as.list(d))), , drop = FALSE]
  }
  a <- sorted(a)
  b <- sorted(b)

  if (!identical(as.character(a$Type), as.character(b$Type))) {
    return(FALSE)
  }

  A <- as.matrix(a[columns])
  B <- as.matrix(b[columns])
  all(abs(A - B) <= 1e-8 * pmax(1, abs(A)))
}


#' @export
removeBatch <- function(x, batches) {
  if (!methods::is(x, "BOCampaign")) {
    stop('"x" must be a BOCampaign.', call. = FALSE)
  }
  if (
    !is.character(batches) ||
    !length(batches) ||
    anyNA(batches) ||
    any(!nzchar(trimws(batches)))
  ) {
    stop('"batches" must contain one or more non-empty batch labels.', call. = FALSE)
  }
  if (!nrow(x@observations)) {
    stop("The campaign contains no completed batches.", call. = FALSE)
  }

  batches <- unique(trimws(batches))
  available <- unique(as.character(x@observations$Batch))
  missing <- setdiff(batches, available)
  if (length(missing)) {
    stop(
      "Unknown batch label(s): ",
      paste(missing, collapse = ", "),
      ". Available: ",
      paste(available, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  x@observations <- x@observations[
    !x@observations$Batch %in% batches,
    ,
    drop = FALSE
  ]
  rownames(x@observations) <- NULL
  x <- .rebuild(x)
  methods::validObject(x)
  x
}


#================== Surrogate Storage ==================

#runSurrogate() fits a batch-aware GP by maximum marginal likelihood and stores it
#under a name; naive = TRUE fits a batch-naive comparison model. Factors are fitted
#on the coded unit cube with one length scale per factor, and the response on its
#standardized scale. Additional ML-II controls passed through ... are hp_start, an
#initial hyperparameter list, and lower/upper bounds ordered as length scale,
#signal, batch and residual variance (defaults c(0.05, 1e-4, 1e-8, 1e-8) and
#c(50, 5, 5, 5)). surrogate() retrieves a stored fit, by default the most recent
#current one.

#' @export
runSurrogate <- function(
    x,
    name,
    kernel = c("matern52", "matern32", "rbf"),
    naive = FALSE,
    n_starts = 8,
    maxit = 200L,
    jitter = 1e-8,
    seed = NULL,
    overwrite = FALSE,
    ...
) {
  if (!methods::is(x, "BOCampaign")) {
    stop('"x" must be a BOCampaign.', call. = FALSE)
  }

  if (missing(name)) {
    stop('"name" is required when fitting a surrogate.', call. = FALSE)
  }
  name <- .validate_name(name)
  overwrite <- .validate_overwrite_name(
    name,
    names(x@surrogates),
    overwrite,
    what = "surrogate"
  )
  naive <- .validate_flag(naive, "naive")
  seed <- .validate_seed(seed)
  kernel <- match.arg(kernel)

  if (nrow(x@observations) < 3) {
    stop("At least three completed runs are required before fitting.", call. = FALSE)
  }

  if (
    x@standardisation$n_sample < 2 ||
    !is.finite(x@standardisation$y_mean) ||
    !is.finite(x@standardisation$y_sd)
  ) {
    stop(
      "At least two Sample responses with non-zero variance are required before fitting.",
      call. = FALSE
    )
  }

  if (!naive && nbatches(x) < 2) {
    warning(
      "A batch-aware surrogate has only one observed batch; the batch and ",
      "residual variance components are weakly identified.",
      call. = FALSE
    )
  }

  y <- .response_coded(x)

  hyperparameters <- .with_seed(
    seed,
    gp_fit_mlii(
      X = x@coded,
      y = y,
      batch = x@observations$Batch,
      kernel = kernel,
      naive = naive,
      n_starts = n_starts,
      jitter = jitter,
      maxit = maxit,
      ...
    )
  )

  fit <- gp(
    X = x@coded,
    y = y,
    batch = x@observations$Batch,
    hp = hyperparameters,
    jitter = jitter
  )

  #Retain the transformation constants belonging to this particular fit.
  fit$y_mean <- x@standardisation$y_mean
  fit$y_sd <- x@standardisation$y_sd
  fit$sign <- x@standardisation$sign
  fit$naive <- naive
  fit$response <- x@response
  fit$response_transform <- x@response_transform
  fit$direction <- x@direction

  x@surrogates[[name]] <- list(
    fit = fit,
    provenance = .provenance_surrogate(x, fit),
    call = match.call()
  )

  methods::validObject(x)
  x
}


#' @export
surrogate <- function(
    x,
    name = NULL,
    strict = FALSE
) {
  if (!methods::is(x, "BOCampaign")) {
    stop('"x" must be a BOCampaign.', call. = FALSE)
  }

  strict <- .validate_flag(strict, "strict")
  name <- .resolve_surrogate_name(x, name)

  entry <- x@surrogates[[name]]

  .check_current(
    entry$provenance,
    x,
    paste0('Surrogate "', name, '"'),
    strict = strict
  )

  entry$fit
}


#Resolve an explicit surrogate name, or select the most recently stored fit
#that still matches the campaign exactly.
.resolve_surrogate_name <- function(x, name = NULL) {
  available <- names(x@surrogates)
  if (!length(available)) {
    stop("The campaign contains no stored surrogates.", call. = FALSE)
  }

  if (is.null(name)) {
    current <- available[vapply(
      x@surrogates,
      function(entry) .is_current(entry$provenance, x),
      logical(1)
    )]
    if (!length(current)) {
      stop(
        "The campaign contains no current surrogate. Fit a new surrogate or ",
        "select a historical fit explicitly. Available: ",
        paste(available, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    return(tail(current, 1L))
  }

  name <- .validate_name(name)
  if (!name %in% available) {
    stop(
      'No surrogate named "', name, '". Available: ',
      paste(available, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  name
}


#================== Campaign Accessors ==================

#Accessors for the design, completed observations, coded design matrix, numbers of
#runs and batches, stored surrogate names, and the best observed sample response.

#' @export
design <- function(x) x@design
#' @export
observations <- function(x) x@observations
#' @export
coded <- function(x) x@coded
#' @export
nruns <- function(x) nrow(x@observations)

#' @export
nbatches <- function(x) {
  if (!nrow(x@observations)) {
    return(0L)
  }
  length(unique(x@observations$Batch))
}

#' @export
surrogateNames <- function(x) names(x@surrogates)


#================== Campaign Summary ==================

.summary_time <- function(provenance) {
  if (inherits(provenance$at, "POSIXt")) provenance$at else as.POSIXct(NA)
}


.summary_surrogates <- function(x) {
  if (!length(x@surrogates)) {
    return(data.frame(
      name = character(), model = character(), kernel = character(),
      batches = integer(), runs = integer(), stored_at = as.POSIXct(character()),
      status = character(), stringsAsFactors = FALSE
    ))
  }

  do.call(rbind, lapply(names(x@surrogates), function(name) {
    entry <- x@surrogates[[name]]
    data.frame(
      name = name,
      model = if (isTRUE(entry$fit$naive)) "naive" else "batch-aware",
      kernel = .gp_kernel_label(entry$fit$hp$kernel),
      batches = entry$provenance$n_batches,
      runs = entry$provenance$n_runs,
      stored_at = .summary_time(entry$provenance),
      status = .provenance_status(entry$provenance, x)$status,
      stringsAsFactors = FALSE
    )
  }))
}


.summary_analyses <- function(x) {
  if (!length(x@analyses)) {
    return(data.frame(
      name = character(), type = character(), surrogate = character(),
      batches = integer(), runs = integer(), stored_at = as.POSIXct(character()),
      status = character(), stringsAsFactors = FALSE
    ))
  }

  do.call(rbind, lapply(names(x@analyses), function(name) {
    entry <- x@analyses[[name]]
    provenance <- entry$provenance
    data.frame(
      name = name,
      type = provenance$analysis,
      surrogate = provenance$surrogate_name,
      batches = provenance$n_batches,
      runs = provenance$n_runs,
      stored_at = .summary_time(provenance),
      status = .provenance_status(provenance, x)$status,
      stringsAsFactors = FALSE
    )
  }))
}


.summary_proposals <- function(x) {
  if (!length(x@proposals)) {
    return(data.frame(
      name = character(), method = character(),
      surrogate = character(),
      batches = integer(), runs = integer(), stored_at = as.POSIXct(character()),
      status = character(), stringsAsFactors = FALSE
    ))
  }

  do.call(rbind, lapply(names(x@proposals), function(name) {
    entry <- x@proposals[[name]]
    provenance <- entry$provenance
    surrogate_name <- provenance$surrogate_name
    if (is.null(surrogate_name)) surrogate_name <- ""
    data.frame(
      name = name,
      method = entry$method,
      surrogate = surrogate_name,
      batches = provenance$n_batches,
      runs = provenance$n_runs,
      stored_at = .summary_time(provenance),
      status = .provenance_status(provenance, x)$status,
      stringsAsFactors = FALSE
    )
  }))
}


#summary() lists the campaign's stored surrogates, analyses and proposals, marking
#each as "current" or "historical".


#' @method summary BOCampaign
#' @export
summary.BOCampaign <- function(object, ...) {
  structure(
    list(
      campaign = list(
        factors = factorNames(object),
        response = object@response,
        transform = object@response_transform,
        direction = object@direction,
        runs = nruns(object),
        batches = nbatches(object)
      ),
      surrogates = .summary_surrogates(object),
      analyses = .summary_analyses(object),
      proposals = .summary_proposals(object)
    ),
    class = "summary.BOCampaign"
  )
}


#' @method print summary.BOCampaign
#' @export
print.summary.BOCampaign <- function(x, ...) {
  campaign <- x$campaign
  cat("BatchBO campaign summary\n")
  cat("  Factors:  ", length(campaign$factors), "\n", sep = "")
  cat("  Response: ", campaign$response, " (", campaign$direction,
      "; ", campaign$transform, ")\n", sep = "")
  cat("  Data:     ", campaign$runs, " runs in ", campaign$batches,
      " batches\n", sep = "")

  print_registry <- function(title, table) {
    cat("\n", title, "\n", sep = "")
    if (!nrow(table)) {
      cat("  None\n")
    } else {
      print(table, row.names = FALSE)
    }
  }

  print_registry("Surrogates", x$surrogates)
  print_registry("Analyses", x$analyses)
  print_registry("Proposals", x$proposals)
  invisible(x)
}


#' @export
bestObserved <- function(x) {
  sample_rows <- which(x@observations$Type == "Sample")

  if (!length(sample_rows)) {
    return(NULL)
  }

  response <- x@observations[[x@response]][sample_rows]
  index <- if (x@direction == "maximize") {
    which.max(response)
  } else {
    which.min(response)
  }

  row <- sample_rows[index]

  list(
    value = response[index],
    row = row,
    batch = x@observations$Batch[row]
  )
}
