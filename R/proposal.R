################################################################################
# Proposal.R -- campaign-facing acquisition and proposed-batch storage
################################################################################


#================== Acquisition Incumbents ==================

#Runs that count as tested conditions: every Sample run, and every reference run
#that lies inside the experimental domain. Out-of-domain references continue to
#condition the fitted GP, but cannot define an attainable incumbent or best.
#Shared by the acquisition incumbent and cumulativeBest().
.eligible_runs <- function(coded, type, include_references = TRUE, tolerance = 1e-8) {
  inside_domain <- apply(
    coded >= -tolerance & coded <= 1 + tolerance,
    MARGIN = 1L,
    FUN = all
  )

  type == "Sample" | (include_references & type == "Ref" & inside_domain)
}


#The acquisition incumbent uses every distinct eligible completed condition.
.acquisition_baseline <- function(x) {
  eligible <- .eligible_runs(x@coded, x@observations$Type)
  unique(x@coded[eligible, , drop = FALSE])
}


#================== Proposal Assembly ==================

#Append the declared number of fixed-reference runs to physical Sample
#conditions. A batch can still be relabelled by editing the Batch column of its
#sheet before addBatch().
.assemble_proposed_batch <- function(x, sample_conditions, batch = NULL) {
  sample_conditions <- as.data.frame(sample_conditions, stringsAsFactors = FALSE)

  sample_batch <- cbind(
    data.frame(
      Type = rep("Sample", nrow(sample_conditions)),
      stringsAsFactors = FALSE
    ),
    sample_conditions[, x@design@label, drop = FALSE]
  )

  if (!is.null(batch)) {
    sample_batch$Batch <- rep(batch, nrow(sample_batch))
    sample_batch <- sample_batch[, c("Batch", setdiff(names(sample_batch), "Batch"))]
  }

  if (!nrow(x@controls)) {
    rownames(sample_batch) <- NULL
    return(sample_batch)
  }

  control_index <- rep(
    seq_len(nrow(x@controls)),
    times = x@control_replicates
  )

  reference_conditions <- x@controls[
    control_index,
    x@design@label,
    drop = FALSE
  ]

  reference_batch <- cbind(
    data.frame(
      Type = rep("Ref", nrow(reference_conditions)),
      stringsAsFactors = FALSE
    ),
    reference_conditions
  )

  if (!is.null(batch)) {
    reference_batch$Batch <- rep(batch, nrow(reference_batch))
    reference_batch <- reference_batch[
      ,
      c("Batch", setdiff(names(reference_batch), "Batch")),
      drop = FALSE
    ]
  }

  batch <- rbind(sample_batch, reference_batch)
  rownames(batch) <- NULL
  batch
}


#Collect optimizer diagnostics attached to the coded proposal matrix.
.proposal_diagnostics <- function(coded_proposal) {
  list(
    value = attr(coded_proposal, "value"),
    convergence = attr(coded_proposal, "convergence"),
    message = attr(coded_proposal, "message"),
    screen_max = attr(coded_proposal, "screen_max"),
    screen_attempts = attr(coded_proposal, "screen_attempts"),
    screen_size = attr(coded_proposal, "screen_size"),
    separation = attr(coded_proposal, "separation")
  )
}


#================== Campaign Proposal ==================

#proposeSpaceFilling() generates the initial conditions from one global Latin
#hypercube or randomized Sobol design, distributed across the requested batches.
#proposeBatch() uses a stored surrogate and a batched acquisition function
#(qnei, kb or ts) to propose the next batch; acquisition controls are passed
#through ... Each proposal is stored under its physical batch label, which is
#also written to the run sheet, and is never silently replaced. proposal()
#retrieves one or more stored sheets by label (details = TRUE for the full entry),
#and writes them to a CSV file when a file name is given.

#' @export
proposeSpaceFilling <- function(
    x,
    batch = NULL,
    n_batches = NULL,
    q,
    method = c("lhs", "sobol"),
    seed = NULL,
    randomize = FALSE,
    reference_order = NULL,
    overwrite = FALSE
) {
  if (!methods::is(x, "BOCampaign")) {
    stop('"x" must be a BOCampaign.', call. = FALSE)
  }

  if (is.null(n_batches)) {
    n_batches <- if (is.null(batch)) 1L else length(batch)
  }
  n_batches <- .validate_count(n_batches, "n_batches")

  #Batch labels: one per batch, or initial_01, initial_02, ...
  if (is.null(batch)) {
    batch <- paste0("initial_", sprintf("%02d", seq_len(n_batches)))
  } else {
    if (length(batch) != n_batches) {
      stop('"batch" must contain one label per batch (', n_batches, ').', call. = FALSE)
    }
    batch <- vapply(batch, .validate_name, character(1), name = "batch", USE.NAMES = FALSE)
    if (anyDuplicated(batch)) {
      stop('"batch" labels must be unique.', call. = FALSE)
    }
  }
  if (nrow(x@observations) && any(batch %in% x@observations$Batch)) {
    stop('Batch label(s) already present in the completed observations: ',
         paste(intersect(batch, x@observations$Batch), collapse = ", "), call. = FALSE)
  }
  for (label in batch) {
    overwrite <- .validate_overwrite_name(label, names(x@proposals), overwrite, what = "proposal")
  }

  q <- .validate_count(q, "q")
  seed <- .validate_seed(seed)
  randomize <- .validate_flag(randomize, "randomize")
  reference_order <- .resolve_reference_order(reference_order, randomize)
  method <- match.arg(method)

  design <- .space_filling_samples(
    x = x,
    n_batches = n_batches,
    q = q,
    method = method,
    seed = seed
  )

  run_order_seeds <- vapply(seq_len(n_batches), function(batch_index) {
    value <- .seed_offset(seed, 1000L + batch_index)
    if (is.null(value)) NA_integer_ else value
  }, integer(1))

  sheets <- lapply(seq_len(n_batches), function(batch_index) {
    batch_label <- batch[batch_index]
    physical_conditions <- .assemble_proposed_batch(
      x = x,
      sample_conditions = design$batches[[batch_index]]
        [, x@design@label, drop = FALSE],
      batch = batch_label
    )

    .prepare_run_sheet(
      conditions = physical_conditions,
      response = x@response,
      randomize = randomize,
      reference_order = reference_order,
      seed = if (is.na(run_order_seeds[batch_index])) {
        NULL
      } else {
        run_order_seeds[batch_index]
      }
    )
  })
  sheets <- unname(sheets)

  settings <- list(
    design_method = method,
    n_batches = n_batches,
    batches = batch,
    q = q,
    seed = seed,
    randomize = randomize,
    reference_order = reference_order,
    run_order_seed = if (is.null(seed)) NULL else run_order_seeds
  )

  #Store each batch under its own label; they share the settings of the global design.
  for (batch_index in seq_len(n_batches)) {
    coded <- design$coded[(batch_index - 1) * q + seq_len(q), , drop = FALSE]

    provenance <- .provenance_space_filling_proposal(x = x, settings = settings)

    x@proposals[[batch[batch_index]]] <- list(
      batch = sheets[[batch_index]],
      batch_label = batch[batch_index],
      coded = coded,
      surrogate = NULL,
      method = "space_filling",
      settings = settings,
      diagnostics = NULL,
      provenance = provenance,
      call = match.call()
    )
  }

  methods::validObject(x)
  x
}


#' @export
proposeBatch <- function(
    x,
    batch,
    surrogate_name = NULL,
    method = c("qnei", "kb", "ts"),
    q = 5L,
    seed = NULL,
    randomize = FALSE,
    reference_order = NULL,
    overwrite = FALSE,
    ...
) {
  if (!methods::is(x, "BOCampaign")) {
    stop('"x" must be a BOCampaign.', call. = FALSE)
  }

  if (missing(batch)) {
    stop('"batch" is required for an adaptive proposal.', call. = FALSE)
  }
  batch <- .validate_name(batch, "batch")
  #The proposal is stored under its batch label.
  overwrite <- .validate_overwrite_name(
    batch,
    names(x@proposals),
    overwrite,
    what = "proposal"
  )
  surrogate_name <- .resolve_surrogate_name(x, surrogate_name)
  q <- .validate_count(q, "q")
  seed <- .validate_seed(seed)
  randomize <- .validate_flag(randomize, "randomize")
  reference_order <- .resolve_reference_order(reference_order, randomize)

  if (
    nrow(x@observations) &&
    batch %in% unique(as.character(x@observations$Batch))
  ) {
    stop(
      'Batch label "',
      batch,
      '" is already present in the completed observations.',
      call. = FALSE
    )
  }

  method <- match.arg(method)
  surrogate_status <- .provenance_status(
    x@surrogates[[surrogate_name]]$provenance,
    x
  )
  if (!surrogate_status$current) {
    source_batches <- x@surrogates[[surrogate_name]]$provenance$n_batches
    warning(
      'Surrogate "', surrogate_name, '" is historical and was fitted using ',
      source_batches, " batch(es); the campaign currently contains ",
      nbatches(x), ".",
      call. = FALSE
    )
  }
  fit <- surrogate(
    x,
    name = surrogate_name,
    strict = FALSE
  )

  sample_rows <- x@observations$Type == "Sample"
  if (!any(sample_rows)) {
    stop(
      "At least one completed Sample condition is required before proposing.",
      call. = FALSE
    )
  }

  X_baseline <- .acquisition_baseline(x)

  settings <- list(...)

  if (length(settings)) {
    if (is.null(names(settings)) || any(!nzchar(names(settings)))) {
      stop(
        "Additional acquisition settings must be supplied as named arguments.",
        call. = FALSE
      )
    }

    reserved <- intersect(
      names(settings),
      c("g", "q", "seed", "X_baseline")
    )

    if (length(reserved)) {
      stop(
        "Do not supply reserved acquisition argument(s) through ...: ",
        paste(reserved, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  }

  acquisition_arguments <- c(
    list(
      g = fit,
      q = q,
      seed = seed
    ),
    settings
  )

  #KB and qNEI also use in-domain references as incumbents. Out-of-domain
  #references remain in the GP likelihood. Thompson sampling has no incumbent.
  if (method %in% c("kb", "qnei")) {
    acquisition_arguments$X_baseline <- X_baseline
  }

  acquisition <- .acquisition_function(method)
  coded_proposal <- do.call(
    acquisition,
    acquisition_arguments
  )

  coded_proposal <- .validate_proposal_matrix(
    coded_proposal,
    q = q,
    D = length(x@design@label)
  )

  colnames(coded_proposal) <- x@design@label
  diagnostics <- .proposal_diagnostics(coded_proposal)
  sample_conditions <- .decode(x@design, coded_proposal)
  physical_conditions <- .assemble_proposed_batch(
    x = x,
    sample_conditions = sample_conditions,
    batch = batch
  )

  run_order_seed <- .seed_offset(seed, 1001L)
  proposed_batch <- .prepare_run_sheet(
    conditions = physical_conditions,
    response = x@response,
    randomize = randomize,
    reference_order = reference_order,
    seed = run_order_seed
  )

  stored_settings <- c(
    list(
      batch = batch,
      method = method,
      q = q,
      seed = seed,
      randomize = randomize,
      reference_order = reference_order,
      run_order_seed = run_order_seed
    ),
    settings
  )

  provenance <- .provenance_proposal(
    x = x,
    surrogate_name = surrogate_name,
    method = method,
    settings = stored_settings
  )

  x@proposals[[batch]] <- list(
    batch = proposed_batch,
    batch_label = batch,
    coded = coded_proposal,
    surrogate = surrogate_name,
    method = method,
    settings = stored_settings,
    diagnostics = diagnostics,
    provenance = provenance,
    call = match.call()
  )

  methods::validObject(x)
  x
}


#================== Proposal Accessors ==================

#Hide an entirely unmeasured response from the ordinary proposal view. Multiple
#initial sheets are combined into one consistently shaped data frame; their
#separate stored sheets remain available through proposal(..., details = TRUE).
.proposal_display <- function(batch, response) {
  if (is.list(batch) && !is.data.frame(batch)) {
    displayed <- lapply(batch, .proposal_display, response = response)
    combined <- do.call(rbind, displayed)
    rownames(combined) <- NULL
    return(combined)
  }

  if (
    is.data.frame(batch) &&
    response %in% names(batch) &&
    all(is.na(batch[[response]]))
  ) {
    batch[[response]] <- NULL
  }

  batch
}

#' @export
proposal <- function(
    x,
    batch,
    file = NULL,
    strict = FALSE,
    details = FALSE
) {
  if (!methods::is(x, "BOCampaign")) {
    stop('"x" must be a BOCampaign.', call. = FALSE)
  }

  strict <- .validate_flag(strict, "strict")
  details <- .validate_flag(details, "details")

  available <- if (length(x@proposals)) paste(names(x@proposals), collapse = ", ") else "none"
  if (missing(batch)) {
    stop('Specify the batch label(s) to retrieve. Available: ', available, ".", call. = FALSE)
  }
  batch <- unique(vapply(batch, .validate_name, character(1), name = "batch", USE.NAMES = FALSE))

  unknown <- setdiff(batch, names(x@proposals))
  if (length(unknown)) {
    stop(
      'No proposal for batch: ',
      paste(unknown, collapse = ", "),
      '. Available: ',
      available,
      ".",
      call. = FALSE
    )
  }

  entries <- x@proposals[batch]

  for (label in batch) {
    .check_current(
      entries[[label]]$provenance,
      x,
      paste0('Proposal "', label, '"'),
      strict = strict
    )
  }

  #Write the requested sheets, empty response column included, to one CSV file.
  if (!is.null(file)) {
    sheet <- do.call(rbind, unname(lapply(entries, `[[`, "batch")))
    rownames(sheet) <- NULL
    .write_run_sheet(sheet, file)
  }

  if (details) {
    if (length(entries) == 1L) entries[[1L]] else entries
  } else {
    .proposal_display(lapply(entries, `[[`, "batch"), x@response)
  }
}


#' @export
proposalNames <- function(x) names(x@proposals)
