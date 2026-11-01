################################################################################
# provenance.R -- campaign-state snapshots and historical-result detection
################################################################################


#================== Campaign State ==================

#Collect the campaign information that determines the fitted surrogate and the
#set of observations eligible to define the acquisition-function incumbent.
.model_state <- function(x) {
  list(
    design = list(
      label = x@design@label,
      lower = x@design@lower,
      upper = x@design@upper,
      transform = x@design@transform
    ),
    response = x@response,
    response_transform = x@response_transform,
    direction = x@direction,
    coded = x@coded,
    response_coded = .response_coded(x),
    batch = as.character(x@observations$Batch),  #character(0) for an empty campaign
    type = as.character(x@observations$Type)
  )
}


#Collect the declared reference conditions separately. Changing these does not
#invalidate an existing surrogate, but it can invalidate a stored proposal.
.control_state <- function(x) {
  list(
    conditions = x@controls,
    replicates = x@control_replicates
  )
}


#Store an exact snapshot of the campaign state. For the intended campaign sizes,
#direct comparison is sufficiently small and avoids an additional hash dependency.
#Each constructor below adds its own kind.
.provenance <- function(x) {
  list(
    state = .model_state(x),
    controls = .control_state(x),
    n_batches = length(unique(x@observations$Batch)),
    n_runs = nrow(x@observations),
    at = Sys.time()
  )
}


#================== Derived-result Provenance ==================

#Store the fitted quantities that determine the surrogate's predictions. The
#campaign snapshot records the data and response standardisation used by the fit.
.surrogate_state <- function(fit) {
  list(
    hp = fit$hp[c(
      "ell",
      "sigma_f2",
      "sigma_b2",
      "sigma2",
      "kernel"
    )],
    beta0 = fit$beta0,
    jitter = fit$jitter,
    y_mean = fit$y_mean,
    y_sd = fit$y_sd,
    sign = fit$sign,
    naive = fit$naive
  )
}


#Attach the fitted surrogate state to the campaign snapshot used for fitting.
.provenance_surrogate <- function(x, fit) {
  provenance <- .provenance(x)
  provenance$kind <- "surrogate"
  provenance$surrogate <- .surrogate_state(fit)
  provenance
}


#Attach the generating surrogate's name and fitted state and the acquisition
#settings to the campaign snapshot from which a proposal was calculated. The
#proposed conditions themselves are stored with the proposal.
.provenance_proposal <- function(x, surrogate_name, method, settings) {
  provenance <- .provenance(x)
  provenance$kind <- "proposal"
  provenance$surrogate_name <- surrogate_name
  provenance$surrogate <- x@surrogates[[surrogate_name]]$provenance$surrogate
  provenance$method <- method
  provenance$settings <- settings
  provenance
}


#A space-filling proposal involves no surrogate: only the campaign snapshot and
#the design settings are recorded.
.provenance_space_filling_proposal <- function(x, settings) {
  provenance <- .provenance(x)
  provenance$kind <- "space_filling_proposal"
  provenance$settings <- settings
  provenance
}


#An analysis interprets one fitted surface, so it inherits the snapshot of the
#campaign the surrogate was fitted to (including its fitted state), plus the
#surrogate's name and the analysis settings.
.provenance_analysis <- function(x, surrogate_name, kind, settings) {
  provenance <- x@surrogates[[surrogate_name]]$provenance
  provenance$kind <- "analysis"
  provenance$at <- Sys.time()
  provenance$analysis <- kind
  provenance$surrogate_name <- surrogate_name
  provenance$settings <- settings
  provenance
}


#================== Current and Historical Results ==================

#Identify whether a stored result belongs to the current campaign state. Any
#other result remains a valid historical snapshot; nothing is deleted merely
#because data were appended, replaced, or removed.
.provenance_status <- function(provenance, x) {
  reasons <- character(0)

  if (!identical(provenance$state, .model_state(x))) {
    reasons <- c(reasons, "the campaign data or response specification differs")
  }

  if (
    provenance$kind %in% c("proposal", "space_filling_proposal") &&
    !identical(provenance$controls, .control_state(x))
  ) {
    reasons <- c(reasons, "the declared reference conditions differ")
  }

  #Proposals and analyses also depend on their surrogate's stored fit.
  if (provenance$kind %in% c("proposal", "analysis")) {
    surrogate_name <- provenance$surrogate_name

    if (!surrogate_name %in% names(x@surrogates)) {
      reasons <- c(
        reasons,
        paste0('surrogate "', surrogate_name, '" is no longer stored')
      )
    } else {
      current_surrogate <-
        x@surrogates[[surrogate_name]]$provenance$surrogate

      if (!identical(provenance$surrogate, current_surrogate)) {
        reasons <- c(
          reasons,
          paste0('surrogate "', surrogate_name, '" differs from the stored fit')
        )
      }
    }
  }

  current <- !length(reasons)
  list(
    current = current,
    status = if (current) "current" else "historical",
    reason = if (current) NULL else paste(unique(reasons), collapse = "; ")
  )
}


#List stored analyses of one type.
.analysis_names <- function(x, kind) {
  if (!methods::is(x, "BOCampaign")) {
    stop('"x" must be a BOCampaign.', call. = FALSE)
  }

  names(x@analyses)[vapply(
    x@analyses,
    function(entry) identical(entry$provenance$analysis, kind),
    logical(1)
  )]
}


#Retrieve and check a stored analysis of the requested type.
.analysis_entry <- function(
    x,
    name,
    kind,
    strict
) {
  strict <- .validate_flag(strict, "strict")
  available <- .analysis_names(x, kind)  #Also checks that x is a BOCampaign

  if (is.null(name)) {
    if (!length(available)) {
      stop('The campaign contains no stored ', kind, " analysis.", call. = FALSE)
    }
    current <- available[vapply(
      x@analyses[available],
      function(entry) .is_current(entry$provenance, x),
      logical(1)
    )]
    if (!length(current)) {
      stop(
        "The campaign contains no current ", kind,
        " analysis. Specify a historical analysis by name. Available: ",
        paste(available, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    name <- tail(current, 1L)
  } else {
    name <- .validate_name(name)
  }

  if (!name %in% available) {
    listed <- if (length(available)) paste(available, collapse = ", ") else "none"
    stop(
      'No ', kind, ' analysis named "', name,
      '". Available: ', listed, ".",
      call. = FALSE
    )
  }

  entry <- x@analyses[[name]]
  .check_current(
    entry$provenance,
    x,
    paste0(kind, ' analysis "', name, '"'),
    strict = strict
  )
  entry
}


#A stored result is current when all of its dependencies remain unchanged.
.is_current <- function(provenance, x) {
  .provenance_status(provenance, x)$current
}


#Explicitly named historical results are returned normally. Strict retrieval is
#available for scripts that require an exact match to the current campaign.
.check_current <- function(
    provenance,
    x,
    what,
    strict = FALSE
) {
  status <- .provenance_status(provenance, x)

  if (status$current) {
    return(invisible(TRUE))
  }

  if (strict) {
    stop(
      what,
      " is historical because ",
      status$reason,
      ".",
      call. = FALSE
    )
  }

  invisible(FALSE)
}
