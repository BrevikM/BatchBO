################################################################################
# AllClasses.R -- class definitions and validity
################################################################################

#DesignSpace stores the declared factor names, physical bounds, and coding
#transformations. Construct it with design_space() rather than by modifying
#slots directly.
#' @exportClass DesignSpace
setClass(
  "DesignSpace",
  slots = c(
    label     = "character",
    lower     = "numeric",
    upper     = "numeric",
    transform = "character"  # "linear" | "log10"
  ),
  prototype = list(
    label = character(0), lower = numeric(0), upper = numeric(0),
    transform = character(0)
  )
)

setValidity("DesignSpace", function(object) {
  e <- character(0)
  D <- length(object@label)
  if (D < 1L) return("A design space needs at least one variable.")

  lens <- c(lower = length(object@lower), upper = length(object@upper),
            transform = length(object@transform))
  bad <- names(lens)[lens != D]
  if (length(bad)) {
    e <- c(e, sprintf("Slot(s) %s must have length %d to match `label`.",
                      paste(bad, collapse = ", "), D))
    return(e)  # nothing below is safe if lengths disagree
  }

  if (anyDuplicated(tolower(object@label))) {
    e <- c(e, "Factor labels must be unique, ignoring capitalization.")
  }
  if (any(!nzchar(trimws(object@label)))) {
    e <- c(e, "Factor labels must be non-empty.")
  }
  #Same reserved list as design_space() (validation.R)
  if (any(tolower(object@label) %in% tolower(.reserved_factor_names()))) {
    e <- c(e, "Factor labels must not use package-reserved output names.")
  }
  if (!all(object@transform %in% c("linear", "log10"))) {
    e <- c(e, "`transform` must be \"linear\" or \"log10\".")
  }
  if (any(!is.finite(object@lower)) || any(!is.finite(object@upper))) {
    e <- c(e, "`lower` and `upper` must be finite.")
  } else if (any(object@lower >= object@upper)) {
    at <- object@label[object@lower >= object@upper]
    e <- c(e, sprintf("`lower` must be strictly below `upper` for: %s.",
                      paste(at, collapse = ", ")))
  }
  lg <- object@transform == "log10"
  if (any(lg) && any(object@lower[lg] <= 0)) {
    at <- object@label[lg & object@lower <= 0]
    e <- c(e, sprintf("log10 variables need a strictly positive `lower`: %s.",
                      paste(at, collapse = ", ")))
  }

  if (length(e)) e else TRUE
})


setMethod("show", "DesignSpace", function(object) {
  cat("Design space with", length(object@label), "factor(s):\n")
  print(as.data.frame(object), row.names = FALSE)
  invisible(object)
})


#BOCampaign is the resumable campaign object: the design, response specification
#and reference condition (fixed after construction), the completed observations
#and their coded form, the response standardisation, and the named surrogates,
#analyses and proposals. Completed data are added through addBatch().
#' @exportClass BOCampaign
setClass(
  "BOCampaign",
  slots = c(
    design             = "DesignSpace",
    response           = "character",  # column name in the lab data
    response_transform = "character",  # "identity" | "log"
    direction          = "character",  # "maximize" | "minimize"

    controls           = "data.frame", # one fixed reference, physical units
    control_replicates = "integer",    # injections of the fixed reference

    observations       = "data.frame", # one row per run; grows via addBatch()
    coded              = "matrix",     # n x D in [0,1]; derived from the above
    standardisation    = "list",       # y_mean, y_sd, sign -- recomputed each add

    surrogates         = "list",
    analyses           = "list",
    proposals          = "list"
  ),
  prototype = list(
    response = NA_character_, response_transform = "identity",
    direction = "maximize",
    controls = data.frame(), control_replicates = integer(0),
    observations = data.frame(), coded = matrix(numeric(0), 0L, 0L),
    standardisation = list(y_mean = NA_real_, y_sd = NA_real_, sign = 1),
    surrogates = list(), analyses = list(), proposals = list()
  )
)

setValidity("BOCampaign", function(object) {
  e <- character(0)
  D <- length(object@design@label)
  n <- nrow(object@observations)

  # --- response specification ------------------------------------------------
  if (
    length(object@response) != 1L ||
    is.na(object@response) ||
    !nzchar(trimws(object@response))
  ) {
    e <- c(e, "`response` must be a single non-empty column name.")
  } else if (
    tolower(object@response) %in%
      tolower(c(
        object@design@label,
        .structural_column_names()
      ))
  ) {
    e <- c(e, "`response` must not duplicate a factor or reserved column name.")
  }
  if (!identical(length(object@response_transform), 1L) ||
      !object@response_transform %in% c("identity", "log")) {
    e <- c(e, "`response_transform` must be \"identity\" or \"log\".")
  }
  if (!identical(length(object@direction), 1L) ||
      !object@direction %in% c("maximize", "minimize")) {
    e <- c(e, "`direction` must be \"maximize\" or \"minimize\".")
  }

  # --- the alignment contract ------------------------------------------------
  if (nrow(object@coded) != n) {
    e <- c(e, sprintf("`coded` has %d rows but `observations` has %d.",
                      nrow(object@coded), n))
  }
  if (ncol(object@coded) != D) {
    e <- c(e, sprintf("`coded` has %d columns but the design declares %d inputs.",
                      ncol(object@coded), D))
  } else if (n > 0L && !identical(colnames(object@coded), object@design@label)) {
    e <- c(e, "`coded` column names must match the factor labels, in order.")
  }
  # Everything below indexes `coded` by vectors derived from `observations`,
  # so a shape mismatch must stop here or it surfaces as an opaque subscript
  # error instead of the real problem.
  if (length(e)) return(e)

  # --- observations schema ---------------------------------------------------
  if (n > 0L) {
    need <- c("Batch", "Type", object@design@label, object@response)
    miss <- setdiff(need, names(object@observations))
    if (length(miss)) {
      e <- c(e, sprintf("`observations` is missing column(s): %s.",
                        paste(miss, collapse = ", ")))
    } else {
      if (!all(object@observations$Type %in% c("Sample", "Ref"))) {
        e <- c(e, "`Type` must contain only \"Sample\" or \"Ref\".")
      }
      if (anyNA(object@observations[, need])) {
        e <- c(e, "`observations` contains NA in a required column.")
      }
      if (any(!nzchar(trimws(as.character(object@observations$Batch))))) {
        e <- c(e, "Batch labels in `observations` must be non-empty.")
      }
      numeric_columns <- c(object@design@label, object@response)
      if (!all(vapply(
        object@observations[, numeric_columns, drop = FALSE],
        is.numeric,
        logical(1)
      ))) {
        e <- c(e, "Factor and response columns in `observations` must be numeric.")
      } else if (any(!is.finite(as.matrix(
        object@observations[, numeric_columns, drop = FALSE]
      )))) {
        e <- c(e, "Factor and response columns in `observations` must be finite.")
      }
      if (
        identical(object@response_transform, "log") &&
        is.numeric(object@observations[[object@response]]) &&
        any(object@observations[[object@response]] <= 0)
      ) {
        e <- c(e, "A log-transformed response must be strictly positive.")
      }
      log_factors <- object@design@label[
        object@design@transform == "log10"
      ]
      if (
        length(log_factors) &&
        all(vapply(
          object@observations[, log_factors, drop = FALSE],
          is.numeric,
          logical(1)
        )) &&
        any(as.matrix(object@observations[, log_factors, drop = FALSE]) <= 0)
      ) {
        e <- c(e, "log10-transformed factor values in `observations` must be positive.")
      }
      if (any(!is.finite(object@coded))) {
        e <- c(e, "`coded` contains non-finite values.")
      } else {
        # Sample rows must lie in the box the acquisition searches. Reference
        # rows need not: a control is frequently an established standard method
        # that predates the design space, and it is an anchor for the batch
        # effect rather than a candidate. The optimiser never proposes outside
        # [0, 1], so an out-of-box control cannot mislead it.
        smp <- object@observations$Type == "Sample"
        if (any(smp)) {
          bad <- object@coded[smp, , drop = FALSE] < -1e-8 |
                 object@coded[smp, , drop = FALSE] > 1 + 1e-8
          if (any(bad)) {
            out <- object@design@label[apply(bad, 2L, any)]
            e <- c(e, sprintf("Sample rows have coded values outside [0, 1] for: %s.",
                              paste(out, collapse = ", ")))
          }
        }
      }
    }
  }

  # --- controls --------------------------------------------------------------
  if (nrow(object@controls)) {
    if (nrow(object@controls) != 1L) {
      e <- c(e, "`controls` must contain exactly one fixed reference condition.")
    }
    miss <- setdiff(object@design@label, names(object@controls))
    if (length(miss)) {
      e <- c(e, sprintf("`controls` is missing column(s): %s.",
                        paste(miss, collapse = ", ")))
    }
    if (length(object@control_replicates) != 1L) {
      e <- c(e, "`control_replicates` must contain one replicate count.")
    } else if (any(object@control_replicates < 1L)) {
      e <- c(e, "`control_replicates` must be at least 1.")
    }
    if (!length(miss)) {
      if (!all(vapply(
        object@controls[, object@design@label, drop = FALSE],
        is.numeric,
        logical(1)
      ))) {
        e <- c(e, "Factor columns in `controls` must be numeric.")
      } else if (any(!is.finite(as.matrix(
        object@controls[, object@design@label, drop = FALSE]
      )))) {
        e <- c(e, "Factor columns in `controls` must be finite.")
      }
      log_factors <- object@design@label[
        object@design@transform == "log10"
      ]
      if (
        length(log_factors) &&
        all(vapply(
          object@controls[, log_factors, drop = FALSE],
          is.numeric,
          logical(1)
        )) &&
        any(as.matrix(object@controls[, log_factors, drop = FALSE]) <= 0)
      ) {
        e <- c(e, "log10-transformed factor values in `controls` must be positive.")
      }
    }
  }

  # --- derived slots ---------------------------------------------------------
  for (sl in c("surrogates", "analyses", "proposals")) {
    v <- slot(object, sl)
    if (length(v)) {
      nm <- names(v)
      if (is.null(nm) || any(!nzchar(nm))) {
        e <- c(e, sprintf("Every entry in `%s` must be named.", sl))
      } else if (anyDuplicated(nm)) {
        e <- c(e, sprintf("`%s` names must be unique.", sl))
      }
      if (!all(vapply(v, function(z) is.list(z) && !is.null(z$provenance),
                      logical(1)))) {
        e <- c(e, sprintf("Every entry in `%s` must carry a `provenance` record.", sl))
      }
    }
  }

  if (length(e)) e else TRUE
})


setMethod("show", "BOCampaign", function(object) {
  n_batches <- if (nrow(object@observations)) {
    length(unique(object@observations$Batch))
  } else {
    0L
  }

  cat("BatchBO campaign\n")
  cat("  Factors:            ", length(object@design@label), "\n", sep = "")
  cat("  Factor names:       ", paste(object@design@label, collapse = ", "), "\n", sep = "")
  cat(
    "  Response:           ", object@response,
    " (", object@direction, "; ", object@response_transform, ")\n",
    sep = ""
  )
  if (nrow(object@controls)) {
    cat(
      "  References:         ", nrow(object@controls),
      " condition, ", object@control_replicates,
      " injection(s) per batch\n",
      sep = ""
    )
  } else {
    cat("  References:         none\n")
  }
  cat("  Completed runs:     ", nrow(object@observations), "\n", sep = "")
  cat("  Completed batches:  ", n_batches, "\n", sep = "")
  cat("  Stored surrogates:  ", length(object@surrogates), "\n", sep = "")
  cat("  Stored analyses:    ", length(object@analyses), "\n", sep = "")
  cat("  Stored proposals:   ", length(object@proposals), "\n", sep = "")

  if (length(object@surrogates)) {
    cat("  Surrogate names:    ", paste(names(object@surrogates), collapse = ", "), "\n", sep = "")
  }
  if (length(object@proposals)) {
    cat("  Proposal names:     ", paste(names(object@proposals), collapse = ", "), "\n", sep = "")
  }
  if (length(object@analyses)) {
    cat("  Analysis names:     ", paste(names(object@analyses), collapse = ", "), "\n", sep = "")
  }

  invisible(object)
})
