################################################################################
# RunSheet.R -- physical run ordering and laboratory files
################################################################################


#================== Reference Placement ==================

#Resolve the reference-placement policy. Without run-order randomization,
#references are presented last unless the user explicitly requests a
#deterministic spread or bracketed arrangement.
.resolve_reference_order <- function(reference_order, randomize) {
  randomize <- .validate_flag(randomize, "randomize")

  if (is.null(reference_order)) {
    return(if (randomize) "spread" else "last")
  }

  if (
    !is.character(reference_order) ||
    length(reference_order) != 1L ||
    is.na(reference_order)
  ) {
    stop(
      '"reference_order" must be one of "spread", "random", "last", or "bracket".',
      call. = FALSE
    )
  }

  reference_order <- match.arg(
    reference_order,
    c("spread", "random", "last", "bracket")
  )

  if (!randomize && identical(reference_order, "random")) {
    stop(
      '"reference_order = \'random\'" requires "randomize = TRUE".',
      call. = FALSE
    )
  }

  reference_order
}


#Select one position from each of equally sized strata spanning a run. With a
#single reference, use the middle third so that "spread" does not accidentally
#place the only anchor at an endpoint.
.spread_reference_positions <- function(total, n_reference, randomize, seed) {
  if (!n_reference) {
    return(integer(0))
  }

  if (n_reference == 1L) {
    lower <- floor(total / 3) + 1L
    upper <- max(lower, ceiling(2 * total / 3))

    if (randomize) {
      return(.with_seed(seed, sample.int(upper - lower + 1L, 1L) + lower - 1L))
    }

    return(as.integer(floor((lower + upper) / 2)))
  }

  boundaries <- floor(seq(0, total, length.out = n_reference + 1L))

  positions <- vapply(seq_len(n_reference), function(i) {
    lower <- boundaries[i] + 1L
    upper <- boundaries[i + 1L]

    if (randomize) {
      as.integer(.with_seed(
        .seed_offset(seed, i),
        sample.int(upper - lower + 1L, 1L) + lower - 1L
      ))
    } else {
      as.integer(floor((lower + upper) / 2))
    }
  }, integer(1))

  sort(positions)
}


#Identify the positions occupied by references under the requested policy.
.reference_positions <- function(
    n_sample,
    n_reference,
    reference_order,
    randomize,
    seed
) {
  total <- n_sample + n_reference

  if (!n_reference) {
    return(integer(0))
  }

  switch(
    reference_order,
    last = seq.int(total - n_reference + 1L, total),
    random = sort(.with_seed(seed, sample.int(total, n_reference))),
    spread = .spread_reference_positions(
      total = total,
      n_reference = n_reference,
      randomize = randomize,
      seed = seed
    ),
    bracket = {
      if (n_reference == 1L) {
        as.integer(ceiling(total / 2))
      } else if (n_reference == 2L) {
        c(1L, total)
      } else {
        interior <- .spread_reference_positions(
          total = total - 2L,
          n_reference = n_reference - 2L,
          randomize = randomize,
          seed = seed
        )
        c(1L, interior + 1L, total)
      }
    }
  )
}


#================== Run-Sheet Assembly ==================

#Turn physical Sample and Ref conditions into the sheet executed in the lab.
#Sample order and reference placement use separate streams derived from the
#single seed supplied to the high-level proposal function.
.prepare_run_sheet <- function(
    conditions,
    response,
    randomize = FALSE,
    reference_order = NULL,
    seed = NULL
) {
  conditions <- as.data.frame(conditions, stringsAsFactors = FALSE)
  randomize <- .validate_flag(randomize, "randomize")
  seed <- .validate_seed(seed)
  reference_order <- .resolve_reference_order(reference_order, randomize)

  if (!"Type" %in% names(conditions)) {
    stop('Run-sheet conditions require a "Type" column.', call. = FALSE)
  }

  if (!"Batch" %in% names(conditions)) {
    stop('Run-sheet conditions require a "Batch" column.', call. = FALSE)
  }

  if (!all(conditions$Type %in% c("Sample", "Ref"))) {
    stop('Run-sheet "Type" must contain only "Sample" or "Ref".', call. = FALSE)
  }

  reserved <- intersect(c("Run", response), names(conditions))
  if (length(reserved)) {
    stop(
      "Run-sheet conditions already contain reserved column(s): ",
      paste(reserved, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  sample_index <- which(conditions$Type == "Sample")
  reference_index <- which(conditions$Type == "Ref")

  if (randomize && length(sample_index) > 1L) {
    sample_index <- .with_seed(
      .seed_offset(seed, 1L),
      sample(sample_index, length(sample_index))
    )
  }

  if (randomize && length(reference_index) > 1L) {
    reference_index <- .with_seed(
      .seed_offset(seed, 2L),
      sample(reference_index, length(reference_index))
    )
  }

  reference_positions <- .reference_positions(
    n_sample = length(sample_index),
    n_reference = length(reference_index),
    reference_order = reference_order,
    randomize = randomize,
    seed = .seed_offset(seed, 3L)
  )

  order_index <- integer(nrow(conditions))
  order_index[reference_positions] <- reference_index
  order_index[setdiff(seq_len(nrow(conditions)), reference_positions)] <-
    sample_index

  sheet <- conditions[order_index, , drop = FALSE]
  rownames(sheet) <- NULL
  sheet$Run <- seq_len(nrow(sheet))
  sheet[[response]] <- NA_real_

  leading <- c("Run", "Batch", "Type")
  sheet <- sheet[, c(leading, setdiff(names(sheet), leading)), drop = FALSE]

  attr(sheet, "reference_order") <- reference_order
  attr(sheet, "randomized") <- randomize
  sheet
}


#================== Laboratory Files ==================

#Write one run sheet to a CSV file, replacing an existing file of that name.
#Blank response cells are written as empty rather than the literal text "NA".
.write_run_sheet <- function(sheet, file) {
  if (
    !is.character(file) ||
    length(file) != 1L ||
    is.na(file) ||
    !nzchar(trimws(file))
  ) {
    stop('"file" must be one non-empty path.', call. = FALSE)
  }

  file <- trimws(file)
  if (!dir.exists(dirname(file))) {
    stop(
      "Run-sheet output directory does not exist: ",
      dirname(file),
      ".",
      call. = FALSE
    )
  }

  utils::write.csv(sheet, file = file, row.names = FALSE, na = "")
  invisible(file)
}
