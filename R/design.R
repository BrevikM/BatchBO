################################################################################
#Declaring factors, and coding to / from [0, 1]^D
################################################################################


#================== Defining factors ==================

#def_factor() defines one continuous experimental factor in physical units.
#design_space() combines one or more factors into the design used throughout a
#campaign. Internally, factors are encoded to [0, 1] before GP fitting and
#acquisition optimization. With transform = "log10", bounds and conditions stay
#in physical units; only the internal coding is logarithmic.

#' @export
def_factor <- function(label, lower, upper, transform = c('linear','log10')){
  transform <- match.arg(transform) #Identify transform
  label <- .validate_factor_names(label)
  
  bounds <- .validate_bounds(
    lower = lower,
    upper = upper,
    length_required = 1,
    name = "factor bounds",
    positive = identical(transform, "log10")
  )
  
  #Return structure
  structure(
    list(label = label,
         lower = bounds$lower,
         upper = bounds$upper,
         transform = transform),
    class = 'BO_Factor'
  ) 
}


#================== Defining experimental space ==================

#' @export
design_space <- function(...){
  vars <- list(...) #Collect the arbitrary number of experimental factors in a list
  
  
  #Define helper to extract particular aspects from each factor
  get_content <- function(fac, asp) {
    fac[[asp]]
  }
  
  #Ensure factors are entered
  if (length(vars) == 0){stop('The design space requires at least one experimental factor.', call. = FALSE)}
  
  #If none of the variables inherit the BO_Factor, they will return false, which will stop the process
  if (!all(vapply(vars, inherits, logical(1), 'BO_Factor'))){
    stop('Experimental factors must be defined by def_factor().', call. = FALSE)
  }
  
  #Ensure that factor labels are valid campaign column names.
  labels <- .validate_factor_names(
    vapply(vars, get_content, character(1), asp = "label")
  )
  
  #Define the S4 object
  new("DesignSpace",
      label      = labels,
      lower      = vapply(vars, get_content, numeric(1), asp = "lower"),
      upper      = vapply(vars, get_content, numeric(1), asp = "upper"),
      transform  = vapply(vars, get_content, character(1), asp = "transform")
  )
}


#================== Defining encoder ==================
.encode <- function(design, raw_df){
  #Identify potentially missing factors
  missing_factors <- setdiff(design@label, colnames(raw_df))
  
  if (length(missing_factors) > 0) {
    stop(
      "Missing factor column(s): ",
      paste(missing_factors, collapse = ", "),
      ".",
      call. = FALSE #Hide function identity from error call
    )
  }
  
  #Extract the design space matrix and the experimental boundaries
  X <- as.matrix(raw_df[, design@label, drop = FALSE]) #Keep drop false to ensure the matrix never converts to a vector
  
  low_encode <- design@lower
  high_encode <- design@upper
  trns_status <- which(design@transform == 'log10') #Select factors to be transformed
  
  #Transform data
  if (length(trns_status) > 0){
    low_encode[trns_status] <- log10(low_encode[trns_status])
    high_encode[trns_status] <- log10(high_encode[trns_status])
    
    
    # Transform observations and bounds to the log10 domain
    X[, trns_status] <- log10(
      X[, trns_status, drop = FALSE]
    )
  }
  
  
  
  #Encode the design to [0,1] by (X-Xmin)/(Xmax-Xmin)
  #Margin = 2 corresponds to performing an operation on a column, a = 1 on a row
  X <- sweep(X, MARGIN = 2, low_encode, '-') #Subtract
  X <- sweep(X, MARGIN = 2, high_encode - low_encode, '/') #Divide
  
  #Dinnames requires rownames (i.e. null) and colnames, i.e. design@label
  dimnames(X) <- list(NULL, design@label)
  X
}

#================== Defining decoder ==================
.decode <- function(design, X){
  X <- as.matrix(X)
  coded <- X
  low_encode <- design@lower
  high_encode <- design@upper
  trns_status <- which(design@transform == 'log10') #Select factors to be transformed 
  
  #If data have been transformed, their encoding levels must be established again
  if (length(trns_status) > 0){
    low_encode[trns_status] <- log10(low_encode[trns_status])
    high_encode[trns_status] <- log10(high_encode[trns_status])
  }
  
  X <- sweep(X, MARGIN = 2, high_encode - low_encode, '*')
  X <- sweep(X, MARGIN = 2, low_encode, '+')
  
  
  #Reverse transformed data
  if (length(trns_status) > 0){
    X[, trns_status] <- 10^(
      X[, trns_status, drop = FALSE]
    )
  }
  
  #Arithmetic in the physical transformation can place an encoded endpoint a
  #few floating-point units beyond its declared bound. Preserve genuinely
  #out-of-domain values, but return numerical boundary values exactly.
  tolerance <- 1e-8
  for (factor_index in seq_along(design@label)) {
    at_lower <- coded[, factor_index] >= -tolerance &
      coded[, factor_index] <= tolerance
    at_upper <- coded[, factor_index] >= 1 - tolerance &
      coded[, factor_index] <= 1 + tolerance
    
    X[at_lower, factor_index] <- design@lower[factor_index]
    X[at_upper, factor_index] <- design@upper[factor_index]
  }
  
  #Return the decoded data as a data.frame
  X <- as.data.frame(X)
  colnames(X) <- design@label
  X
}


#================== Space-Filling Designs ==================

#Generate a randomized Latin hypercube in the coded unit cube. Every marginal
#contains exactly one point in each of the n equally sized strata.
.random_lhs <- function(n, D, seed = NULL) {
  .with_seed(seed, {
    X <- matrix(NA, nrow = n, ncol = D)
    
    for (j in seq_len(D)) {
      X[, j] <- (sample.int(n) - stats::runif(n)) / n
    }
    
    X
  })
}


#Generate the Sample conditions for a response-free starting design. One global
#design is constructed before its rows are randomly distributed across the
#requested number of batches. Storage and physical run ordering are handled by
#proposeSpaceFilling(), which also validates the inputs.
.space_filling_samples <- function(
    x,
    n_batches = 1,
    q,
    method = c("lhs", "sobol"),
    seed = NULL
) {
  method <- match.arg(method)
  
  D <- length(x@design@label)
  n_samples <- n_batches * q
  
  coded <- switch(
    method,
    lhs = .random_lhs(
      n = n_samples,
      D = D,
      seed = seed
    ),
    sobol = .sbl(
      n = n_samples,
      D = D,
      lower = rep(0, D),
      upper = rep(1, D),
      seed = seed
    )
  )
  
  #Random allocation reduces the risk that a region of the global design is
  #systematically associated with an early or late laboratory batch.
  allocation <- .with_seed(
    .seed_offset(seed, 1),
    sample.int(n_samples)
  )
  coded <- coded[allocation, , drop = FALSE]
  
  batches <- lapply(seq_len(n_batches), function(batch_index) {
    first <- (batch_index - 1) * q + 1
    rows <- seq.int(first, first + q - 1)
    sample_batch <- cbind(
      data.frame(
        Type = rep("Sample", q),
        stringsAsFactors = FALSE
      ),
      .decode(x@design, coded[rows, , drop = FALSE])
    )
    
    rownames(sample_batch) <- NULL
    sample_batch
  })
  
  list(
    coded = coded,
    batches = unname(batches)
  )
}


#================== Inspecting experimental space ==================

#factorNames() returns the declared factor labels of a DesignSpace or BOCampaign.

#' @export
factorNames <- function(x) {
  if (methods::is(x, "BOCampaign")) {
    x <- x@design
  }
  
  if (!methods::is(x, "DesignSpace")) {
    stop('"x" must be a DesignSpace or BOCampaign.', call. = FALSE)
  }
  
  x@label
}


#as.data.frame() returns one row per factor with its physical bounds and coding.

#' @export
setMethod(
  "as.data.frame",
  "DesignSpace",
  function(x, row.names = NULL, optional = FALSE, ...) {
    data.frame(
      factor = x@label,
      lower = x@lower,
      upper = x@upper,
      transform = x@transform,
      row.names = row.names,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }
)