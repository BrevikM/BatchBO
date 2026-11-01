################################################################################
# Randomness.R -- reproducible seeded work without changing the caller's RNG
################################################################################


#Evaluate code under an optional local seed. An explicit seed is reproducible
#and leaves the caller's random-number state exactly as it was on entry.
.with_seed <- function(seed, code) {
  seed <- .validate_seed(seed)

  if (is.null(seed)) {
    return(force(code))
  }

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    previous_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }

  on.exit(
    if (had_seed) {
      assign(".Random.seed", previous_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    },
    add = TRUE
  )

  set.seed(seed)
  force(code)
}


#Offset a supplied seed for a later Sobol screen without integer overflow.
.seed_offset <- function(seed, offset) {
  if (is.null(seed)) {
    return(NULL)
  }

  as.integer((as.double(seed) + as.double(offset)) %% .Machine$integer.max)
}
