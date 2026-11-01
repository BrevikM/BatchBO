#================== Sobol Sequence ==================

.sbl <- function(n, D, lower, upper, seed = NULL){
  #Define the Sobol matrix for acquisition function optimization. 
  X_sbl <- .with_seed(
    seed,
    matrix(
      qrng::sobol(n, d = D, randomize = 'digital.shift'),
      nrow = n,
      ncol = D
    )
  )
  
  #Rescale from Sobol [0,1] to acquisition optimizer's settings
  X_sbl <- sweep(X_sbl, MARGIN = 2, upper - lower, '*')
  X_sbl <- sweep(X_sbl, MARGIN = 2, lower, '+')
  
  X_sbl
} 

#================== Acqusition Optimizer ==================
#A Sobol screen identifies promising starting points, which are expanded with
#larger sets if no positive values are found. Automatic expansion is disabled for TS

.screen_acquisition <- function(
    objective, screen, lower, upper,
    max_screen_factor = 4, seed = NULL
) {
  n_screen <- nrow(screen)
  
  #Evaluate the acquisition function on the initial Sobol screen
  values <- objective(screen)
  attempts <- 1
  
  #If no screening values are positive and rescreening is enabled, evaluate
  #progressively larger Sobol screens.
  if (max(values) <= 0 && max_screen_factor > 1) {
    for (factor in seq.int(2, max_screen_factor)) {
      #Update the seed so that each Sobol screen explores new conditions.
      replacement_seed <- .seed_offset(seed, factor - 1)
      
      screen <- .sbl(
        n = as.integer(n_screen*factor),
        D = length(lower),
        lower = lower,
        upper = upper,
        seed = replacement_seed
      )
      
      values <- .validate_screen_values(
        objective(screen),
        nrow(screen)
      )
      attempts <- attempts + 1
      
      #Stop enlarging the screen once a positive value has been found.
      if (max(values) > 0) {
        break
      }
    }
  }
  
  list(
    screen = screen,
    values = values,
    attempts = attempts
  )
}




#================== Acqusition Optimizer ==================
#Maximize the acquisition function over the bounded domain using the Sobol screen. 

.optimize_acquisition <- function(
    objective, screen, lower, upper,
    restarts, maxit,
    max_screen_factor = 4, seed = NULL
) {
  screening <- .screen_acquisition(
    objective = objective,
    screen = screen,
    lower = lower,
    upper = upper,
    max_screen_factor = max_screen_factor,
    seed = seed
  )
  
  screen <- screening$screen
  values <- screening$values
  
  
  #Select the best screened points as starting locations for further optimization
  order_best <- order(values, decreasing = TRUE)
  top_indices <- order_best[seq_len(restarts)]
  top <- screen[top_indices, , drop = FALSE]
  
  #Convert each candidate vector into the one-row matrix expected by the
  #acquisition function. Negation converts maximization into minimiztion expected by optim()
  fn <- function(x) {
    -.validate_objective_value(
      objective(matrix(data = x, nrow = 1))
    )
  }
  
  fits <- lapply(seq_len(restarts), function(i) {
    stats::optim(
      par = top[i, ],
      fn = fn,
      method = "L-BFGS-B",
      lower = lower,
      upper = upper,
      control = list(maxit = maxit)
    )
  })
  
  #The smallest minimized value corresponds to the largest acquisition value.
  best_fit <- fits[[
    which.min(vapply(fits, function(fit) fit$value, numeric(1)))
  ]]
  
  list(
    x = best_fit$par,
    value = -best_fit$value,
    convergence = best_fit$convergence,
    message = best_fit$message,
    screen_attempts = screening$attempts,
    screen_size = nrow(screen),
    screen_max = max(values)
  )
}


#================== Construct Batches Greedily ==================
.greedy_cont <- function(
    g,
    q,
    state,
    make_acquisition,
    update_state,
    n_screen = 4096,
    restarts = 12,
    maxit = 500,
    lower = NULL,
    upper = NULL,
    seed = NULL,
    max_screen_factor = 4
) {
  D <- ncol(g$X) #Define the number of factors

  controls <- .validate_acquisition_controls(
    n_screen = n_screen,
    restarts = restarts,
    maxit = maxit,
    lower = lower,
    upper = upper,
    D = D,
    max_screen_factor = max_screen_factor,
    seed = seed
  )
  n_screen <- controls$n_screen
  restarts <- controls$restarts
  maxit <- controls$maxit
  lower <- controls$lower
  upper <- controls$upper
  max_screen_factor <- controls$max_screen_factor
  seed <- controls$seed
  
  #Create a common Sobol set as the basis for acquisition maximization.
  screen <- .sbl(
    n = n_screen,
    D = D,
    lower = lower,
    upper = upper,
    seed = seed
  )
  
  #Create a q times D matrix of NAs for proposal storage.
  X <- matrix(NA, q, D)
  
  #Assign placeholder values for the proposal diagnostics defined below.
  value <- numeric(q)
  screen_max <- numeric(q)
  
  convergence <- integer(q)
  screen_attempts <- integer(q)
  screen_size <- integer(q)
  
  optimizer_message <- character(q)
  
  #Construct the batch one proposal at a time, i.e. greedily.
  for (j in seq_len(q)) {
    objective <- make_acquisition(state, j)
    
    result <- .optimize_acquisition(
      objective = objective,
      screen = screen,
      lower = lower,
      upper = upper,
      restarts = restarts,
      maxit = maxit,
      max_screen_factor = max_screen_factor,
      seed = seed
    )
    
    X[j, ] <- result$x
    value[j] <- result$value
    convergence[j] <- result$convergence
    screen_max[j] <- result$screen_max
    screen_attempts[j] <- result$screen_attempts
    screen_size[j] <- result$screen_size
    
    optimizer_message[j] <- .optimizer_message(
      result$convergence,
      result$message
    )
    
    #Update the acquisition state before selecting the next batch member.
    state <- update_state(state, result$x, j)
  }
  
  #Store the optimized acquisition value for each proposed condition.
  attr(X, "value") <- value
  
  #Store the L-BFGS-B convergence code for each proposed condition.
  attr(X, "convergence") <- convergence
  
  #Store the L-BFGS-B termination message for each proposed condition.
  attr(X, "message") <- optimizer_message
  
  #Store the largest acquisition value found during the final Sobol screen.
  attr(X, "screen_max") <- screen_max
  
  #Store the number of Sobol-screening attempts used for each proposal.
  attr(X, "screen_attempts") <- screen_attempts
  
  #Store the number of candidate points in the final Sobol screen.
  attr(X, "screen_size") <- screen_size
  
  #Store the smallest length scale-normalized Euclidean distance between batch points.
  attr(X, "separation") <- if (q > 1) {
    min(
      stats::dist(
        sweep(X, 2, g$hp$ell, "/")
      )
    )
  } else {
    NA
  }
  
  X
}


#================== Acquisition Dispatch ==================

.acquisition_function <- function(
    method = c("kb", "ts", "qnei")
) {
  method <- match.arg(method)
  
  switch(
    method,
    kb = acq_kb,
    ts = acq_ts,
    qnei = acq_qnei
  )
}
