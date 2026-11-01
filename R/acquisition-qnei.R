#================== qNEI Step Preparation ==================

.qnei_step <- function(
    g,
    X_baseline,
    selected,
    innovations
) {
  #At greedy position j, compare the new candidate against
  #every distinct condition evaluated before the current batch and
  #proposals 1,...,j-1 already selected for the pending batch.
  incumbent_set <- rbind(X_baseline,selected)
  
  m <- nrow(incumbent_set)
  
  #Obtain p(f_C | y), the joint posterior at the incumbent locations.
  joint <- gp_joint(g,incumbent_set)
  U_incumbent <- joint$chol 
  
  #Use the first m innovation rows to generate draws at the incumbent
  #locations. Columns correspond to independent Monte Carlo simulations.
  incumbent_innovations <- innovations[seq_len(m), ,drop = FALSE]
  
  #The master innovations matrix contains one row for every baseline condition
  #and every possible proposal. The next unused row is reserved for the
  #candidate's conditionally independent residual component in each simulation.
  candidate_innovation <- innovations[m + 1, ]
  
  #For every Monte Carlo simulation s, draw:
  #f_C^(s) = mu_C + U^T e^(s) for every draw s (columns)
  incumbent_draws <- joint$mean + crossprod(U_incumbent, incumbent_innovations)
  
  #Find the best latent response currently available in every simulation.
  #The result contains one incumbent value for each Monte Carlo draw.
  current_best <- apply(incumbent_draws, MARGIN = 2, FUN = max)
  
  #Calculate the covariance between every observed training row
  #and every incumbent location:
  #
  #K_XC = sigma_f^2 R(X, incumbent_set).
  #
  #g$X retains duplicated measurements because V is defined over every
  #observation, whereas X_baseline contains only distinct latent locations.
  K_observed_incumbent <- g$hp$sigma_f2*.corr(g$X,incumbent_set,g$hp$ell,g$hp$kernel)
  
  #Apply the inverse training covariance:
  #
  #V^{-1} K_XC.
  #
  #The existing Cholesky factor g$U is reused, so V is neither reconstructed nor
  #factorized here. This term remains constant while the candidate location moves.
  V_inverse_K_observed_incumbent <- backsolve(
    g$U,forwardsolve(
      t(g$U),K_observed_incumbent)
  )
  
  #Return everything that remains fixed while Sobol screening and L-BFGS-B
  #evaluate possible locations for proposal j.
  list(
    incumbent_set = incumbent_set,
    U_incumbent = U_incumbent,
    incumbent_innovations = incumbent_innovations,
    candidate_innovation = candidate_innovation,
    current_best = current_best,
    V_inverse_K_observed_incumbent = V_inverse_K_observed_incumbent
  )
}



.qnei_candidate_draws <- function(
    g,
    X,
    step
) {
  X <- as.matrix(X)
  
  K_candidate_observed <- g$hp$sigma_f2*.corr(X, g$X, g$hp$ell, g$hp$kernel)
  
  candidate_mean <- g$beta0 + drop(K_candidate_observed %*% g$a)
  projected_candidate <- forwardsolve(t(g$U), t(K_candidate_observed))
  candidate_variance <- pmax(g$hp$sigma_f2 - colSums(projected_candidate^2), 0)
  candidate_incumbent_covariance <- g$hp$sigma_f2*.corr(X, step$incumbent_set,
      g$hp$ell,g$hp$kernel) - K_candidate_observed %*% step$V_inverse_K_observed_incumbent
  
  #Map the incumbent innovations into their correlated contribution to the
  #candidate draw.
  A <- t(forwardsolve(t(step$U_incumbent),
      t(candidate_incumbent_covariance))
  )
  
  candidate_conditional_sd <- sqrt(
    pmax(
      candidate_variance +
        g$jitter -
        rowSums(A^2),
      0
    )
  )
  
  #Return candidate xi (row i) and its value in each of the S draws/columns
  candidate_mean + A %*% step$incumbent_innovations +
    outer(candidate_conditional_sd, step$candidate_innovation)
}


.qnei_objective <- function(g,step){
  function(X) {
    candidate_draws <- .qnei_candidate_draws(
      g = g,
      X = X,
      step = step
    )
    #For every candidate and every draw, it computes how much the candidate beats;
    #the value to beat in that draw, and 0 if it doesn't
    improvement <- pmax(
      sweep(candidate_draws, MARGIN = 2, STATS = step$current_best, FUN = "-"),
      0)
    
    rowMeans(improvement) #Averages each row across the S columns, 
    #a set of these are then used for further optimization
  }
}


#================== qNEI Batch Acquisition ==================
acq_qnei <- function(
    g,
    q,
    S_n = 256,
    n_screen = 4096,
    restarts = 12,
    maxit = 500,
    lower = NULL,
    upper = NULL,
    seed = NULL,
    max_screen_factor = 4,
    X_baseline = NULL
) {
  S_n <- .validate_count(S_n, "S_n")
  D <- ncol(g$X)
  
  #Measured conditions in the set to beat; by default every distinct training location
  if (is.null(X_baseline)) X_baseline <- unique(g$X)
  
  n_baseline <- nrow(X_baseline)
  
  .with_seed(seed, {
    innovations <- matrix(rnorm((n_baseline + q) * S_n),
      nrow = n_baseline + q,
      ncol = S_n
    )
    
    .greedy_cont(
      g = g,
      q = q,
      state = matrix(numeric(0), nrow = 0, ncol = D),
      make_acquisition = function(selected, j) {
        step <- .qnei_step(
          g = g,
          X_baseline = X_baseline,
          selected = selected,
          innovations = innovations
        )
        
        .qnei_objective(
          g = g,
          step = step
        )
      },
      update_state = function(selected, x, j) {
        rbind(selected, x)
      },
      n_screen = n_screen,
      restarts = restarts,
      maxit = maxit,
      lower = lower,
      upper = upper,
      seed = seed,
      max_screen_factor = max_screen_factor
    )
  })
}