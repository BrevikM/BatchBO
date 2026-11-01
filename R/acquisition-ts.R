#================== Spectral Kernel Representations ==================
.rff_spectral <- function(F_n, D, ell, kernel) {
  
  #Generate a matrix of standard normal draws, one row per Fourier feature
  frequencies <- matrix(
    rnorm(F_n * D),
    nrow = F_n,
    ncol = D)
  
  #For the RBF kernel, the multivariate normal corresponds to the spectral distribution
  #For the Materns, they must be scaled to multivariate t-distributions
  if (kernel != "rbf") {
    nu <- switch(
      kernel,
      matern32 = 3,
      matern52 = 5)
    
    #One chi-squared draw is shared across all dimensions of each feature.
    #This produces the multivariate-t spectrum of the radial Matern kernel.
    eta <- stats::rchisq(F_n, df = nu)
    
    frequencies <- sweep(
      frequencies,
      MARGIN = 1,
      STATS = sqrt(nu/eta),
      FUN = "*")
  }
  
  #Apply the fitted ARD length scales.
  sweep(
    frequencies,
    MARGIN = 2,
    STATS = ell,
    FUN = "/")
}



#================== Gaussian-process posterior path ==================

gp_path <- function(g, F_n) {
  X <- as.matrix(g$X)
  n <- nrow(X)
  D <- ncol(X)
  hp <- g$hp
  
  frequencies <- .rff_spectral(
    F_n = F_n,
    D = D,
    ell = hp$ell,
    kernel = hp$kernel
  )
  
  phase <- runif(F_n, 0, 2*pi)
  weight <- rnorm(F_n)
  
  #Draw an approximate prior path from the latent kernel.
  prior_path <- function(Xnew) {
    Xnew <- as.matrix(Xnew)
    
    #F*m matrix of features cos(omega_m^T x_i + psi_m), one row per feature
    features <- cos(frequencies %*% t(Xnew) + phase)
    
    #Evaluate the prior draw at every point in Xnew
    g$beta0 + sqrt(2*hp$sigma_f2/F_n)*drop(crossprod(features, weight))
  }
  
  #Simulate the measurements the prior draw would have produced: one offset per
  #batch, shared by its runs, plus an independent residual per run (jitter
  #included so the simulated covariance matches V)
  batch_levels <- unique(as.character(g$batch))
  batch_index <- match(as.character(g$batch), batch_levels)
  
  batch_draw <- rnorm(length(batch_levels), mean = 0, sd = sqrt(hp$sigma_b2))
  residual_draw <- rnorm(n, mean = 0, sd = sqrt(hp$sigma2 + g$jitter))
  prior_observation <- prior_path(X) + batch_draw[batch_index] + residual_draw
  
  #Find the Matheron correction weights, V^{-1}(y - y_tilde)
  correction <- backsolve(g$U,
                          forwardsolve(t(g$U), g$y - prior_observation))
  
  #Return the batch-free latent posterior path.
  function(Xnew) {
    Xnew <- as.matrix(Xnew)
    cross_covariance <- hp$sigma_f2*.corr(Xnew, X, hp$ell, hp$kernel)
    prior_path(Xnew) + drop(cross_covariance%*%correction)
  }
}





#================== Thompson Sampling ==================
acq_ts <- function(
    g,
    q,
    F_n = 1024,
    n_screen = 4096,
    restarts = 12,
    maxit = 500,
    lower = NULL,
    upper = NULL,
    seed = NULL
) {
  F_n <- .validate_count(F_n, "F_n")

  .with_seed(seed, {
    #Draw one independent posterior latent path per batch member.
    paths <- lapply(
      seq_len(q),
      function(j) gp_path(g, F_n = F_n)
    )

    .greedy_cont(
      g = g,
      q = q,
      state = paths,

      #Proposal j maximizes posterior path j.
      make_acquisition = function(state, j) {
        state[[j]]
      },

      update_state = function(state, x, j) {
        state
      },

      n_screen = n_screen,
      restarts = restarts,
      maxit = maxit,
      lower = lower,
      upper = upper,
      seed = seed,

      #A Thompson path can legitimately be negative everywhere;
      #so rescreening based on positivity is disabled.
      max_screen_factor = 1
    )
  })
}
