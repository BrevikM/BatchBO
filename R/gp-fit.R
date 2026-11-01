################################################################################
#Fit a batch-aware (mixed-effects) Gaussian Process
################################################################################


#================== Define the list of available Kernels ==================
.corr_scalar <- function(d, kernel = c('matern52', 'matern32', 'rbf')) {
  kernel <- match.arg(kernel)
  
  switch(
    kernel,
    
    matern52 = {
      t <- sqrt(5)*d
      (1+t+t^2/3)*exp(-t)
    },
    
    matern32 = {
      t <- sqrt(3)*d
      (1+t)*exp(-t)
    },
    
    rbf = {
      exp(-0.5*d^2)
    }
  )
}


#================== Find the correlation between data points ==================
.corr <- function(X1, X2, ell, kernel) {
  #The correlation depends on the Euclidean distance after each
  #coordinate has been standardized by its corresponding length scale.
  #||a-b||^2 = ||a||^2 + ||b||^2 - 2a^Tb

  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  
  A <- sweep(X1, 2, ell, '/') #Standardize each column with the length-scale vector
  B <- sweep(X2, 2, ell, '/')
  
  normA2 <- rowSums(A^2) #Square each element and sum each row to a single scalar
  normB2 <- rowSums(B^2) #This corresponds to the ||a||^2 and ||b||^2 
  
  #Compute ||a||^2 + ||b||^2 for every pair of rows (i,j) from A and B.
  pair_norm <- outer(normA2, normB2, '+') 
  
  #Find a_i^Tb_j
  pair_dot <- tcrossprod(A, B) #Equivalent to A %*% t(B)
  
  d2 <- pair_norm-2*pair_dot 
  
  d <- sqrt(pmax(d2, 0)) #Truncate numerical negative values to zero, then obtain distances
  
  .corr_scalar(d, kernel = kernel)
}


#================== Batch membership ==================
.S <- function(batch) {
  batch <- as.character(batch)
  outer(batch, batch, "==")*1 #Compare every pair of labels. The *1 turns TRUE/FALSE into 1/0
}

#================== Covariance Matrix ==================
.V <- function(X, batch, hp, jitter){
  V <- hp$sigma_f2*.corr(X,X, hp$ell,hp$kernel)+ #Kernel covariance
    hp$sigma_b2*.S(batch)+ #Batch effects
    diag(hp$sigma2+jitter, nrow(X)) #Experimental noise
  V
}



#================== Gaussian Process ==================
#The implementation largely follows that detailed in Rasmussen and Williams (2006)
gp <- function(X, y, batch, hp, jitter = 1e-8, beta0 = NULL){
  n <- nrow(X)
  
  stopifnot(
    length(y) == n,
    length(batch) == n
  )
  
  #Decompose the covariance matrix using the Cholesky decomposition, V = U^TU
  U <- chol(.V(X, batch, hp, jitter))
  
  #To fit the GP, we need two matrix-products involving inverses, V^{-1}*1 and V^{-1}*y.
  ones <- rep(1, n) 
  #Since V = U^T U, we find V^{-1}*1 = (U^T U)^{-1}1 = U^{-1}(U^T)^{-1}1.
  #Forwardsolve applies (U^T)^{-1}1, and backsolve subsequently applies U^{-1}.
  Vi1 <- backsolve(U,
                    forwardsolve(t(U), ones)
)
  
  
  #Likewise, V^{-1}y = (U^T U)^{-1}y = U^{-1}(U^T)^{-1}y.
  #forwardsolve computes (U^T)^{-1}y, and backsolve subsequently applies U^{-1}.
  Viy <- backsolve(U, 
                   forwardsolve(t(U), y))
  
  #Unless beta0 has already been specified, estimate it 
  if (is.null(beta0)){
    beta0 <- sum(Vi1*y)/sum(Vi1)
  }
  
 
  #Define the terms for the marginal likelihood
  a <- Viy - beta0 * Vi1  #The weighted residual, a = V^{-1}y - B0V^{-1}1 = V^{-1}(y - B01).
  quad <- sum((y-beta0)*a) #(y - B01)^T * V^{-1} * (y - B01)
  logdet <- 2*sum(log(diag(U))) #Since V = U^T U,log|V| = log|U^T U| = 2 log|U| = 2 Sum_i log(U_ii).
  
  #Evaluate the Gaussian negative log marginal likelihood.
  nll <- 0.5*(quad+logdet+n*log(2*pi))
  
  structure(
    list(
      X = X,
      y = y,
      batch = batch,
      hp = hp,
      jitter = jitter,
      U = U,
      beta0 = beta0,
      a = a,
      quad = quad,
      logdet = logdet,
      nll = nll
    ), 
    class = 'GP'
  )
}


#================== Practitioner-Facing GP Summary ==================

.gp_kernel_label <- function(kernel) {
  switch(
    kernel,
    rbf = "RBF",
    matern32 = "Matern 3/2",
    matern52 = "Matern 5/2",
    as.character(kernel)
  )
}


.gp_convergence_label <- function(code) {
  if (!length(code) || is.na(code)) return("Not optimized")
  switch(
    as.character(code),
    "0"  = "Converged",
    "1"  = "Iteration limit reached (code 1)",
    "51" = "WARNING: code 51",
    "52" = "ERROR: code 52",
    paste0("code ", code)
  )
}


#Print and summarize a fitted Gaussian-process surrogate returned by
#surrogate(). Printing the fit gives a compact model overview.
#summary() adds variance components, ARD length scales, and ML-II diagnostics. 

#' @method print GP
#' @export
print.GP <- function(x, ...) {
  model <- summary(x)$model

  cat(model$model, " Gaussian Process\n", sep = "")
  cat("  Kernel:       ", model$kernel, "\n", sep = "")
  cat("  Observations: ", model$n_runs, "\n", sep = "")
  cat("  Batches:      ", model$n_batches, "\n", sep = "")
  cat("  Factors:      ", model$n_factors, "\n", sep = "")
  cat("  ML-II:        ", model$convergence, "\n", sep = "")

  invisible(x)
}


#' @method summary GP
#' @export
summary.GP <- function(object, scaled = TRUE, ...) {
  #Define hyperparameters and factor names
  hp <- object$hp 
  D <- ncol(object$X)
  factor_names <- colnames(object$X)
  starts <- hp$fit_diagnostics

  #Report length scales and potential boundary hits
  length_scales <- data.frame(
    factor = factor_names,
    length_scale = hp$ell,
    at_bound = unname(hp$at_bound[paste0("ell_", seq_len(D))])
  )

  #Return variance components. By default in standardized units.
  variance <- c(hp$sigma_f2, if (object$naive) NA else hp$sigma_b2, hp$sigma2)
  if (!scaled) variance <- variance*object$y_sd^2

  variance_components <- data.frame(
    component = c("Response Surface", "Between-Batch", "Residual"),
    variance = variance,
    sd = sqrt(variance),
    at_bound = unname(hp$at_bound[c("sigma_f2", "sigma_b2", "sigma2")])
  )
  
  #Final model
  model <- list(
    model = if (object$naive) "Batch-Naive" else "Batch-Aware",
    kernel = .gp_kernel_label(hp$kernel),
    response = object$response,
    direction = object$direction,
    log_response = object$response_transform == "log",
    n_runs = nrow(object$X),
    n_batches = length(unique(object$batch)),
    n_factors = D,
    nll = object$nll,
    convergence = .gp_convergence_label(hp$convergence)
  )

  #Optimizer diagnostics
  optimization <- list(
    n_starts = nrow(starts),
    n_converged = sum(starts$convergence == 0, na.rm = TRUE),
    message = hp$message
  )

  structure(
    list(
      model = model,
      variance_components = variance_components,
      length_scales = length_scales,
      optimization = optimization,
      starts = starts,
      scaled = scaled
    ),
    class = "summary.GP"
  )
}


#' @method print summary.GP
#' @export
print.summary.GP <- function(x, digits = 4, ...) {
  model <- x$model
  optimization <- x$optimization

  cat(model$model, " Gaussian Process\n\n", sep = "")
  cat("Data\n")
  cat("  Observations: ", model$n_runs, "\n", sep = "")
  cat("  Batches:      ", model$n_batches, "\n", sep = "")
  cat("  Factors:      ", model$n_factors, "\n\n", sep = "")

  cat("Model\n")
  cat("  Kernel:       ", model$kernel, "\n", sep = "")
  cat("  Response:     ", model$response, " (", model$direction, ")\n", sep = "")
  cat("  Log response: ", model$log_response, "\n", sep = "")
  cat("  NLL:          ", format(model$nll, digits = digits), "\n\n", sep = "")

  units <- if (x$scaled) {
    "standardized"
  } else if (model$log_response) {
    paste0("log ", model$response)
  } else {
    model$response
  }
  cat("Variance components (", units, " scale)\n", sep = "")
  print(x$variance_components, row.names = FALSE, digits = digits)

  cat("\nLength scales\n")
  print(x$length_scales, row.names = FALSE, digits = digits)

  cat("\nML-II optimization\n")
  cat("  Starts:       ", optimization$n_starts, "\n", sep = "")
  cat("  Converged:    ", optimization$n_converged, "\n", sep = "")
  cat("  Selected:     ", model$convergence, "\n", sep = "")

  invisible(x)
}


#================== Gaussian Process Prediction ==================
gp_predict <- function(g, Xstar){
  Xstar <- as.matrix(Xstar) #Point(s) to be predicted 
  K_star <- g$hp$sigma_f2 * .corr(g$X, Xstar, g$hp$ell, g$hp$kernel) #Training-prediction cross-covariance
  
  #Since V = U^T U, solve U^T W = K_star, giving W = (U^T)^{-1}K_star.
  #Then W^T W = K_star^T V^{-1} K_star
  W <- forwardsolve(t(g$U), K_star)
  
  list(
    mean = g$beta0 + drop(crossprod(K_star, g$a)), #mu* = B_0 + K_star^T V^{-1}(y - B_0 1), since a = V^{-1}(y - B_0 1), this reduces to B_0 + K_star^T a.
    var = pmax(g$hp$sigma_f2 - colSums(W * W), 0) #colSums evaluates the diagonals of the variance term, since only point predictions are considered
  )
}

#================== Posterior cross-covariance ==================
#Posterior covariance of the batch-free surface between the points in X1 and X2:
#C(X1, X2) = K(X1, X2) - K(X1, X)V^{-1}K(X, X2).
.posterior_cross_covariance <- function(g, X1, X2) {
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  
  K12 <- g$hp$sigma_f2 * .corr(
    X1, X2, g$hp$ell, g$hp$kernel
  )
  
  K1n <- g$hp$sigma_f2 * .corr(
    X1, g$X, g$hp$ell, g$hp$kernel
  )
  
  Kn2 <- g$hp$sigma_f2 * .corr(
    g$X, X2, g$hp$ell, g$hp$kernel
  )
  
  #Since V = U^T U, W1^T W2 = K(X1, X) V^{-1} K(X, X2)
  W1 <- forwardsolve(t(g$U), t(K1n))
  W2 <- forwardsolve(t(g$U), Kn2)
  
  K12 - crossprod(W1, W2)
}


#================== Joint posterior at a set of points ==================
#Posterior of the batch-free surface at the points in Xstar, given the data in g:
#f(Xstar) | y ~ N(mu, C). Returns mu, C and the Cholesky factor of C.
gp_joint <- function(g, Xstar) {
  Xstar <- as.matrix(Xstar)
  
  #Joint posterior covariance among the points: their cross-covariance with themselves
  C_star <- .posterior_cross_covariance(g, Xstar, Xstar)
  
  #Restore exact symmetry lost through floating-point arithmetic
  #and add a small diagonal jitter before Cholesky factorization.
  C_star <- (
    C_star + t(C_star)
  )/2 + diag(g$jitter, nrow(C_star))
  
  #Return summary, with the Cholesky factor for later acquisition function usage
  list(
    mean = gp_predict(g, Xstar)$mean,
    covariance = C_star,
    chol = chol(C_star)
  )
}


#================== Gaussian Process MLII Estimation==================

#Constrain hyperparameters to be positive by a log-transform
.pack <- function(hp, naive = FALSE) {
  
  if (naive) {
    log(c(
      hp$ell,
      hp$sigma_f2,
      hp$sigma2
    ))
  } else {
    log(c(
      hp$ell,
      hp$sigma_f2,
      hp$sigma_b2,
      hp$sigma2
    ))
  }
}

#Reverse the log-transform post hyperparameter optimization, and introduce the 
#kernel into the returned object
.unpack <- function(p, D, kernel, naive = FALSE) {
  
  if (naive) {
    list(
      ell = exp(p[seq_len(D)]),
      sigma_f2 = exp(p[D + 1]),
      sigma_b2 = 1e-8,
      sigma2 = exp(p[D + 2]),
      kernel = kernel
    )
  } else {
    list(
      ell = exp(p[seq_len(D)]),
      sigma_f2 = exp(p[D + 1]),
      sigma_b2 = exp(p[D + 2]),
      sigma2 = exp(p[D + 3]),
      kernel = kernel
    )
  }
}


#Keep numerically boundary-valued warm starts within the box supplied to
#L-BFGS-B, while continuing to reject genuinely incompatible guesses.
.bound_mlii_start <- function(p, lower, upper, tolerance = 1e-8) {
  scale <- pmax(1, abs(lower), abs(upper))
  outside <- pmax(lower - p, p - upper, 0)

  if (
    any(!is.finite(p)) ||
    any(outside > tolerance * scale)
  ) {
    stop(
      'The supplied hyperparameter guess lies outside the current optimization bounds. ',
      'Update "lower" and/or "upper" if you want to use this starting point.',
      call. = FALSE
    )
  }

  pmin(pmax(p, lower), upper)
}


#Empirical Bayes Maximum (Marginal) Likelihood Estimation
gp_fit_mlii <- function(
    X,
    y,
    batch,
    kernel,
    naive = FALSE,
    n_starts = 8,
    hp_start = NULL,
    jitter = 1e-8,
    maxit = 200,
    lower = c(0.05, 1e-4, 1e-8, 1e-8),
    upper = c(50, 5, 5, 5)
) {
  X <- as.matrix(X)
  D <- ncol(X)


  controls <- .validate_mlii_controls(
    n_starts = n_starts,
    maxit = maxit,
    jitter = jitter,
    lower = lower,
    upper = upper
  )
  n_starts <- controls$n_starts
  maxit <- controls$maxit
  jitter <- controls$jitter
  lower <- controls$lower
  upper <- controls$upper
  
  #Define hyperparameter constraints
  if (naive) {
    
    lo <- log(c(
      rep(lower[1], D),  #The D length scales
      lower[2],          #Signal variance
      lower[4]           #Residual noise
    ))
    
    hi <- log(c(
      rep(upper[1], D),
      upper[2],
      upper[4]
    ))
    
  } else {
    
    lo <- log(c(
      rep(lower[1], D),  #The D length scales
      lower[2],          #Signal variance
      lower[3],          #Batch variance
      lower[4]           #Residual noise
    ))
    
    hi <- log(c(
      rep(upper[1], D),
      upper[2],
      upper[3],
      upper[4]
    ))
    
  }
  
  
  
  #This is the function optimized by optim() to find our GP hyperparameters.
  #For a candidate set of (log) hyperparameters, it returns the two quantities
  #optim() needs, namely the function value to be minimized (the NLL) and its gradient.

  cached_parameters <- NULL
  cached_evaluation <- NULL

  evaluate_mlii <- function(p) {
    #Check if the same parameters were proposed again, if so simply return those
    if (identical(as.numeric(p), cached_parameters)) return(cached_evaluation)

    #Return the NLL and gradient; if the fit fails, a large NLL and zero gradient.  
    evaluation <- tryCatch(
      {
        fit <- gp(
          X,
          y,
          batch,
          .unpack(p, D, kernel, naive),
          jitter
        )

        list( #Returned if fit succeeds
          value = fit$nll,
          gradient = .mlii_gradient(fit, naive = naive)
        )
      },
      error = function(e) {
        list( #Returned if fit fails
          value = 1e6,
          gradient = rep(0, length(p))
        )
      }
    )


    cached_parameters <<- as.numeric(p) #The <<- is a "super-assignment", updating variables one level above the function
    cached_evaluation <<- evaluation 
    evaluation
  }

  #Number of parameters optimized in the current model.
  P <- if (naive) D + 2 else D + 3
  
  #Draw the starting points: one per row, uniformly within the (log-scale) bounds.
  starts <- matrix(runif(n_starts * P, lo, hi), nrow = n_starts, byrow = TRUE)
  
  #If a hyperparameter guess is provided, use it as the first starting point.
  if (!is.null(hp_start)) starts[1, ] <- .bound_mlii_start(.pack(hp_start, naive = naive), lo, hi)
  
  best <- NULL
  best_value <- Inf
  records <- vector('list', n_starts) #Initialize fit diagnostics.
  
  
  
  for (i in seq_len(n_starts)){#For each candidate starting point, s_i:
    
    s <- starts[i,]
    initial_value <- evaluate_mlii(s)$value
    
    #Optimize the marginal likelihood from this start; if optim() itself fails,
    #record the unoptimized start with the error message
    fit <- tryCatch(
      stats::optim(
        s,
        fn = function(p) evaluate_mlii(p)$value,
        gr = function(p) evaluate_mlii(p)$gradient,
        method = 'L-BFGS-B',
        lower = lo,
        upper = hi,
        control = list(maxit = maxit)
      ),
      error = function(e) {
        list(
          par = s,
          value = initial_value,
          convergence = NA,
          message = conditionMessage(e)
        )
      }
    )
    
    #Record the results from each starting condition
    records[[i]] <- data.frame(
      start = i,
      initial_nll = initial_value,
      final_nll = fit$value,
      convergence = fit$convergence,
      message = .optimizer_message(fit$convergence, fit$message)
    )
    
    #Keep start result if it is the current best 
    if (isTRUE(fit$value < best_value)) {
      best_value <- fit$value
      best <- fit
    }
  }
  
  #Ensure that at least one start produced a valid fit (1e6 marks a failed evaluation).
  if (best_value >= 1e6) stop('ML-II fitting failed at every starting point.', call. = FALSE)
  
  #Stack the optimization diagnostics from all starting points into one data frame.
  diagnostics <- do.call(rbind, records)

  
  #Unpack the optimal hyperparameters.
  hp <- .unpack(best$par, D, kernel, naive = naive)

  parameter_names <- if (naive) {
    c(paste0("ell_", seq_len(D)), "sigma_f2", "sigma2")
  } else {
    c(paste0("ell_", seq_len(D)), "sigma_f2", "sigma_b2", "sigma2")
  }

  #Find parameters of the solution that lies on the optimization boundary
  tolerance <- 1e-6
  at_lower <- abs(best$par - lo) <= tolerance
  at_upper <- abs(best$par - hi) <= tolerance
  at_bound <- at_lower | at_upper #Check if a parameter is either at the lower or upper bound
  names(at_bound) <- parameter_names
  if (naive) at_bound["sigma_b2"] <- FALSE #Not estimated in the naive model
  
  #Append optimization results
  hp$nll <- best_value
  hp$fit_diagnostics <- diagnostics
  hp$convergence <- best$convergence
  hp$message <- .optimizer_message(best$convergence, best$message)
  hp$at_bound <- at_bound
  
  if (!isTRUE(best$convergence == 0)) {
    warning(
      "The selected ML-II solution did not report normal convergence: ",
      hp$message,
      call. = FALSE
    )
  }
  
  if (any(at_bound)) {
    warning(
      "The selected ML-II solution reached an optimization bound for: ",
      paste(names(at_bound)[at_bound], collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  
  hp
}
