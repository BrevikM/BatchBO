################################################################################
# Gradients.R -- analytic derivatives for GP hyperparameter optimization
################################################################################


#================== Radial-Kernel Derivative Factor ==================

#For the implemented kernels, the derivative of the correlation with respect to the
#length-scale-normalized Euclidean distance d can be written in the form -d h(d).
.kernel_gradient_factor <- function(distance, kernel) {
  
  switch( #Return h(d)
    kernel,
    matern52 = {
      sqrt5_distance <- sqrt(5)*distance
      (5/3)*(1+sqrt5_distance)*exp(-sqrt5_distance)
    },
    matern32 = {
      3*exp(-sqrt(3)*distance)
    },
    rbf = {
      exp(-0.5*distance^2)
    }
  )
}


#================== Correlation Length-Scale Derivatives ==================

#Derivative of the correlation matrix R with respect to each log length scale.
#dR_ij / dlog(ell_k) = h(d_ij) (Delta_ijk / ell_k)^2
.corr_log_ell_derivative <- function(X, ell, kernel) {
  X <- as.matrix(X)
  n <- nrow(X)
  D <- ncol(X)
  
  #Scaled differences in each factor, Delta_ijk / ell_k, summed into the squared
  #Euclidean distance d_ij^2 = sum_k (Delta_ijk / ell_k)^2
  scaled_difference <- array(0, dim = c(n, n, D))
  squared_distance <- matrix(0, n, n)
  for (k in seq_len(D)) {
    scaled_difference[, , k] <- outer(X[, k], X[, k], "-") / ell[k]
    squared_distance <- squared_distance + scaled_difference[, , k]^2
  }
  
  h <- .kernel_gradient_factor(sqrt(squared_distance), kernel) #h(d_ij) for every pair
  
  #One derivative matrix per factor: h(d_ij) (Delta_ijk / ell_k)^2
  dR <- array(0, dim = c(n, n, D))
  for (k in seq_len(D)) {
    dR[, , k] <- h * scaled_difference[, , k]^2
  }
  
  dR
}


#================== Profiled ML-II Gradient ==================

#Differentiate the profiled negative log marginal likelihood with respect to
#the log hyperparameters used by .pack() and .unpack(). The constant mean needs
#no separate derivative because it is profiled at its GLS stationary point.
.mlii_gradient <- function(g, naive = FALSE) {
  X <- as.matrix(g$X)
  hp <- g$hp
  D <- ncol(X)

  V_inverse <- chol2inv(g$U)
  trace_weight <- V_inverse - tcrossprod(g$a)

  correlation <- .corr(
    X,
    X,
    hp$ell,
    hp$kernel
  )
  correlation_derivative <- .corr_log_ell_derivative(
    X,
    hp$ell,
    hp$kernel
  )

  length_scale_gradient <- vapply(
    seq_len(D),
    function(k) {
      0.5*hp$sigma_f2*sum(#Find the trace sum
        trace_weight*correlation_derivative[, , k]
      )
    },
    numeric(1)
  )

  signal_gradient <- 0.5*hp$sigma_f2*sum(trace_weight*correlation)
  residual_gradient <- 0.5*hp$sigma2*sum(diag(trace_weight))
  batch_gradient <- 0.5*sum(trace_weight*hp$sigma_b2*.S(g$batch))

  if (naive) {
    return(c(
      length_scale_gradient,
      signal_gradient,
      residual_gradient
    ))
  }

  c(
    length_scale_gradient,
    signal_gradient,
    batch_gradient,
    residual_gradient
  )
}
