
#================== Expected Improvement ==================
.ei <- function(mu, s, fbest){
  s <- pmax(s, 1e-12) #Avoid potential zero-division
  z <- (mu - fbest)/s
  
  (mu-fbest)*pnorm(z) + s*dnorm(z)
}



#================== Kriging-Believer prediction ==================
#Posterior of the batch-free surface at the points in X, after treating the
#fantasies as measured at their predicted values. Returns an unchanged mean and
#the reduced variance s_j^2(x) = s^2(x) - c_j(x)^T C_j^{-1} c_j(x).

#The fantasy quantities do not depend on the candidate, so they are computed
#once per batch member rather than at every evaluation of the acquisition.
.kb_fantasy_state <- function(g, fantasies) {
  if (!nrow(fantasies)) {
    return(NULL)
  }

  #Cholesky factor of the posterior covariance among the exact fantasy locations.
  fantasy_chol <- gp_joint(g, fantasies)$chol

  #V^{-1/2} K(X, fantasies), reused for every candidate's posterior cross-covariance.
  fantasy_training <- forwardsolve(
    t(g$U),
    g$hp$sigma_f2 * .corr(g$X, fantasies, g$hp$ell, g$hp$kernel)
  )

  list(fantasies = fantasies, chol = fantasy_chol, training = fantasy_training)
}


.kb_predict <- function(g, X, fantasy_state) {
  prediction <- gp_predict(g, X)

  #Before any fantasy has been selected, return the ordinary GP prediction.
  if (is.null(fantasy_state)) {
    return(prediction)
  }

  #Posterior covariance between candidates and fantasy locations,
  #C(X, F) = K(X, F) - K(X, X_n) V^{-1} K(X_n, F).
  candidate_training <- forwardsolve(
    t(g$U),
    g$hp$sigma_f2 * .corr(g$X, X, g$hp$ell, g$hp$kernel)
  )
  candidate_fantasy_covariance <- g$hp$sigma_f2 * .corr(
    X, fantasy_state$fantasies, g$hp$ell, g$hp$kernel
  ) - crossprod(candidate_training, fantasy_state$training)

  #Condition the candidate variance on the exact fantasy observations.
  W <- forwardsolve(
    t(fantasy_state$chol),
    t(candidate_fantasy_covariance)
  )

  prediction$var <- pmax(
    prediction$var - colSums(W^2),
    0
  )

  prediction
}


acq_kb <- function(
    g,
    q,
    n_screen = 4096,
    restarts = 12,
    maxit = 500,
    lower = NULL,
    upper = NULL,
    seed = NULL,
    max_screen_factor = 4,
    X_baseline = NULL
) {
  D <- ncol(g$X)


  #Measured conditions defining the incumbent; by default every distinct training location
  if (is.null(X_baseline)) X_baseline <- unique(g$X)
  
  state <- list(
    fantasies = matrix(numeric(0), nrow = 0, ncol = D),
    fbest = max(gp_predict(g, X_baseline)$mean)
  )
  
  .greedy_cont(
    g = g,
    q = q,
    state = state,
    
    #EI with the unchanged mean, the variance reduced by the fantasies, and f*
    make_acquisition = function(state, j) {
      fantasy_state <- .kb_fantasy_state(g, state$fantasies)
      function(X) {
        prediction <- .kb_predict(g = g, X = X, fantasy_state = fantasy_state)
        .ei(mu = prediction$mean, s = sqrt(prediction$var), fbest = state$fbest)
      }
    },
    
    #Believe the prediction at the new pick: add it as a fantasy and raise f* if it is higher
    update_state = function(state, x, j) {
      Xnew <- matrix(x, nrow = 1)
      state$fantasies <- rbind(state$fantasies, Xnew)
      state$fbest <- max(state$fbest, gp_predict(g, Xnew)$mean)
      state
    },
    
    n_screen = n_screen,
    restarts = restarts,
    maxit = maxit,
    lower = lower,
    upper = upper,
    seed = seed,
    max_screen_factor = max_screen_factor
  )
}
