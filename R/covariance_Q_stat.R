cov_Q <- function(Q.sim.stats) {
  # variance-covariance matrix of Q_sim statistics for each combination subset
  # 4 estimator case assumed for purposes of this but could automate using
  # choose(n, n-1) + ... choose(n.est, n-2)
  n.sub <- length(Q.sim.stats)
  C.Q.sim <- matrix(NA, nrow = n.sub, ncol = n.sub)
  # name matrix using the combination names used in Q.sim.stats
  colnames(C.Q.sim) = rownames(C.Q.sim) = names(Q.sim.stats)

  # iterate through each element of the matrix and compute the covariance
  # this can likely be made more efficient in future code through apply funcs
  for (s in 1:n.sub) {
    sub.name <- colnames(C.Q.sim)[s]
    # this ran before with original hard coding and no use of lapply
    #print(paste("SUBSET:", sub.name))

    C.Q.sim[s, ] <- cov(Q.sim.stats[[s]], data.frame(Q.sim.stats))
    #print("################################################################")

  }
  return(C.Q.sim)
}
