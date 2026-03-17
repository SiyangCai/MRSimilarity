#'
#'@export


calculate_similarity_Q_stat <- function(SE,
                                    cor_mat,
                                    M,
                                    obs.ests,
                                    Q.stat.out = TRUE,
                                    alpha.nom = 0.05,
                                    Bonf.corr = 1,
                                    Bonf.corr.strict) {
  ###
  # Compute both forms of the similarity statistic from a given covariance
  # matrix.
  ###

  Bonf.corr <- sort(Bonf.corr)

  # compute the corrected threshold under the Bonferroni method (taken as 0.05
  # unless specified otherwise)
  alpha.Bonf <- alpha.nom / Bonf.corr
  # strict Bonferroni correction that is overly conservative
  alpha.Bonf.strict <- alpha.nom / Bonf.corr.strict
  # compute number of estimators captured in C
  E <- length(SE)

  # rescaled covariance matrix
  SE = diag(as.numeric(SE))
  C = SE %*% cor_mat %*% SE

  # variance of each estimate type ; convert into first order weights
  V = diag(C)
  W = 1 / V
  # store effect estimates in Ests
  # store mean bootstrap estimates in Ests2
  Ests = Ests2 = matrix(nrow = 1, ncol = E)
  # store the observed effect estimates
  Ests[1,] = obs.ests
  Ests2[1,] = apply(M, 2, mean)
  # weighted mean of effect estimates
  MU = sum(W * Ests) / sum(W)
  # bias of estimates (beta_j - beta_hat)
  d = Ests - MU

  # matrix multiplication to calculate the generalised Q-statistic
  Q.1 = d %*% solve(C) %*% t(d)
  Q.1.pval <-  1 - pchisq(Q.1, E - 1)

  Q.sim.vec = Q.1
  p.sim.vec = Q.1.pval

  # if Q statistics aren't requested a null object will be compiled for output
  Q.sim.df <- NULL

  if (Q.stat.out) {
    Q.sim.df <- data.frame(Q.stat = Q.sim.vec, p.sim = p.sim.vec)
    rownames(Q.sim.df) <- c("Q statistics")
  }

  # compile a matrix of vectors to store boolean results of nPCs thresholds
  Bonf.corr.M <- matrix(NA, ncol = length(Bonf.corr), nrow = 1)

  colnames(Bonf.corr.M) =
    paste("PC", as.character(Bonf.corr), sep = "")
  rownames(Bonf.corr.M) <- c("Q statistics")

  for (i in Bonf.corr) {
    Bonf.ind <- paste("PC", i, sep = "")
    Bonf.corr.M[, Bonf.ind] <- (p.sim.vec < alpha.nom / i)

  }

  alpha.nom.vec <- (p.sim.vec < alpha.nom)
  Bonf.corr.strict.vec <- (p.sim.vec < alpha.Bonf.strict)

  out <-
    list(
      p.sim.vec = p.sim.vec,
      Bonf.corr.vec = Bonf.corr.M,
      alpha.nom.vec = alpha.nom.vec,
      Bonf.corr.strict.vec = Bonf.corr.strict.vec,
      Q.sim.df = Q.sim.df
    )

  return(out)

}
