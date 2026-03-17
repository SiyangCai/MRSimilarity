PCA_scenarios <- function(SE, M.boot, ests.M, threshold = 0.95) {

  # extract size of bootstrap matrix
  n.boot <- nrow(M.boot)
  n.est <- ncol(M.boot)
  # set up objects needed to run MR_similarity_analysis Q_sim loops

  # create a list to store the combination matrix of each j value
  combs <- NULL

  # compute covariance matrix of bootstrap estimates - will be treated as
  # constant within PCA implementation
  C.boot = cov(M.boot)

  # new objects needed for PCA implementation
  # list to store Q_simi values for all combination of estimates
  Q.stat.i <- NULL


  # add in loop to calculate the variance covariance matrix of the different
  # Q_sim statistics

  for (j in 2:n.est) {

    if(j != 2 && j != n.est){
      next
    }

    # stack combination sets for each j
    combs[[j]] <- gtools::combinations(n.est, j, colnames(M.boot))
    # iterate through each matrix and compute the corresponding similarity
    # statistics
    for (k in 1:nrow(combs[[j]])) {
      ind.vals <- combs[[j]][k, ]
      print(ind.vals)
      comb.name.full <- paste(combs[[j]][k, ], collapse = "--")

      ## TEMP CODING FOR TRAINING - NEED TO EXPAND TO ALL POSSIBLE INPUTS
      # replace full estimator names with prefix
      # for full package store a list that can be easily accessed
      comb.name <- gsub("Inverse variance weighted", "IVW", comb.name.full)
      comb.name <- gsub("MR Egger", "MRE", comb.name)
      comb.name <- gsub("Weighted median", "WM", comb.name)
      comb.name <- gsub("Weighted mode", "MBE", comb.name)

      print("###########################")
      print(comb.name)
      print("###########################")

      # objects will later form an output data frame of Q_sim values for each
      # subset
      #Q.pvals.names <- c(Q.pvals.names, comb.name)
      #Q.pvals.full <- rbind(Q.pvals.full, Q.pvals)

      # reset object used for n.boot Q_simi values derived bu subsetting through
      # M.boot
      Q.stat.i.sub <- NULL
      for (l in 1:n.boot) {
        # select new causal effect estimates
        #print("PRENAMING")
        ests.M.sub <- matrix(nrow = 1, ncol = n.est)
        #print("ests.M.sub created")
        #print(l)
        ests.M.sub <- as.matrix(M.boot[l, ])
        dim(ests.M.sub) <- c(1, n.est)
        #print(ests.M.sub)
        colnames(ests.M.sub) <- colnames(M.boot)
        #print("colnames assigned")
        Q.sim.res <- calculate_similarity_Q_stat(SE[ind.vals], cor(M.boot)[ind.vals, ind.vals],
                                             M.boot[, ind.vals],
                                             ests.M.sub[1, ind.vals],
                                             Q.stat.out = TRUE)
        #Q.pvals.sub <- Q.sim.res$p.sim
        # only interested in Q_sim1 form for this area
        #print(Q.sim.res)
        Q.stat.i.sub[l] <- Q.sim.res$Q.sim.df["Q.1", "Q.stat"]
      }


      # name variance-covariance matrix using comb.name

      # think about how to iterate through and compute covariance calculations

      # populate Q.stat.i for each subset combination of estimators
      Q.stat.i[[comb.name]] <- Q.stat.i.sub
      # CHECKPOINTS
      #print(k)
      #print(comb.name)
      #print(Q.stat.i)
      print("END OF COMBINATION")
      print("------------------------------------------------------------------")
    }

    # populate Q sim covariance matrix
    #print(Q.stat.i)
    #print(Q.sim.res)
    #print(Q.stat.i)
    C.Q.cov <- cov_Q(Q.stat.i)
  }

  # apply PCA to covariance matrix - use a custom PCA function built from JB code
  C.PCs <- pca_cov(C.Q.cov, threshold = threshold)

  # Bonf.corr <- ncol(C.PCs$PCs.sub)
  Bonf.corr = C.PCs$minPCs

  output <- list(Q.stat.i = Q.stat.i, C.Q.cov = C.Q.cov, C.PCs = C.PCs,
                 Bonf.corr = Bonf.corr)

  return(output)

}
