pca_cov <- function(M.cov, threshold=0.95) {

  # add arguments to match with prcomp

  # derive PCs
  M.PCA <- prcomp(M.cov)

  # extract cumulative proportion of variance explained across all PCs
  cum.prop <- summary(M.PCA)$importance[3, ]

  # determine the minimum number of PCs needed to exceed the specified threshold
  PCs.subn <- length(cum.prop[cum.prop<=threshold]) + 1

  # extract PCs subset
  PCs.sub <- M.PCA$rotation[, 1:PCs.subn]

  # summary of PC results to output
  PCs.summ <- summary(M.PCA)$importance

  out <- list(PCs.sub = PCs.sub, PCs.summ = PCs.summ, PCs.all = M.PCA,
              minPCs = PCs.subn)

  return(out)

}
