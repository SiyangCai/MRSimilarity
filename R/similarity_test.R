#' Perform similarity test on the GWAS dataset with multiple testing.
#'
#' Given exposure and outcome beta effects, this function calculated the estimates using given listed MR methods from \code{TwoSampleMR} package.
#' An overall hypothesis, followed by group-wise comparisons or pair-wise comparisons as multiple testing will be performed.
#'
#'
#' @param data Dataframe of input GWAS dataset.
#' @param title_text The title of scatter plot of instrument effects and causal effect estimates.
#' @param xlab The name of x-axis of scatter plot of instrument effects and causal effect estimates.
#' @param ylab The name of y-axis of scatter plot of instrument effects and causal effect estimates.
#' @param n.boot The number of iteration of bootstrap to generate correlation matrix of causal estimates.
#' @param parametric_bootstrap Whether to use parametric or non-parametric bootstrap to generate correlation matrix of estimates. The default is FALSE.
#' @param groupwise Whether in multiple testing, group-wise comparisons are performed. The default is FALSE to only generate pair-wise comparisons.
#' @param MR.cols Column names that are required for exposure and outcome. If they are not included, this function will generate by itself.
#' @param methods A list of names of available MR methods that are elements from \code{TwoSampleMR:mr_method_list()}. The default includes
#' Inverse variance weighted (IVW), MR Egger (MRE), Weighted median (WM), Weighted mode (MBE) and MR GRIP (MG).
#' @param use_pca_gao Whether to use PCA on scenarios to find significant levels for multiple testing, as an improved method of Bonferroni correction. If using default as FALSE, all significant levels are nominal level of 0.05.
#' @param seed Random seed for reproducibility.
#' @param output.M.boot Whether to include the bootstrap estimates matrix in the final output. The default is FALSE.
#'
#' @import plotly
#' @import gtools
#'
#' @return A list of results and diagnostics contains:
#'   \itemize{
#'     \item{\code{res} The MR estimates results from listed methods.}
#'     \item{\code{test.stats} Test statistics for presence of heterogenity and weak instrument bias.}
#'     \item{\code{corr} The correlation matrix of causal estimates from bootstrap.}
#'     \item{\code{cov.mat} The rescaled covariance of causal estimates using the standard errors from original methods.}
#'     \item{\code{plot} scatter plot of instrument effects and causal effect estimates.}
#'     \item{\code{Q.stats.df} The similarity Q-statistics for all comparison scenarios.}
#'     \item{\code{Q.pvals.df} The p-values of similarity Q-statistics for all comparison scenarios.}
#'     \item{\code{Q.pvals.Bonf.df} Whether the similarity Q-statistics are significant on the signifcant level using adjusted Bonferroni.}
#'     \item{\code{alpha.Bonf} The final significant level when \code{use_pca_gao} is enabled.}
#'   }
#'
#' @examples
#' data(LDL_CAD)
#'
#' # n.boot = 100 is used here for quick demonstration.
#' # Increase n.boot for more reliable results in real analyses.
#' result <- MR_similarity_analysis(LDL_CAD, n.boot = 100)
#' print(result)
#'
#'
#'
#' @export


MR_similarity_analysis <- function(data,
                                   title_text = "Genetic instrument associations of exposure and outcome",
                                   xlab = "SNP-exposure association",
                                   ylab = "SNP-outcome association",
                                   n.boot = 1000,
                                   parametric_bootstrap = FALSE,
                                   groupwise = FALSE,
                                   MR.cols = c("id.exposure", "id.outcome",
                                               "exposure", "outcome",
                                               "mr_keep"),
                                   methods = c("Inverse variance weighted",
                                               "MR Egger", "Weighted median",
                                               "Weighted mode", "MR GRIP"),
                                   use_pca_gao = FALSE,
                                   seed = 42,
                                   PC.threshold = 0.99,
                                   output.M.boot = FALSE) {

  # set seed for reproducibility
  set.seed(seed)

  # mapping list for each MR.col default
  MR.cols.def <- list(id.exposure = "exposure", id.outcome = "outcome",
                      exposure = "exposure", outcome = "outcome",
                      mr_keep = TRUE)

  MR.tab = data

  for (c in MR.cols) {
    if (c %in% colnames(MR.tab) == FALSE) {
      MR.tab[, c] <- MR.cols.def[[c]]
    }
  }


  beta.exposure = MR.tab$beta.exposure
  beta.outcome = MR.tab$beta.outcome
  se.exposure = MR.tab$se.exposure
  se.outcome = MR.tab$se.outcome

  # flip the sign of BetaYG if BetaXG has been flipped
  beta.outcome = beta.outcome * sign(beta.exposure)
  # map BetaXG onto the positive plane
  beta.exposure = abs(beta.exposure)

  # number of genetic instruments
  L <- length(beta.exposure)

  # method calculation
  # first set up for default four estimators only then expand to an automated
  # set up with all TWoSampleMR estimators

  # quick plan: compute mr(d) then extract the rows of interest, and use the
  # MR-Egger specific function to extract the intercept if need be

  # iterate through each method name and convert into the function string
  # stored in TwoSampleMR::mr_method_list
  # convert the obj column to a list then assign the name string using names
  mr_methods <- TwoSampleMR::mr_method_list()

  mr_methods.list <- as.list(mr_methods$obj)
  names(mr_methods.list) <- mr_methods$name

  methods.vec <- NULL
  for (m in methods) {
    methods.vec <- c(methods.vec, mr_methods.list[[m]])
  }

  print("MR METHODS FUNCTIONS USED:")
  print(methods.vec)

  # later plan: see if it's quicker to compute reduced MR functions for each
  # method once you've got the input arguments figured out
  res <- TwoSampleMR::mr(MR.tab, method_list = methods.vec)

  print("RESULTS GENERATED")


  n.est <- nrow(res)
  beta.intercepts <- rep(0, n.est)

  Egger.int.p = NULL
  GRIP.int.p = NULL

  i = 0
  for (m in res$method) {
    i = i + 1
    if (m == "MR Egger") {
      resEgger <- TwoSampleMR::mr_egger_regression(MR.tab$beta.exposure,
                                                   MR.tab$beta.outcome,
                                                   MR.tab$se.exposure,
                                                   MR.tab$se.outcome)
      beta.intercepts[i] = resEgger$b_i
      # store p-value of intercept term to output with other statistics
      Egger.int.p <- resEgger$pval_i
    }

    if (m == "MR GRIP") {
      resGRIP <- TwoSampleMR::mr_grip(MR.tab$beta.exposure,
                                      MR.tab$beta.outcome,
                                      MR.tab$se.exposure,
                                      MR.tab$se.outcome)
      beta.intercepts[i] = resGRIP$b_i
      # store p-value of intercept term to output with other statistics
      GRIP.int.p <- resGRIP$pval_i
    }
  }


  # append the intercepts onto the results dataframe
  res$beta.intercepts = beta.intercepts

  print("INTERCEPTS ADDED")

  p <- MR_scatter_estimates(beta.exposure = beta.exposure,
                          beta.outcome = beta.outcome,
                          se.exposure = se.exposure, se.outcome = se.outcome,
                          b = res$b, beta.intercepts = beta.intercepts,
                          title_text = title_text, xlab = xlab, ylab = ylab,
                          methods = res$method)

  print("PLOT CREATED")


  # Q statistics
  het_stats <- TwoSampleMR::mr_heterogeneity(MR.tab)

  # Cochran
  Cochran.Q <- het_stats$Q[het_stats$method == "Inverse variance weighted"]
  Cochran.Q.p <- het_stats$Q_pval[het_stats$method == "Inverse variance weighted"]

  # Rucker
  Rucker.Q <- het_stats$Q[het_stats$method == "MR Egger"]

  # difference p-value
  Q.diff.p <- 1 - pchisq(Cochran.Q - Rucker.Q, 1)

  # dilution statistics
  F.bar <- MR_F_bar(beta.exposure = beta.exposure, se.exposure = se.exposure)

  # MR-Egger dilution
  Isq_XG <- TwoSampleMR::Isq(y = beta.exposure, s = se.exposure)


  params.boot <- TwoSampleMR::default_parameters()
  params.boot$nboot <- 2

  M.boot1 <- matrix(NA, nrow = n.boot, ncol = n.est)
  for (i in c(1:n.boot)) {

    if(parametric_bootstrap == FALSE){
      S = sample(L, replace=TRUE)
      MR.tab.boot1 <- MR.tab[S, ]
    }
    else{
      MR.tab.boot1  = MR.tab

      MR.tab.boot1$beta.exposure = rnorm(length(MR.tab$beta.exposure), mean = MR.tab$beta.exposure, sd = MR.tab$se.exposure)
      MR.tab.boot1$beta.outcome = rnorm(length(MR.tab$beta.outcome), mean = MR.tab$beta.outcome, sd = MR.tab$se.outcome)
    }


    # generate necessary estimates
    res.boot1 <- suppressMessages(TwoSampleMR::mr(MR.tab.boot1,
                                                  method_list = methods.vec,
                                                  parameters = params.boot))
    # store results in bootstrap matrix
    M.boot1[i, ] <- res.boot1$b

    # extract the method names for the bootstrap matrix
    if (i == 1) {
      M.boot.names <- res.boot1$method
      colnames(M.boot1) <- M.boot.names
    }
  }

  M.boot2 <- matrix(NA, nrow = n.boot, ncol = n.est)
  for (i in c(1:n.boot)) {

    if(parametric_bootstrap == FALSE){
      S = sample(L, replace=TRUE)
      MR.tab.boot2 <- MR.tab[S, ]
    }
    else{
      MR.tab.boot2  = MR.tab

      MR.tab.boot2$beta.exposure = rnorm(length(MR.tab$beta.exposure), mean = MR.tab$beta.exposure, sd = MR.tab$se.exposure)
      MR.tab.boot2$beta.outcome = rnorm(length(MR.tab$beta.outcome), mean = MR.tab$beta.outcome, sd = MR.tab$se.outcome)
    }


    # generate necessary estimates
    res.boot2 <- suppressMessages(TwoSampleMR::mr(MR.tab.boot2,
                                                  method_list = methods.vec,
                                                  parameters = params.boot))
    # store results in bootstrap matrix
    M.boot2[i, ] <- res.boot2$b

    # extract the method names for the bootstrap matrix
    if (i == 1) {
      M.boot.names <- res.boot2$method
      colnames(M.boot2) <- M.boot.names
    }
  }

  print("BOOTSRAP SAMPLING COMPLETE")


  # compute covariance matrix
  C.boot1 = cov(M.boot2)
  colnames(C.boot1) = rownames(C.boot1) = M.boot.names

  # compute correlation matrix for the output objects
  C1 = cor(M.boot2)

  # calculation

  ests.M <- matrix(nrow = 1, ncol = n.est)
  colnames(ests.M) <- M.boot.names

  for (m in res$method) {
    ests.M[1, m] <- res$b[res$method == m]
  }

  SE = res$se
  names(SE) = M.boot.names



  # compute variance covariance matrix for the output objects
  cov.mat = diag(as.numeric(SE)) %*% C1 %*% diag(as.numeric(SE))

  # implement PCA before calculating Q_sim p-values to report
  Q.sim.PCA <- PCA_scenarios(SE, M.boot = M.boot1, ests.M = ests.M, groupwise = groupwise,
                                         threshold = PC.threshold)

  # extract Bonferroni correction ascertained from PCA of Q_sim cov matrix
  Bonf.corr <- Q.sim.PCA$Bonf.corr

  # extract cumulative variance results across the full PC space
  PCs.summ <- Q.sim.PCA$C.PCs$PCs.summ



  # construct all the combinations of the valid estimators
  # create a list to store the combination matrix of each j value
  combs <- NULL
  # create a dummy object to store each set of similarity values
  Q.pvals.full <- NULL
  # create a dummy vector of names to store the names of each combination
  Q.pvals.names <- NULL
  # create a dummy vector to store the bools of whether each subset has a p_sim
  # value below the Bonferroni threshold
  Q.pvals.Bonf <- NULL
  # dummt object to store Q_sim statistics
  Q.stats.full <- NULL
  for (j in 2:n.est) {
    # stack combination sets for each j
    combs[[j]] <- gtools::combinations(n.est, j, res$method)
    # iterate through each matrix and compute the corresponding similarity
    # statistics
    for (k in 1:nrow(combs[[j]])) {

      if(groupwise == FALSE){
        if(j != 2 && j != n.est){
          next
        }
      }

      ind.vals <- combs[[j]][k, ]

      Q.pvals.res <- calculate_similarity_Q_stat(SE[ind.vals],
                                             C1[ind.vals, ind.vals],
                                             M.boot1[, ind.vals],
                                             ests.M[1, ind.vals],
                                             Bonf.corr = Bonf.corr,
                                             Q.stat.out = TRUE,
                                             Bonf.corr.strict = ifelse(groupwise, sum(choose(n.est, 2:n.est)), choose(n.est,2)))

      # extract p-values
      Q.pvals <- Q.pvals.res$p.sim.vec
      comb.name <- paste(combs[[j]][k, ], collapse = "--")

      Q.pvals.names <- c(Q.pvals.names, comb.name)
      Q.pvals.full <- rbind(Q.pvals.full, Q.pvals)
      Q.pvals.Bonf <- rbind(Q.pvals.Bonf, t(Q.pvals.res$Bonf.corr.vec))
      Q.stats <- t(Q.pvals.res$Q.sim.df$Q.stat)
      Q.stats.full <- rbind(Q.stats.full, Q.stats)
    }
  }

  # output string to summarise findings of PC
  PC.summ.text <- paste(Bonf.corr, " out of a possible ",
                        length(Q.pvals.names),
                        " principal components were needed to explain ",
                        100*PC.threshold,
                        "% of the cumulative variance within the Q_sim variance-covariance matrix. This yields an adjusted significance threshold of ",
                        round(0.05/Bonf.corr, 3),
                        sep = "")

  print("SIMILARITY STATISTIC CALCULATED")

  Q.pvals.df <- data.frame(round(Q.pvals.full, 4))
  rownames(Q.pvals.df) <- Q.pvals.names
  colnames(Q.pvals.df) <- c("Q statistics")

  Q.stats.df <- round(Q.stats.full, 2)
  rownames(Q.stats.df) <- Q.pvals.names
  colnames(Q.stats.df) <- c("Q statistics")

  if(use_pca_gao == TRUE){
    Q.pvals.Bonf.df <- data.frame(round(Q.pvals.Bonf, 4))
    rownames(Q.pvals.Bonf.df) <- Q.pvals.names
    colnames(Q.pvals.Bonf.df) <- c("Q statistics")
  }
  else{
    Q.pvals.Bonf.df = Q.pvals.df
    Q.pvals.Bonf.df$`Q statistics` = as.integer(Q.pvals.Bonf.df$`Q statistics` < 0.05)
  }

  # collect test statistics into a single df
  test.stats <- c("Cochran's Q", "Cochran's Q (p)", "Cochran-Rucker diff. (p)",
                  "F bar", "I^2_XG", "MR-Egger intercept (p)")
  values <- rbind(Cochran.Q, Cochran.Q.p, Q.diff.p, F.bar, Isq_XG, Egger.int.p)
  test.stats.df <- data.frame(test.stat = test.stats, value = round(values, 4))
  rownames(test.stats.df) <- NULL

  print("TEST STATISTICS DF CREATED")

  output <- list(res = res, test.stats = test.stats.df, corr = C1, cov.mat = cov.mat, plot = p,
                 Q.stats.df = Q.stats.df, Q.pvals.df = Q.pvals.df, Q.pvals.Bonf.df = Q.pvals.Bonf.df,
                 alpha.Bonf = 0.05/Bonf.corr)

  if (output.M.boot) {
    output$M.boot = M.boot1
  }

  return(output)
}
