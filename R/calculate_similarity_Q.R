calculate_similarity_Q_stat <- function(
    SE,
    cor_mat,
    M,
    obs.ests,
    Q.stat.out = TRUE,
    alpha.nom = 0.05,
    Bonf.corr = 1,
    Bonf.corr.strict = 26,
    Hotelling.corr = TRUE,
    cov.df = NULL
) {

  Bonf.corr <- sort(Bonf.corr)

  E <- length(SE)

  if (E < 2) {
    stop("At least two estimators are required.")
  }

  if (length(obs.ests) != E) {
    stop("The length of obs.ests must equal the length of SE.")
  }

  if (!all(dim(cor_mat) == c(E, E))) {
    stop("cor_mat must be an E by E correlation matrix.")
  }

  if (ncol(M) != E) {
    stop("The number of columns in M must equal the length of SE.")
  }

  if (any(!is.finite(SE)) || any(SE <= 0)) {
    stop("All standard errors must be positive and finite.")
  }

  # By default, assume that C was estimated using the rows of M
  if (is.null(cov.df)) {
    cov.df <- nrow(M) - 1
  }

  # Strict Bonferroni threshold
  alpha.Bonf.strict <- alpha.nom / Bonf.corr.strict



  S <- diag(as.numeric(SE), nrow = E, ncol = E)

  C <- S %*% cor_mat %*% S


  if (rcond(C) < 1e-12) {
    stop("The covariance matrix is singular or numerically unstable.")
  }


  one <- rep(1, E)

  C.inv.one <- solve(C, one)

  GLS.denominator <- sum(C.inv.one)

  calculate_GLS_mean <- function(estimates) {

    estimates <- as.numeric(estimates)

    # C^{-1} beta
    C.inv.estimates <- solve(C, estimates)

    # (1'C^{-1}beta) / (1'C^{-1}1)
    sum(C.inv.estimates) / GLS.denominator
  }


  Ests <- as.numeric(obs.ests)

  MU <- calculate_GLS_mean(Ests)

  d <- Ests - MU

  Q.1 <- as.numeric(
    crossprod(d, solve(C, d))
  )


  Q.df <- E - 1


  if (Hotelling.corr) {

    if (!is.numeric(cov.df) ||
        length(cov.df) != 1 ||
        !is.finite(cov.df) ||
        cov.df < Q.df) {
      stop(
        paste0(
          "cov.df must be at least E - 1. ",
          "Usually cov.df = number of bootstrap samples - 1."
        )
      )
    }

    Hotelling.df2 <- cov.df - Q.df + 1

    Hotelling.scale <-
      Hotelling.df2 / (Q.df * cov.df)

    F.1 <- Hotelling.scale * Q.1

    Q.1.pval <- stats::pf(
      F.1,
      df1 = Q.df,
      df2 = Hotelling.df2,
      lower.tail = FALSE
    )


  } else {

    F.1 <- NA_real_
    Hotelling.df2 <- NA_real_

    Q.1.pval <- stats::pchisq(
      Q.1,
      df = Q.df,
      lower.tail = FALSE
    )
  }

  p.sim.vec <- Q.1.pval

  Q.sim.df <- NULL

  if (Q.stat.out) {

    if (Hotelling.corr) {

      Q.sim.df <- data.frame(
        Q.stat = Q.1,
        F.stat = F.1,
        df1 = Q.df,
        df2 = Hotelling.df2,
        p.sim = p.sim.vec
      )

    } else {

      Q.sim.df <- data.frame(
        Q.stat = Q.1,
        p.sim = p.sim.vec
      )
    }

    rownames(Q.sim.df) <- "Q statistics"
  }


  Bonf.corr.M <- matrix(
    NA,
    nrow = 1,
    ncol = length(Bonf.corr)
  )

  colnames(Bonf.corr.M) <- paste0("PC", Bonf.corr)
  rownames(Bonf.corr.M) <- "Q statistics"

  for (i in seq_along(Bonf.corr)) {
    Bonf.corr.M[, i] <-
      p.sim.vec < (alpha.nom / Bonf.corr[i])
  }

  alpha.nom.vec <- p.sim.vec < alpha.nom

  Bonf.corr.strict.vec <-
    p.sim.vec < alpha.Bonf.strict

  out <- list(
    p.sim.vec = p.sim.vec,
    Bonf.corr.vec = Bonf.corr.M,
    alpha.nom.vec = alpha.nom.vec,
    Bonf.corr.strict.vec = Bonf.corr.strict.vec,
    Q.sim.df = Q.sim.df
  )

  return(out)
}
