#' Example dataset
#'
#' This summary statistics aims to estimate the causal effect of low-density lipoprotein (LDL) cholesterol on coronary artery disease (CAD) risk, using eight SNPs with strong biological
#' links to LDL cholesterol.
#'
#' @format A data frame with 100 rows and 8 variables:
#' \describe{
#'   \item{rsID}{SNP identifier}
#'   \item{beta.exposure}{Beta effect on exposure}
#'   \item{se.exposure}{Standard error of beta.exposure}
#'   \item{beta.outcome}{Beta effect on outcome}
#'   \item{se.outcome}{Standard error of beta.outcome}
#'   \item{id.exposure}{Exposure identifier}
#'   \item{id.outcome}{Outcome identifier}
#'   \item{mr_keep}{Whether to include this instrument}
#' }
#'
#' @source Available at: https://pmc.ncbi.nlm.nih.gov/articles/PMC6124628/
"LDL_CAD"
