MR_F_bar <- function(beta.exposure, se.exposure) {
  #' Calculate mean F statistic
  #'
  #' Calculate mean F statistic to capture the dilution in the inverse variance weighted estimate for a set of instrument-exposure effects.
  #' @param beta.exposure A vector of instrument-exposure association estimates.
  #' @param se.exposure A vector of standard error estimates for the set of instrument-exposure association estimates.
  #'
  #' @return Returns the F bar statistic value.
  #' @export
  #'
  #' @examples
  #' F.bar <- MR_F_bar(
  #'   beta.exposure = LDL_CAD$beta.exposure,
  #'   se.exposure = LDL_CAD$se.exposure
  #' )
}
