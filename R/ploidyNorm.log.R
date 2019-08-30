
#' Correct HRD metric for ploidy assuming a log-linear relationship
#' 
#' @param HRD vector of numeric values, the HRD metric to normalize. Can be HRD.Score or any individual measure.
#' @param ploidy vector of ploidy values.
#' 
#' @details Regress out the effect of ploidy from HRD 
#' 
#' @return logHRD corrected for ploidy
#' 
#' @export

ploidyNorm.log <- function( HRD, ploidy )
{
  fit <- lm( log(HRD) ~ ploidy )
  return(fit$residuals)
}