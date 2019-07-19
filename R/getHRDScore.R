#' get the standard HRD score
#' 
#' @param seq.dat the preprocessed sequencing data
#' 
#' @details calculates the 'standard' HRD score: the sum of LST, LOH, and NTAI, all normalized. Can optionally
#' return the HRD score rescaled to a range of 0 - 100.
#' 
#' @return a numeric value representing the HRD score
#' 
#' @examples
#'seq.dat <- preprocessSeq( sub01.segments )
#'
#'CN.dat <- getCNT( seq.dat )
#'HRD <- getHRD.Score( seq.dat, CN.dat, sub01.ploidy, scaleTotal = FALSE)
#' @export

# ------------------------------- getHRD.Score ---------------------------------- #
# function to get the total HRD score. this is a simple summation of LST, NTAI, and LOH.
# this is the 'standard' HRD score most commonly seen in publication. this function is
# implemented to automatically calculate this score- eg, select the proper mean/raw criteria
getHRD.Score <- function( seq.dat, CN.dat, ploidy.dat, min.seg.size = 11e06, scaleTotal = FALSE )
  # input: seq.dat, (data.frame) with chromosome, start.pos, end.pos, CNt, alleleA, alleleB;
  #         ploidy.dat (data.frame), the ploidy data
  #         CN.dat (data.frame), copy number data
  #         min.seg.size (integer), minimum segment size used in TAI calculations
  #         scaleTotal (boolean), rescale HRD total to 0-100
  # output: HRD.Score (integer), the total HRD Score
{
  HRD.Score <- getLOH.norm(  seq.dat, ploidy.dat, 15e06 ) + getLST.norm(  seq.dat, ploidy.dat ) + getNTAI.norm( seq.dat, CN.dat, min.seg.size )
  
  # rescale HRD total to 0-100
  if( scaleTotal == TRUE )
  {
    library(scales)
    HRD.Score <- round( rescale(HRD.Score, to = c(0,100), from = range(HRD.Score)), 2 )
  }
  
  return(HRD.Score)
  
  
}
# ------------------------------------------------------------------------------- #
