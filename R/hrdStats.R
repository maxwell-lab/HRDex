#'
#' 
#' @param seq.dat the preprocessed sequencing data
#' @param CN.dat the copy number data
#' @param ploidy.dat the ploidy data
#' @param min.seg.size optional: the minimum segment size (default to 11Mbp)
#' @param scaleTotal optional: rescale the total to range of 0-100
#' @param type optional: type of score to return, one of "sum" or "average"
#' 
#' 
#' @return LOH, NTAIraw, NTAInorm, LST, total
#' 
#' @examples
#'seq.dat <- hrd.stats(  )
#'
#' @export

# ------------------------------------- hrd.stats ------------------------------------- #
# hrd.stats is a function to compute the three HRD metrics (HRD-LOH, HRD-NTAI, and
# HRD-LST), as well as total HRD and mean HRD. this simply wraps up the output in one
# dataframe- you can run all of these steps individually for error checking.

#
# input: seq.dat, (data.frame) with chromosome, start.pos, end.pos, CNt, alleleA, alleleB;
#         ploidy.dat (data.frame), the ploidy data
#         CN.dat (data.frame), copy number data
#         min.seg.size (integer), minimum segment size used in TAI calculations
#         scaleTotal (boolean), rescale HRD total to 0-100
# output: out, a data.frame with HRD metrics
hrd.stats <- function(seq.dat, CN.dat, ref = "grch37")
{
  
  seq.dat <- preprocessHRD(seq.dat, ref)
  
  # use default min seg size for each call
  HRD.NTAIr <- getNTAI.raw( seq.dat )
  HRD.LST  <- getLST(  seq.dat )
  HRD.LOH   <- getLOH(  seq.dat )
  HRD.NTAIm <- getNTAI.norm( seq.dat, CN.dat, ploidy.dat )

  
  out = data.frame(
                   HRD.LOH  = HRD.LOH,
                   HRD.NTAIr = HRD.NTAIr,
                   HRD.NTAIm = HRD.NTAIm,
                   HRD.LST  = HRD.LST
                   )
  
  # HRD.Score can be the total of any 3 metrics, raw or norm- this is the 'standard' score
  out$HRD.Score <- getHRD.Score( seq.dat, CN.dat, ploidy.dat, min.seg.size, scaleTotal = FALSE )
  
  return(out)
  
}
# ------------------------------------------------------------------------------- #
