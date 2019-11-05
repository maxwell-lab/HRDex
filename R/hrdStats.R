

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
hrd.stats <- function(seq.dat, ploidy.dat, CN.dat, min.seg.size = 6e06)
{
  
  seq.dat <- preprocessHRD(seq.dat)
  
  # raw data
  HRD.NTAIr <- getNTAI.raw( seq.dat, min.seg.size )
  HRD.LST  <- getLST(  seq.dat )
  HRD.LOH   <- getLOH(  seq.dat )
  HRD.NTAIm <- getNTAI.norm( seq.dat, CN.dat, ploidy.dat, min.seg.size )

  
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
