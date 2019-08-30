

#' Compute the number of loss of heterozygosity (LOH) events
#'
#' @param seq.dat a data.frame of the sequencing data
#' @param min.seg.size minimum size of the seqgment
#' @param ploidy the estimated ploidy
#' 
#' @details raw LOH is calculated as the number of events where: the segment length is greater than the minimum segment size 
#' (default of 15e06), and the segment size is less than 90% of the length of the entire chromosome.

#' @export

# ------------------------------- getLOH --------------------------------------- #
# function to get LOH 
getLOH <- function(seq.dat, min.seg.size = 15e06)
  # input:
  #   seq.dat (data.frame), the sequencing data (eg, .seqz_segments.txt)
  # output:
  #   HRD.LOH (integer), number of LOH events
{
  # HRD-LOH calculation
  # loss of heterozygosity
  # if B == 0 & s > 15mbp & within chromosome, then HRD-LOH is TRUE
  HRD.LOH <- sum( (seq.dat$seg.len > min.seg.size) & ((seq.dat$seg.len / seq.dat$chr.size) < 0.9) & 
                    (seq.dat$B == 0) )
  
  return(HRD.LOH)
}
# --------------------------------------------------------------------------------- #

