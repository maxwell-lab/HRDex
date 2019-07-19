

#' Compute the number of loss of heterozygosity (LOH) events
#'
#' @param seq.dat a data.frame of the sequencing data
#' @param min.seg.size minimum size of the seqgment
#' @param ploidy the estimated ploidy
#' 
#' @details raw LOH is calculated as the number of events where: the segment length is greater than the minimum segment size 
#' (default of 15e06), and the segment size is less than 90% of the length of the entire chromosome.
#' Normalized LOH is calculated as raw LOH * k, the ploidy correction factor
#' @seealso \code{\link{ploidyCorrectionFactor}}
#' @export

# ------------------------------- getLOH --------------------------------------- #
# function to get LOH 
getLOH.raw <- function(seq.dat, min.seg.size = 15e06)
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


# ------------------------------- getLOH --------------------------------------- #
# function to get LOH 
getLOH.norm <- function(seq.dat, ploidy.dat, min.seg.size = 15e06)
  # input:
  #   seq.dat (data.frame), the sequencing data (eg, .seqz_segments.txt)
  # output:
  #   HRD.LOH (numeric), number of LOH events corrected for ploidy
{
  # HRD-LOH calculation
  # loss of heterozygosity
  # if B == 0 & s > 15mbp & within chromosome, then HRD-LOH is TRUE
  HRD.LOH <- getLOH.raw(seq.dat, min.seg.size) * ploidyCorrectionFactor(ploidy.dat$ploidy.estimate[2])

  return(HRD.LOH)
}
# --------------------------------------------------------------------------------- #
