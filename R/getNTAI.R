#' Compute the number of non-telomeric allelic imbalance (NTAI) events
#' 
#' @param seq.dat the data.frame of sequencing data
#' @param min.seg.size the minimum segment size

#' 
#' @details raw NTAI is calculated as the number of segments where the segment size is greater than the 
#' minimum segment size; allelic imbalance is present; the segment does not cross the centromere; and the
#' segment is not in either of the telomeres. 
#' @return the number of NTAI events
#' 
#' @examples 
#' seq.dat <- preprocessHRD( seq.dat )
#' CN.dat <- getCNt( seq.dat )
#' 
#' ## the number of NTAI events in seq.dat
#' ntai <- getNTAI.raw( seq.dat )
#' 
#' @export


# --------------------------------------------------------------------------------- #
getNTAI.raw <- function(seq.dat, min.seg.size = 11e06)
  # input:
  #   seq.dat (data.frame), the sequencing data (eg, .seqz_segments.txt)
  #   min.seg.size, (integer), the minimum segment size required for analysis
  # output:
  #   HRD.NTAIm (numeric)
{
  # if length of segment is greater than the minimum; and allelic imbalance is present; and
  # segment is not in the telomeres or centromere...
  
  if( !(all( c("AI", "cross.arm", "post.telomere", "pre.telomere") %in% colnames(seq.dat))))
  {
    stop("some required fields missing. run 'preprocessHRD' first.")
  }
  
  HRD.NTAI <- sum(seq.dat$seg.len > min.seg.size & 
                    seq.dat$AI & 
                    !seq.dat$cross.arm & 
                    !seq.dat$post.telomere & 
                    !seq.dat$pre.telomere)
  
  return(HRD.NTAI)
}
# --------------------------------------------------------------------------------- #



#' Compute the number of non-telomeric allelic imbalance (NTAI) events
#' 
#' @param seq.dat the data.frame of sequencing data
#' @param CN.dat the data.frame of copy number information
#' @param min.seg.size the minimum segment size
#' 
#' @details raw NTAI is calculated as the number of segments where the segment size is greater than the 
#' minimum segment size; allelic imbalance is present; the segment does not cross the centromere; and the
#' segment is not in either of the telomeres. NTAI is normalized by removing all main copy number segments.
#' 
#' @return the number of NTAI events
#' 
#' @examples 
#' seq.dat <- preprocessSeq( seq.dat )
#' CN.dat <- getCNt( seq.dat )
#' 
#' ## the number of NTAI events in seq.dat
#' ntai <- getNTAI.norm( seq.dat, CN.dat)
#' 
#' @export


# ------------------------------------------------------------------------------ #
getNTAI.norm <- function(seq.dat, CN.dat, min.seg.size = 11e06 )
  # input:
  #   seq.dat (data.frame), the sequencing data (eg, .seqz_segments.txt)
  #   min.seg.size, (integer), the minimum segment size required for analysis
  # output:
  #   HRD.NTAIm (numeric)
{
  # create an index of main.CN segments; these get removed to normalize TAI
  rm.ind <- c()
  
  for( i in 1:length(CN.dat$chromosome))
  {
    rm.ind <- c(rm.ind, which(seq.dat$chromosome == CN.dat$chromosome[i] & seq.dat$CNt == CN.dat$main.CN[i]))
  }
  
  # normalize
  HRD.NTAI <- getNTAI.raw(seq.dat[-rm.ind,], min.seg.size)
  
  return(HRD.NTAI)
}
# --------------------------------------------------------------------------------- #