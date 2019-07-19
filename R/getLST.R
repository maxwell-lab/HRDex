#' Compute the number of large state transition (LST) events
#' 
#' @param seq.dat the data.frame of sequencing data
#' @param ploidy.dat the data.frame of ploidy data
#' 
#' @details raw LST is calculated as the number of segments where the gap between is < 3Mbp, each adjacent 
#' segment is > 10Mbp, and the segments do not cross the centromere. NTAI is normalized by k, the ploidy
#' correction factor
#' 
#' @return the number of LST events
#' 
#' @export

# functions to get large state transition (LST) events, both raw and corrected for ploidy

# ---------------------------------- getLST.raw ----------------------------------- #
getLST.raw <- function(seq.dat)
  # input:
  #   seq.dat (data.frame), the sequencing data (eg, .seqz_segments.txt)
  # output:
  #   HRD.LST (numeric), number of LST events
{
  
  # length of gap between segments
  seq.dat$brk.len <- 0
  seq.dat$LST <- FALSE
  
  
  n.segs <- dim(seq.dat)[1]
  
  # if the gap between 2 segments is < 3mbp, and each adjacent segment is > 10mbp, and the segment does not
  # cross the centromere, it is an LST
  for( i in 1:(n.segs - 1))
  {
    seq.dat$brk.len[i + 1] <- seq.dat$start.pos[i + 1] - seq.dat$end.pos[i]
    seq.dat$LST[i] <- seq.dat$brk.len[i] < 3e06 &
      seq.dat$cross.arm[i] == FALSE & 
      seq.dat$seg.len[i] > 10e06 & 
      seq.dat$seg.len[i + 1] > 10e06
  }
  
  HRD.LST <- sum(seq.dat$LST)
  return(HRD.LST)
  
}
# ------------------------------------------------------------------------------- #



# ---------------------------------- getLST.norm --------------------------------------- #
getLST.norm <- function(seq.dat, ploidy.dat)
  # input:
  #   seq.dat (data.frame), the sequencing data (eg, .seqz_segments.txt)
  # output:
  #   HRD.LST (numeric), large state transition
{
  
  HRD.LSTm = getLST.raw(seq.dat) * ploidyCorrectionFactor(ploidy.dat$ploidy.estimate[2])
  return(HRD.LSTm)
  
}
# ------------------------------------------------------------------------------- #
