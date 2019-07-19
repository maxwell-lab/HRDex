
# functions to get telomeric allelic imbalance (TAI) without including main CNt segments
# both raw and ploidy corrected

# ------------------------------- getTAI.raw ----------------------------------------- #
getTAI.raw <- function(seq.dat, min.seg.size = 11e06)
  # input:
  #   seq.dat (data.frame), the sequencing data (eg, .seqz_segments.txt)
  #   min.seg.size, (integer), the minimum segment size required for analysis
  # output:
  #   HRD.TAIm (numeric)
{
  
  # seq.dat$s: length of the segment
  # if length of segment is greater than the minimum; and allelic imbalance is present; and
  # segment is in the telomeres; and is not on centromere...
  HRD.TAI <- sum(seq.dat$s > min.seg.size & seq.dat$AI & !seq.dat$cross.arm 
                 & (seq.dat$post.telomere | seq.dat$pre.telomere))
  
  return(HRD.TAI)
}
# ------------------------------------------------------------------------------- #


# ------------------------------- getTAI.norm ----------------------------------- #
getTAI.norm <- function(seq.dat, CN.dat, min.seg.size = 11e06)
  # input:
  #   seq.dat (data.frame), the sequencing data (eg, .seqz_segments.txt)
  #   min.seg.size, (integer), the minimum segment size required for analysis
  # output:
  #   HRD.TAIm (numeric)
{
  
  # seq.dat$s: length of the segment
  # if length of segment is greater than the minimum; and allelic imbalance is present; and
  # segment is in the telomeres; and is not on centromere...
  
  # create an index of main.CN segments; these get removed to normalize TAI
  rm.ind <- c()
  
  for( i in 1:length(CN.dat$chromosome))
  {
    rm.ind <- c(rm.ind, which(seq.dat$chromosome == CN.dat$chromosome[i] & seq.dat$CNt == CN.dat$main.CN[i]))
  }
  
  # normalized
  HRD.TAI <- getTAI.raw(seq.dat[-rm.ind,], 11e06)
  
  return(HRD.TAI)
}
# ------------------------------------------------------------------------------- #