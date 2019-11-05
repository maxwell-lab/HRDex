
#' Combine two chromosomal segments into one
#'
#' @param seq.dat a data.frame containing the sequencing data
#' @param max.brk.len only combine segments where the distance between the two is less than max.brk.len (in mbp)
#' @param ap.calc set to true if calculating aneuploidy, which uses different criteria than HRD
#' 
#' @details For HRD: if CNt, A, and B are equal in two adjacent segments, and break length < max.brk.len, combine
#' the two segments into one. For aneuploidy; if aneuploidy status is equal in two adjacent segments, combine into one.
#' @return a truncated data frame with the combined segments
#' @examples
#' ## There are 4 segments on chromosome 14
#' seq.dat <- sub01.segments[ sub01.segments$chromosome == "chr14",]
#' seq.dat
#' 
#' ## The gap between segments 3 and 4 is very small; copynumber, A, and B are identical in both segments
#' ## this is a false split. It will be combined into one segment. 
#' seq.dat <- combineSeg(seq.dat)
#' seq.dat
#' 
#' ## the data now has 3 segments.
#' @export

# ------------------------------------------------------------------------------------- #
# some long segments are falsely read as multiple shorter segments
# recombine these into one segment
combineSeg <- function(seq.dat, max.brk.len = 3e06, ap.calc = FALSE)
  # input:
  #   seq.dat (data.frame), the sequencing data (eg, .seqz_segments.txt)
  #   max.brk.len (integer), the maximum size of the break (gap between two segments)
  #       combine segments where break.len < max.brk.len
  #   ap.calc (logical), are you calculating aneuploidy?
  # rules for combining segments are different for aneuploidy and HRD-
  # for HRD: CNt, A, and B must be equal in two adjacent alleles with break length < 3Mbp to combine
  # for aneuploidy, no limit on breaklength, and only aneuploidy needs to be equal in adjacent segments 
  # to combine
  #
  # output:
  #   seq.dat (data.frame), the sequencing data with segments joined and centromere segments removed  
  
{
  seq.dat$brk.len <- 0
  
  n.segs <- dim(seq.dat)[1]
  rm.ind <- c()
  
  for( i in 1:(n.segs - 1))
  {
    seq.dat$brk.len[i + 1] <- seq.dat$start.pos[i + 1] - seq.dat$end.pos[i]
    
    # for aneuploidy, join two adjacent segments if ap status is the same and gap is less than max break length
    if( ap.calc == TRUE )
    {
      if( seq.dat$a.stat[i + 1] == seq.dat$a.stat[i] & 
          seq.dat$brk.len[i + 1] > 0 &
          seq.dat$brk.len[i + 1] < max.brk.len)
      {
        seq.dat$start.pos[i + 1] <- seq.dat$start.pos[i]
        rm.ind <- c(rm.ind, i)
      }
      
    # hrd segment joining
    } else 
      {
        # if break length is < max.brk.len; & CNt1 == CNt2 & A1 == A2 & B1 == B2
        # then combine the two segments
        if( seq.dat$brk.len[i + 1] < max.brk.len & 
            seq.dat$brk.len[i + 1] > 0 &
            seq.dat$CNt[i + 1] == seq.dat$CNt[i] &
              seq.dat$A[i + 1] == seq.dat$A[i] &
                seq.dat$B[i + 1] == seq.dat$B[i]   )
      
        {
          seq.dat$start.pos[i + 1] <- seq.dat$start.pos[i]
          rm.ind <- c(rm.ind, i)
        }
      }
  }
  
  # remove one of the segments from each pair that was joined
  if( !is.null(rm.ind) )
  {
    seq.dat <- seq.dat[-rm.ind,]
  }
 
  return(seq.dat)
  
  
}
# ------------------------------------------------------------------------------------- #


