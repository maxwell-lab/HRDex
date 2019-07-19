getLOH.seqlengths <- function( seq.dat, min.seg.size = 15e06 )
{
  s <-  seq.dat$s[ (seq.dat$s > min.seg.size) & ((seq.dat$s / seq.dat$chr.size) < 0.9) & 
                       (seq.dat$B == 0) ] 
  if( length(s) == 0 )
  {
    return(0)
  }
  
  return(s)
}

getNTAI.seqlengths <- function( seq.dat, CN.dat, min.seg.size = 11e06 )
{
  rm.ind <- c()
  
  for( i in 1:length(CN.dat$chromosome))
  {
    rm.ind <- c(rm.ind, which(seq.dat$chromosome == CN.dat$chromosome[i] & seq.dat$CNt == CN.dat$main.CN[i]))
  }
  
  seq.dat <- seq.dat[-rm.ind,]
  s <- seq.dat$s[ seq.dat$s > min.seg.size & 
                      seq.dat$AI & 
                      !seq.dat$cross.arm & 
                      !seq.dat$post.telomere & 
                      !seq.dat$pre.telomere ] 
  if( length(s) == 0 )
  {
    return(0)
  }
  
  return(s)
}