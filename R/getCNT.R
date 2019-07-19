#' get the main copynumber state for each chromosome
#' 
#' @param seq.dat the data.frame of sequencing data
#' 
#' @details Get the main copynumber state for each chromosome. Calculate the total proportion of
#' the chromosome associated with each unique copynumber state; return the copynumber state that is most
#' represented.
#' 
#' @return a data.frame containing the main copynumber state for each input chromosome
#' 
#' @examples 
#' seq.dat <- preprocessSeq( seq.dat )
#' CN.dat <- getCNt( seq.dat )
#' 
#' @export


# ----------------------------------- getCNT ------------------------------------------ #
# setup CNt values
getCNt <- function( seq.dat )
  # input:
  #   seq.dat (data.frame), the sequencing data (eg, .seqz_segments.txt)
  # out:
  #   CN.dat (data.frame), then number
{
  main.CN <- rep(0, length(unique(seq.dat$chromosome)))
  
  ct1 <- 0
  for( chr in unique(seq.dat$chromosome))
  {
    ct1 <- ct1 + 1
    
    # reduce data to one chromosome at a time
    dat <- seq.dat[seq.dat$chromosome == chr,]
    
    # remove any NAs or this will crash
    dat <- dat[!is.na(dat$CNt),]
    CNt.frac <- rep(0, length(unique(dat$CNt)))
    ct2 <- 0
    
    # get the fraction of the chromosome association with each copynumber state
    for( CNt in unique(dat$CNt) )
    {
      ct2 <- ct2 + 1
      CNt.frac[ct2] <- sum(dat$frac.chr[dat$CNt == CNt])
    }
    
    # record values + associated CNt
    out <- data.frame(CNt=unique(dat$CNt), CNt.frac=CNt.frac)
    
    # get copynumber state most prominent in the chromosome
    main.CN[ct1] <- out$CNt[which(CNt.frac == max(CNt.frac))]
    
  }
  
  CN.dat <- data.frame(chromosome=unique(seq.dat$chromosome), main.CN=main.CN)
  return(CN.dat)
}
# ------------------------------------------------------------------------------- #
