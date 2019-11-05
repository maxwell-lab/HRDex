#' Compute the aneuploidy score for a given chromosome
#' 
#' @param seq.dat the data.frame of sequencing data
#' @param ploidy.dat the data.frame of ploidy data
#' @param chr the chromosome of interest
#' 
#' @details Aneuploidy status is defined as amplified if main copynumber > ploidy; neutral if main copy number and ploidy
#' are equal; and deleted if main copy number < ploidy. Segments are combined if the aneuploidy status of two adjacent
#' segments is equal, and the segments do not cross the centromere.
#' 
#' @seealso \code{\link{combineSeg}}
#' @export

# ------------------------------ getAneuploidy ------------------------------------------ #
getAneuploidy <- function ( seq.dat, ploidy.dat, chr, ref = "grch37", max.brk.len = 3e06 )
# input: seq.dat, (data.frame) with chromosome, start.pos, end.pos, CNt, alleleA, alleleB;
#        ploidy.dat (data.frame), the ploidy data
#        chr (integer or character), the chromosome of interest
#
# output: out (data.frame), p and q arm aneuploid status and chr
#         dat (data.frame), the seq data updated with aneuploid information, this is used
#           for plotting
# dont use preprocess seg- we dont want to combine or eliminate anything
  
{
  if( ref == "grch37" )
  {
    ref.dat <- grch37.ref.dat
  } else
    if( ref == "grch38" )
    {
      ref.dat <- grch38.ref.dat
    } else
    {
      print(paste(ref, "is not a valid reference genome.", sep = " "))
      stop("select one of: grch37, grch38")
    }
  
  # check that chromosome is formatted as "chr1"
  if( substr(chr, 1, 1) != "c")
  {
    chr <- paste("chr", chr, sep = "")
  }
  
  if( substr(seq.dat$chromosome[1], 1, 1) != "c")
  {
    seq.dat$chromosome <- paste("chr", seq.dat$chromosome, sep = "")
  }
  
  if( nchar(chr) > 5 | nchar(chr) < 4)
  {
    stop("incorrect chr format")
  }
  
  dat <- seq.dat[ seq.dat$chromosome == chr, ]
  
  # brk.len is a variable used in HRD calculation
  # it isnt (currently) used in aneuploidy but the combineSeg function
  # returns breaklength as part of the df. need to attach it to the main data structure here
  # to prevent throwing a warning
  if( !("brk.len" %in% colnames(dat)))
  {
    dat$brk.len <- 0
  }
  
  if( dim(dat)[1] == 0 )
  {
    stop("failed to subset by chr")
  }
  
  ploidy <- round(ploidy.dat$ploidy.estimate[2])
  
  # centromere for this chromosome
  ref.cent.start <- ref.dat$centromere.start[ ref.dat$chromosome == chr]
  ref.cent.end <- ref.dat$centromere.end[ ref.dat$chromosome == chr]
  
  if( chr %in% c("chr13", "chr14", "chr15", "chr21", "chr22"))
  {
    dat$arm <- "c"
    dat$arm[ !in.centromere( dat$start.pos, dat$end.pos, ref.cent.start, ref.cent.end )] <- "p"
  } else
  {
    dat$arm = "c"
    dat$arm[dat$end.pos <= ref.cent.start] <- "p"
    dat$arm[dat$start.pos >= ref.cent.end] <- "q"
  }
  
  # remove any segments that cross into the centromere
  ind <- mapply(in.centromere, dat$start.pos, dat$end.pos, ref.cent.start, ref.cent.end)
  if( any(ind))
  {
    dat <- dat[!ind,]
  }
  
  if( dim(dat)[1] == 0 )
  {
    print("no segments remaining")
    return(list(NULL, NULL))
  }
  
  # setup aneuploidy status
  dat$a.stat <- "na"
  dat$a.stat <- as.factor(dat$a.stat)
  levels(dat$a.stat) <- c("Amplified", "Neutral", "Deleted")
  
  dat$a.stat[ dat$CNt > ploidy] <- levels(dat$a.stat)[1]
  dat$a.stat[ dat$CNt == ploidy] <- levels(dat$a.stat)[2]
  dat$a.stat[ dat$CNt < ploidy] <- levels(dat$a.stat)[3]
  

  # combine adjacent segments IFF anuploidy status is the same in both
  if( dim(dat)[1] > 1)
  {
    i = 1
    repeat {
      if( dat$a.stat[i] == dat$a.stat[i + 1])
      {
        # Inf means no limit on break length- we ignore it in this case
        tmp <- combineSeg( dat[i:(i+1),], max.brk.len, ap.calc = TRUE)
        
        # dont join segments across centromere
        if( !in.centromere(tmp$start.pos, tmp$end.pos, ref.cent.start, ref.cent.end))
        {
          dat[i + 1,] <- tmp
          dat <- dat[-i,]
        } else
        {
          # move on to the next segment
          i = i + 1
        }
      } else
      {
        i = i + 1
      }
      
      # end loop when at the last segment
      if( i == dim(dat)[1])
      {
        break
      }
    }
  }
 
 
  p.arm.len <- ref.cent.start
  q.arm.len <- ref.dat$chr.size[ref.dat$chromosome == chr] - ref.cent.end
  dat$s <- dat$end.pos - dat$start.pos
  
  p.aneuploid = q.aneuploid = FALSE
  
  # p arm
  # some chrs have no P arm
  p.dat <- dat[dat$end.pos <= ref.cent.start,]
  

  
  
  if( dim(p.dat)[1] > 0)
  {
    # only look at long arm for chr13, 14, 15, 21, 22
    if( !(chr %in% c("chr13", "chr14", "chr15", "chr21", "chr22")))
    {
      if( dim(p.dat)[1] == 0)
      {
        print(paste("no segments in p-arm of ", chr, sep = ""))
      } else
      {
        if( any( p.dat$a.stat %in% c("Amplified", "Deleted")))
        {
          
          len <- p.dat$s[ which(p.dat$s == max(p.dat$s[p.dat$a.stat %in% c("Amplified", "Deleted")])) ]
          
          # length of the altered region must be > 80% of the total arm to count as an alteration
          if( (len / p.arm.len) > .8) { p.aneuploid = TRUE }
        }
        
        
      }
    }
  }
  
  # q arm
  q.dat <- dat[dat$start.pos >= ref.cent.end,]
  
  if( dim(q.dat)[1] > 0)
  {
    if( any( q.dat$a.stat %in% c("Amplified", "Deleted")))
    {
      len <- q.dat$s[ which(q.dat$s == max(q.dat$s[q.dat$a.stat %in% c("Amplified", "Deleted")])) ]
      if( (len / q.arm.len) > .8) { q.aneuploid = TRUE }
    }
  }

  out <- data.frame( p.aneuploid = p.aneuploid, q.aneuploid = q.aneuploid, chr = chr)
  return(list(out, dat))
}
# ------------------------------------------------------------------------------- #



#' 
#' @param seq.dat the data.frame of sequencing data
#' @param ploidy.dat the data.frame of ploidy data
#' @param include.X, include chromosome (logical
#' 
#' @details get aneuploidy status across the whole genome
#' 
#' @seealso \code{\link{getAneuploidy}}
#' @export
# ------------------------------- getAneuploidyGenome --------------------------- #
# function to get aneuploidy scores for the whole genome
getAneuploidyGenome <- function( seq.dat, ploidy.dat, include.X = FALSE )
# indput: seq.dat (data.frame), the raw sequencing data
#         ploidy.dat (data.frame), the ploidy data
# output: list of
#         out1 (data.frame), p and q arm aneuploid status
#         out2 (data.frame), the seq data updated with aneuploid information, this is used
#           for plotting
{
  out1 <- c()
  out2 <- c()
  
  chrs <- seq(1:22)
  
  if( include.X == TRUE )
  {
    chrs <- c(chrs, "X")
  }
  
  for( i in chrs)
  {
    
    x = getAneuploidy( seq.dat, ploidy.dat, i)
    out1 <- rbind( out1, x[[1]] )
    out2 <- rbind( out2, x[[2]] )
  }
  
  return(list(out1, out2))
}
# ------------------------------------------------------------------------------- #



#' 
#' @param dat, data.frame containing aneuploidy scores
#'
#' @details description goes here
#' 
#' @seealso \code{\link{getAneuploidy}}
#' @export
# ---------------------------- getAneuploidyScore ------------------------------- #
# function to get total aneuploidy score (out of 39)
getAneuploidyScore <- function( dat )
  # input: dat (data.frame), seq data that been run through getAneuploidy
  # output: score (integer), score in the range of 0-39
{
  if( any(!( c("p.aneuploid", "q.aneuploid") %in% colnames(dat))))
  {
    print("missing p.aneuploid and q.aneuploid columns in input data")
    print("use the summary data- first element of the output of getAneuploidy")
    stop("missing input")
  }
  return( sum( dat$p.aneuploid == TRUE ) + sum( dat$q.aneuploid == TRUE) )
}
# ------------------------------------------------------------------------------- #
