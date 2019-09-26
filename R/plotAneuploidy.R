
#' plot aneuploidy


# -------------------------- plotAneuploidy -------------------------------------- #
# function to make a plot of the aneuploid data
# input: dat (data.frame), the seq data that has been run through getAneuploidy, so it
#         contains aneuploid status. this can be either a single chromosome, or multiple
#
# output: p1 (ggplot.object), a ggplot object, plotted outside of this function so image parameters
# can be adjusted
plotAneuploidy <- function( dat, ref = "grch37", X.include = FALSE )
{
  i.int = 1
  
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
 
  
  chr.labs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
                "chr20", "chr21", "chr22", "chrX")


  if( X.include == FALSE)
  {
    ref.dat <- ref.dat[1:22,]
    chr.labs <- chr.labs[1:22]
  }
  
  levels(dat$chromosome) <- levels(ref.dat$chromosome)
  levels(dat$a.stat) <- c("Amplifed", "Neutral", "Deleted")
  # if there is more than one chromosome, start/end positions of the segments and reference points
  # need to be adjusted to their point along whole genome (as opposed to their point in the chromosome)
  # this adds the appropriate offset
  # ---
  ref.tmp <- ref.dat
  chr.size = ref.dat$chr.size
  
  
  if(length(unique(dat$chromosome)) > 1)
  {
    chr.lab.pos <- rep(0, dim(ref.dat)[1])
    chr.lab.pos[1] <- chr.size[1] / 2
    
    for(i in ref.tmp$chromosome)
    {
      
      if( i == "chrX" )
      {
        i.int <- 23
      } else
      {
        i.int <- as.integer(gsub("chr", "", i))
      }
      
      # reference coordinates are relative to the start of the chromosome
      # this code changes the reference point for ALL chromosomes to be the
      # start of the genome
      if(i.int > 1)
      {
        dat$start.pos[dat$chromosome == i] <- dat$start.pos[dat$chromosome == i] + sum(ref.tmp$chr.size[1:(i.int-1)])
        dat$end.pos[dat$chromosome == i] <- dat$end.pos[dat$chromosome == i] + sum(ref.tmp$chr.size[1:(i.int-1)])
        ref.tmp$centromere.start[ref.tmp$chromosome == i] <- ref.tmp$centromere.start[ref.tmp$chromosome == i] + sum(ref.tmp$chr.size[1:(i.int-1)])
        ref.tmp$centromere.end[ref.tmp$chromosome == i] <- ref.tmp$centromere.end[ref.tmp$chromosome == i] + sum(ref.tmp$chr.size[1:(i.int-1)])
        chr.size[i.int] <- sum(ref.tmp$chr.size[1:i.int])
        chr.lab.pos[i.int] <- sum(ref.tmp$chr.size[1:(i.int - 1)]) + (ref.tmp$chr.size[i.int] / 2)
      }
      
    }
  }
  # ---
  
  y <- rep(0, dim(dat)[1])
  y[ dat$a.stat == "Amplified" ] <- 1
  y[ dat$a.stat == "Deleted" ] <- -1
  
  if( length(unique(dat$chromosome)) > 1) 
  {
    chr = ref.tmp$chromosome
  } else
  {
    chr = unique(dat$chromosome)
  }
  

  p1 <- ggplot( ) + 
    geom_segment(data = dat, 
                 aes( x = dat$start.pos, y = y, xend = dat$end.pos, yend = y, colour = a.stat), size = 3) +
    scale_y_continuous(breaks = seq(from = -4, to = 4, by = 0.25), limits = c(-4,4)) +
    scale_color_manual(values = c("green", "blue", "red")) +
    geom_vline(xintercept = ref.tmp$centromere.start[ref.tmp$chromosome %in% chr], linetype = "dashed") +
    geom_vline(xintercept = ref.tmp$centromere.end[ref.tmp$chromosome %in% chr], linetype = "dashed") +
    xlab("Genomic Position") +
    ylab("Status") +
    ggtitle( "Aneuploidy") +
    labs(color = "Alteration Status", drop = FALSE) +
    theme_classic( ) +
    theme( axis.text.y = element_blank(), axis.ticks.y = element_blank() ) 
    
#    geom_vline(xintercept = chr.size[1:23], color = "red", linetype = "dashed") +
    if( i.int > 1)
    {
      p1 <- p1 + geom_rect(aes(xmin = 0, xmax = chr.size[1], ymin = -Inf, ymax = Inf), alpha = 0.1, fill = "blue") +
        geom_rect(aes(xmin = chr.size[2], xmax = chr.size[3], ymin = -Inf, ymax = Inf), alpha = 0.1, fill = "blue") +
        geom_rect(aes(xmin = chr.size[4], xmax = chr.size[5], ymin = -Inf, ymax = Inf), alpha = 0.1, fill = "blue") +
        geom_rect(aes(xmin = chr.size[6], xmax = chr.size[7], ymin = -Inf, ymax = Inf), alpha = 0.1, fill = "blue") +
        geom_rect(aes(xmin = chr.size[8], xmax = chr.size[9], ymin = -Inf, ymax = Inf), alpha = 0.1, fill = "blue") +
        geom_rect(aes(xmin = chr.size[10], xmax = chr.size[11], ymin = -Inf, ymax = Inf), alpha = 0.1, fill = "blue") +
        geom_rect(aes(xmin = chr.size[12], xmax = chr.size[13], ymin = -Inf, ymax = Inf), alpha = 0.1, fill = "blue") +
        geom_rect(aes(xmin = chr.size[14], xmax = chr.size[15], ymin = -Inf, ymax = Inf), alpha = 0.1, fill = "blue") +
        geom_rect(aes(xmin = chr.size[16], xmax = chr.size[17], ymin = -Inf, ymax = Inf), alpha = 0.1, fill = "blue") +
        geom_rect(aes(xmin = chr.size[18], xmax = chr.size[19], ymin = -Inf, ymax = Inf), alpha = 0.1, fill = "blue") +
        geom_rect(aes(xmin = chr.size[20], xmax = chr.size[21], ymin = -Inf, ymax = Inf), alpha = 0.1, fill = "blue")
      p1 <- p1 + scale_x_discrete(drop = FALSE, limits = chr.lab.pos[ ref.tmp$chromosome %in% chr ], 
                                  labels = chr.labs )
    }
    
    if( X.include == TRUE)
    {
      p1 <- p1 + geom_rect(aes(xmin = chr.size[22], xmax = chr.size[23], ymin = -Inf, ymax = Inf), alpha = 0.1, fill = "blue")
    }
    
     #  geom_vline(xintercept = ref.tmp$chr.size, color = "white", linetype = "dashed")
    # scale_fill_discrete(drop = FALSE) +
    if( i.int == 1)
    {
      p1 <- p1 + xlim(0, ref.tmp$chr.size[ref.tmp$chromosome == chr])
    }
  
  return(p1)
 
}
# ------------------------------------------------------------------------------- #
