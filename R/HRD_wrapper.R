#!/usr/bin/env Rscript
library(HRDex)
library(argparse)

args = commandArgs(trailingOnly = TRUE)

if( length(args) < 2)
{
		stop("need to provide 1 argument: seq data file, reference genome (grch37 or grch38)")
}

seq.dat <- read.table(args[1], header = TRUE)
ref <- args[2]
sample.id <- args[3]

#column names for facets
#csv text or csv

seq.dat <- preprocessHRD( seq.dat, ref )#needs required fields
#Doies not check for A or B NA values, which facet may have.
CN.dat <- getCNt(seq.dat )
HRD <- getHRD.Score( seq.dat, CN.dat, scaleTotal = FALSE)

#HRD.NTAIr <- getNTAI.raw( seq.dat )
HRD.LST  <- getLST(  seq.dat )
HRD.LOH   <- getLOH(  seq.dat )
HRD.NTAIm <- getNTAI.norm( seq.dat, CN.dat )


cat(paste(sample.id,HRD.LOH,HRD.NTAIm,HRD.LST,HRD,sep=","))
