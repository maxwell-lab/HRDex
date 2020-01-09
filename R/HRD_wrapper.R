#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

if( length(args) < 2)
{
		stop("need to provide 1 argument: seq data file, reference genome (grch37 or grch38)")
}

seq.dat <- read.table(args[1], header = TRUE)
ref <- args[2]

seq.dat <- preprocessHRD( seq.dat, ref )
CN.dat <- getCNT(seq.dat )
HRD <- getHRD.Score( seq.dat, CN.dat, scaleTotal = FALSE)

print(HRD)
