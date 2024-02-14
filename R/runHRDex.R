#!/usr/bin/env Rscript
library(HRDex)
library(argparse)

p <- ArgumentParser()
p$add_argument('-i', '--infile', help = 'Path to input BED file.')
p$add_argument('-o', '--outfile', help = "Path to output file.")
p$add_argument('--tumor', help = 'Tumor id.')
p$add_argument('--build', help = 'Genome build; {grch37,grch38}.')

args <- p$parse_args()

# Reading the BED file and specifying column names
seq.dat <- read.table(args$infile, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(seq.dat) <- c("chromosome", "start.pos", "end.pos", "name", "score", "strand")

# Extracting A and B values from the 'name' column and calculating CNt
seq.dat$A <- as.numeric(sub(".*;A(\\d+);.*", "\\1", seq.dat$name))
seq.dat$B <- as.numeric(sub(".*;B(\\d+).*", "\\1", seq.dat$name))
seq.dat$CNt <- seq.dat$A + seq.dat$B

# Ensuring chromosome format is consistent for analysis
seq.dat$chromosome <- as.character(seq.dat$chromosome)
seq.dat$chromosome <- ifelse(seq.dat$chromosome == "23", "X", ifelse(seq.dat$chromosome == "24", "Y", seq.dat$chromosome))

ref <- tolower(args$build)
sample.id <- args$tumor

# Now seq.dat includes the required fields:
seq.dat <- preprocessHRD(seq.dat, ref) # Make sure preprocessHRD is ready to accept this structure

CN.dat <- getCNt(seq.dat)

# Computing HRD scores
HRD.score <- getHRD.Score(seq.dat, CN.dat, scaleTotal = FALSE)
HRD.LST <- getLST(seq.dat)
HRD.LOH <- getLOH(seq.dat)
HRD.NTAIr <- getNTAI.raw(seq.dat)
HRD.NTAIm <- getNTAI.norm(seq.dat, CN.dat)

# Preparing output
output <- data.frame(TumorID = args$tumor, HRD.NTAIr = HRD.NTAIr, HRD.NTAIm = HRD.NTAIm, HRD.LOH = HRD.LOH, HRD.LST = HRD.LST, HRD.score = HRD.score)
write.csv(output, file = args$outfile, row.names = FALSE, col.names = TRUE)
