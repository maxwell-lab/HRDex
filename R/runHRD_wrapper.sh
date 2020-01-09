#!/bin/bash

# get the number of input parameters
NPARAM=$#

# if user enters less than the required number of arguments, print the usage
if [ $NPARAM -lt 2 ]
then
	echo ""
	echo "USAGE :: "
	echo "./runHRD_wrapper.sh SEQFILENAME REF ... "
	echo "SEQFILENAME is the file containing seq data"
	echo "REF is the reference genome, grch37 or grch38"
	echo ""
	exit
fi

Rscript runHRD.R $1 $2
