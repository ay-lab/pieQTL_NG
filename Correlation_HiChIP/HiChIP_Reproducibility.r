#!/usr/bin/env Rscript

source('Header.r')

#============================
# schematic script to compute correlation between a pair of HiChIP contact maps
# as provided in the FitHiChIP output files "*.interactions_FitHiC.bed"

# computed for a single chromosome, say 'chr1'
#============================

currChr <- 'chr1'

##==== FitHiChIP file 1 (containing significance of all contacts) - the name is "*.interactions_FitHiC.bed"
allloop_inpfile1 <- 'File1_FitHiChIP.interactions_FitHiC.bed'

##==== FitHiChIP file 2 (containing significance of all contacts) - the name is "*.interactions_FitHiC.bed"
allloop_inpfile2 <- 'File2_FitHiChIP.interactions_FitHiC.bed'
			
# extract current chromosome specific interactions for the first file
currchr_allloop_inpfile1 <- paste0('allloop_inpfile1_', currChr, '.bed')
system(paste0("awk \'{if ((NR>1) && ($1==\"", currChr, "\")) {print $0}}\' ", allloop_inpfile1, " | cut -f1-7 > ", currchr_allloop_inpfile1))

# extract current chromosome specific interactions for the second file
currchr_allloop_inpfile2 <- paste0('allloop_inpfile2_', currChr, '.bed')
system(paste0("awk \'{if ((NR>1) && ($1==\"", currChr, "\")) {print $0}}\' ", allloop_inpfile2, " | cut -f1-7 > ", currchr_allloop_inpfile2))
			
# union of significant loops and their contact counts in both samples
commonloopfile <- paste0('Commonloop_', currChr, '.bed')
system(paste("bedtools pairtopair -a", currchr_allloop_inpfile1, "-b", currchr_allloop_inpfile2, "| cut -f1-7,14 > ", commonloopfile))

currdf <- data.table::fread(commonloopfile, header=F)
correlation_val <- cor(currdf[, (ncol(currdf) - 1)], currdf[, ncol(currdf)])
cat(sprintf("\n correlation of HiChIP contact maps : %s ", correlation_val))

