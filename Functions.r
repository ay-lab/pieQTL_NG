#!/usr/bin/env Rscript

source('Header.r')

#===================================================
# custom logarithm function
# parameters: 
# mult: if -1, -log is computed
#===================================================
CustomLog <- function(InpVec, base=10, mult=1) {
	idx <- which(InpVec > 0)
	if (length(idx) > 0) {
		InpVec[idx] <- ((log(InpVec[idx]) / log(base)) * mult)	
	}
	newlist=list(vec=InpVec)
	return (newlist)
}

#=========================================
# function to get the number of fields for a particular file
GetNumFields <- function(inpfile) {
	nfield <- as.integer(system(paste("cat", inpfile, "| tail -n 1 | awk \'{print NF}\' -"), intern = TRUE))
	return(nfield)
}

#=============================================================
# function to extract loops for specific chromosome from an input interaction file
ExtractChrData <- function(InpFile, chrName, OutFile=NULL, header=TRUE, dist=c(-1,-1), mid=FALSE) {
	if (is.null(OutFile)) {
		OutFile <- paste0(dirname(InpFile), "/temp_Chr_data.bed")
	}

	# process the distance thresholds
	# and insert in two variables
	if ((dist[1] > 0) & (dist[2] > 0) & (dist[2] > dist[1])) {
		distthrlow <- dist[1]
		distthrhigh <- dist[2]		
	} else {
		distthrlow <- -1
		distthrhigh <- -1
	}

	# condition based on using gzipped input file
	# or plain text file

	if (file_ext(InpFile) == "gz") {
		if (header == TRUE) {
			if (mid == TRUE) {
				if ((distthrlow > 0) & (distthrhigh > 0)) {
					system(paste0("zcat ", InpFile, " | awk -v dt=", distthrlow, " -v dth=", distthrhigh, " \' function abs(v) {return v < 0 ? -v : v} {if (NR>1 && $1==\"", chrName, "\" && $3==\"", chrName, "\" && (abs($2-$4)>=dt) && (abs($2-$4)<=dth)) {print $0}}\' -  > ", OutFile))
				} else {
					system(paste0("zcat ", InpFile, " | awk \' {if (NR>1 && $1==\"", chrName, "\" && $3==\"", chrName, "\") {print $0}}\' -  > ", OutFile))
				}
			} else {
				if ((distthrlow > 0) & (distthrhigh > 0)) {
					system(paste0("zcat ", InpFile, " | awk -v dt=", distthrlow, " -v dth=", distthrhigh, " \' function abs(v) {return v < 0 ? -v : v} {if (NR>1 && $1==\"", chrName, "\" && $4==\"", chrName, "\" && (abs($2-$5)>=dt) && (abs($2-$5)<=dth)) {print $0}}\' -  > ", OutFile))
				} else {
					system(paste0("zcat ", InpFile, " | awk \' {if (NR>1 && $1==\"", chrName, "\" && $4==\"", chrName, "\") {print $0}}\' -  > ", OutFile))
				}
			}
		} else {
			if (mid == TRUE) {
				if ((distthrlow > 0) & (distthrhigh > 0)) {
					system(paste0("zcat ", InpFile, " | awk -v dt=", distthrlow, " -v dth=", distthrhigh, " \' function abs(v) {return v < 0 ? -v : v} {if ($1==\"", chrName, "\" && $3==\"", chrName, "\" && (abs($2-$4)>=dt) && (abs($2-$4)<=dth)) {print $0}}\' -  > ", OutFile))
				} else {
					system(paste0("zcat ", InpFile, " | awk \' {if ($1==\"", chrName, "\" && $3==\"", chrName, "\") {print $0}}\' -  > ", OutFile))						
				}
			} else {
				if ((distthrlow > 0) & (distthrhigh > 0)) {
					system(paste0("zcat ", InpFile, " | awk -v dt=", distthrlow, " -v dth=", distthrhigh, " \' function abs(v) {return v < 0 ? -v : v} {if ($1==\"", chrName, "\" && $4==\"", chrName, "\" && (abs($2-$5)>=dt) && (abs($2-$5)<=dth)) {print $0}}\' -  > ", OutFile))
				} else {
					system(paste0("zcat ", InpFile, " | awk \' {if ($1==\"", chrName, "\" && $4==\"", chrName, "\") {print $0}}\' -  > ", OutFile))
				}
			}
		}
	} else {
		if (header == TRUE) {
			if (mid == TRUE) {
				if ((distthrlow > 0) & (distthrhigh > 0)) {
					system(paste0("cat ", InpFile, " | awk -v dt=", distthrlow, " -v dth=", distthrhigh, " \' function abs(v) {return v < 0 ? -v : v} {if (NR>1 && $1==\"", chrName, "\" && $3==\"", chrName, "\" && (abs($2-$4)>=dt) && (abs($2-$4)<=dth)) {print $0}}\' -  > ", OutFile))
				} else {
					system(paste0("cat ", InpFile, " | awk \' {if (NR>1 && $1==\"", chrName, "\" && $3==\"", chrName, "\") {print $0}}\' -  > ", OutFile))
				}
			} else {
				if ((distthrlow > 0) & (distthrhigh > 0)) {
					system(paste0("cat ", InpFile, " | awk -v dt=", distthrlow, " -v dth=", distthrhigh, " \' function abs(v) {return v < 0 ? -v : v} {if (NR>1 && $1==\"", chrName, "\" && $4==\"", chrName, "\" && (abs($2-$5)>=dt) && (abs($2-$5)<=dth)) {print $0}}\' -  > ", OutFile))
				} else {
					system(paste0("cat ", InpFile, " | awk \' {if (NR>1 && $1==\"", chrName, "\" && $4==\"", chrName, "\") {print $0}}\' -  > ", OutFile))
				}
			}
		} else {
			if (mid == TRUE) {
				if ((distthrlow > 0) & (distthrhigh > 0)) {
					system(paste0("cat ", InpFile, " | awk -v dt=", distthrlow, " -v dth=", distthrhigh, " \' function abs(v) {return v < 0 ? -v : v} {if ($1==\"", chrName, "\" && $3==\"", chrName, "\" && (abs($2-$4)>=dt) && (abs($2-$4)<=dth)) {print $0}}\' -  > ", OutFile))
				} else {
					system(paste0("cat ", InpFile, " | awk \' {if ($1==\"", chrName, "\" && $3==\"", chrName, "\") {print $0}}\' -  > ", OutFile))						
				}
			} else {
				if ((distthrlow > 0) & (distthrhigh > 0)) {
					system(paste0("cat ", InpFile, " | awk -v dt=", distthrlow, " -v dth=", distthrhigh, " \' function abs(v) {return v < 0 ? -v : v} {if ($1==\"", chrName, "\" && $4==\"", chrName, "\" && (abs($2-$5)>=dt) && (abs($2-$5)<=dth)) {print $0}}\' -  > ", OutFile))
				} else {
					system(paste0("cat ", InpFile, " | awk \' {if ($1==\"", chrName, "\" && $4==\"", chrName, "\") {print $0}}\' -  > ", OutFile))
				}
			}
		}
	}

}	# end function

#=========================================
# function to get the number of lines for a particular file
GetNumLines <- function(inpfile) {
	nline <- as.integer(system(paste("cat", inpfile, "| wc -l"), intern = TRUE))
	return(nline)
}


#=============================================================
# function to compute and report the overlapping loops
OverlapLoop <- function(Inpdata1, Inpdata2, boundary=1, offset=0, uniqov=TRUE, IDX=FALSE) {

	# first compute the overlap between Inpdata1[,1:3],  Inpdata2[,1:3]
	# and between Inpdata1[,4:6],  Inpdata2[,4:6]
	ov1 <- as.data.frame(findOverlaps(GRanges(Inpdata1[,1], IRanges(Inpdata1[,2]+boundary-offset, Inpdata1[,3]-boundary+offset)),GRanges(Inpdata2[,1], IRanges(Inpdata2[,2]+boundary-offset, Inpdata2[,3]-boundary+offset))))
	ov2 <- as.data.frame(findOverlaps(GRanges(Inpdata1[,4], IRanges(Inpdata1[,5]+boundary-offset, Inpdata1[,6]-boundary+offset)),GRanges(Inpdata2[,4], IRanges(Inpdata2[,5]+boundary-offset, Inpdata2[,6]-boundary+offset))))
	overlap_uniq_mat <- ov1[unique(which(paste(ov1[,1], ov1[,2], sep=".") %in% paste(ov2[,1], ov2[,2], sep="."))),]

	# then compute the overlap between Inpdata1[,4:6],  Inpdata2[,1:3]
	# and between Inpdata1[,1:3],  Inpdata2[,4:6]
	# this is because input data can have first interval coordinate > second interval coordinate
	ov1A <- as.data.frame(findOverlaps(GRanges(Inpdata1[,4], IRanges(Inpdata1[,5]+boundary-offset, Inpdata1[,6]-boundary+offset)),GRanges(Inpdata2[,1], IRanges(Inpdata2[,2]+boundary-offset, Inpdata2[,3]-boundary+offset))))
	ov2A <- as.data.frame(findOverlaps(GRanges(Inpdata1[,1], IRanges(Inpdata1[,2]+boundary-offset, Inpdata1[,3]-boundary+offset)),GRanges(Inpdata2[,4], IRanges(Inpdata2[,5]+boundary-offset, Inpdata2[,6]-boundary+offset))))
	overlap_uniq_mat_1 <- ov1A[unique(which(paste(ov1A[,1], ov1A[,2], sep=".") %in% paste(ov2A[,1], ov2A[,2], sep="."))),]

	# now concatenate these two outcomes to get the final overlap statistics
	if (uniqov == TRUE) {
		ov_idx_file1 <- unique(c(overlap_uniq_mat[,1], overlap_uniq_mat_1[,1]))
		ov_idx_file2 <- unique(overlap_uniq_mat[,2], overlap_uniq_mat_1[,2])		
	} else {
		ov_idx_file1 <- c(overlap_uniq_mat[,1], overlap_uniq_mat_1[,1])
		ov_idx_file2 <- c(overlap_uniq_mat[,2], overlap_uniq_mat_1[,2])
	}
	nonov_idx_file1 <- setdiff(seq(1, nrow(Inpdata1)), ov_idx_file1)
	nonov_idx_file2 <- setdiff(seq(1, nrow(Inpdata2)), ov_idx_file2)

	# return the overlapping and non-overlapping set of indices
	if (IDX == FALSE) {
		newList <- list(A_AND_B = ov_idx_file1, B_AND_A = ov_idx_file2, A_MINUS_B = nonov_idx_file1, B_MINUS_A = nonov_idx_file2, A_AND_B.df = Inpdata1[ov_idx_file1, ], B_AND_A.df = Inpdata2[ov_idx_file2, ], A_MINUS_B.df = Inpdata1[nonov_idx_file1, ], B_MINUS_A.df = Inpdata2[nonov_idx_file2, ])
	} else {
		newList <- list(A_AND_B = ov_idx_file1, B_AND_A = ov_idx_file2, A_MINUS_B = nonov_idx_file1, B_MINUS_A = nonov_idx_file2)
	}
	return(newList)
}


#=============================================================
# function to compute and report the overlapping 1D segments 
Overlap1D <- function(Inpdata1, Inpdata2, boundary=1, offset=0, uniqov=TRUE, chrWise=FALSE) {

	if (chrWise == FALSE) {
		ov1 <- as.data.frame(findOverlaps(GRanges(Inpdata1[,1], IRanges(Inpdata1[,2]+boundary-offset, Inpdata1[,3]-boundary+offset)),GRanges(Inpdata2[,1], IRanges(Inpdata2[,2]+boundary-offset, Inpdata2[,3]-boundary+offset))))
		if (uniqov == TRUE) {
			ov_idx_file1 <- unique(ov1[,1])
			ov_idx_file2 <- unique(ov1[,2])		
		} else {
			ov_idx_file1 <- ov1[,1]
			ov_idx_file2 <- ov1[,2]
		}
		nonov_idx_file1 <- setdiff(seq(1, nrow(Inpdata1)), ov_idx_file1)
		nonov_idx_file2 <- setdiff(seq(1, nrow(Inpdata2)), ov_idx_file2)

	} else {
		# get the chromosome list from the given peaks
		ChrList <- union(unique(Inpdata1[,1]), unique(Inpdata2[,1]))

		ov_idx_file1 <- c()
		ov_idx_file2 <- c()
		nonov_idx_file1 <- c()
		nonov_idx_file2 <- c()
		
		# process individual chromosomes in a loop
		for (idx in (1:length(ChrList))) {
			CurrChr <- ChrList[idx]
			# peak indices of first data corresponding to the current chromosome
			Inpdata1_CurrChr_IdxSet <- which(Inpdata1[,1] == CurrChr)
			# peak indices of second data corresponding to the current chromosome
			Inpdata2_CurrChr_IdxSet <- which(Inpdata2[,1] == CurrChr)
			# proceed if both data have some entries for the current chromosome
			if ((length(Inpdata1_CurrChr_IdxSet) == 0) | (length(Inpdata2_CurrChr_IdxSet) == 0)) {
	 			next
			}
			# extract the input datasets limited to the current chromosome
			Inpdata1_CurrChr <- Inpdata1[Inpdata1_CurrChr_IdxSet, ]
			Inpdata2_CurrChr <- Inpdata2[Inpdata2_CurrChr_IdxSet, ]
			
			# compute the overlap
			ov1 <- as.data.frame(findOverlaps(GRanges(Inpdata1_CurrChr[,1], IRanges(Inpdata1_CurrChr[,2]+boundary-offset, Inpdata1_CurrChr[,3]-boundary+offset)),GRanges(Inpdata2_CurrChr[,1], IRanges(Inpdata2_CurrChr[,2]+boundary-offset, Inpdata2_CurrChr[,3]-boundary+offset))))
			
			# overlapping / non-overlapping indices for both categories
			if (uniqov == TRUE) {
				temp_ov_idx_file1 <- unique(ov1[,1])
				temp_ov_idx_file2 <- unique(ov1[,2])		
			} else {
				temp_ov_idx_file1 <- ov1[,1]
				temp_ov_idx_file2 <- ov1[,2]
			}
			temp_nonov_idx_file1 <- setdiff(seq(1, nrow(Inpdata1_CurrChr)), temp_ov_idx_file1)
			temp_nonov_idx_file2 <- setdiff(seq(1, nrow(Inpdata2_CurrChr)), temp_ov_idx_file2)

			# indices with respect to original loop indices
			# (i.e. with respect to Inpdata1 and Inpdata2)
			if (length(temp_ov_idx_file1) > 0) {
				ov_idx_file1_currchr <- Inpdata1_CurrChr_IdxSet[temp_ov_idx_file1]	
				# append in the final list
				ov_idx_file1 <- c(ov_idx_file1, ov_idx_file1_currchr)
			}
			
			if (length(temp_ov_idx_file2) > 0) {
				ov_idx_file2_currchr <- Inpdata2_CurrChr_IdxSet[temp_ov_idx_file2]
				# append in the final list
				ov_idx_file2 <- c(ov_idx_file2, ov_idx_file2_currchr)
			}
			
			if (length(temp_nonov_idx_file1) > 0) {
				nonov_idx_file1_currchr <- Inpdata1_CurrChr_IdxSet[temp_nonov_idx_file1]
				# append in the final list
				nonov_idx_file1 <- c(nonov_idx_file1, nonov_idx_file1_currchr)
			}
			
			if (length(temp_nonov_idx_file2) > 0) {
				nonov_idx_file2_currchr <- Inpdata2_CurrChr_IdxSet[temp_nonov_idx_file2]
				# append in the final list
				nonov_idx_file2 <- c(nonov_idx_file2, nonov_idx_file2_currchr)
			}
		}	# end chromosomes loop
	}

	# return the overlapping and non-overlapping set of indices
	newList <- list(A_AND_B = ov_idx_file1, B_AND_A = ov_idx_file2, A_MINUS_B = nonov_idx_file1, B_MINUS_A = nonov_idx_file2)
	return(newList)

}

#=========================================
# function to read contacts from 3C data
ReadContactData <- function(inpfile, headerInp=TRUE, chrNUM=FALSE, chrSET=c(1,4), mid=FALSE, midCOL=c(2), binSize=5000) {
    
    # handle gzipped file
	if (file_ext(inpfile) == "gz") {
    	outDF <- read.table(gzfile(inpfile), header=headerInp, sep="\t", stringsAsFactors=F)
    } else {
        outDF <- read.table(inpfile, header=headerInp, sep="\t", stringsAsFactors=F)
    }

    # handle chromosome intervals midpoints 
    if ((mid==TRUE) & (length(midCOL) > 0) & (binSize > 0)) {
        n <- ncol(outDF)
        for (i in (1:n)) {
            if (i %in% midCOL) {
                s1 <- outDF[,i] - (binSize/2)
                e1 <- outDF[,i] + (binSize/2)
                if (i == 1) {
                    ModDF <- cbind.data.frame(s1, e1)
                } else {
                    ModDF <- cbind.data.frame(ModDF, s1, e1)
                }
            } else {
                if (i == 1) {
                    ModDF <- outDF[,i]
                } else {
                    ModDF <- cbind.data.frame(ModDF, outDF[,i])
                }
            }
        }
        # reassign the modified data frame
        outDF <- ModDF
    }

    # handle chromosome numbers instead of name
    if (chrNUM==TRUE) {
        for (i in (1:length(chrSET))) {
            colno <- chrSET[i]      
            outDF[,colno] <- paste0('chr', as.character(outDF[,colno])) 
        }
    }
        
    return (outDF)
}
