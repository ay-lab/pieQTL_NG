#!/usr/bin/env Rscript

source('Header.r')
source('Functions.r')

##============ input parameters
##============ user needs to edit these values
LoopFile <- 'FitHiChIP.interactions_FitHiC_Q0.01.bed'
TSSFile <- 'gEncode_Genes_Complete_NEW_TSS.gtf'
ChIPSeqFile <- 'macs2_peaks.narrowPeak'

##==== eQTL list file - assumed that the chromosome is mentioned in the first column and the SNP location is mentioned in the 2nd column
EQTLFile <- 'eQTL_list.bed'
EQTLChrCol <- 1
EQTLSNPCol <- 2

##=== offset (bp) between TSS and an interacting bin, to consider that as promoter bin
OffsetValue <- 5000

##=== offset (bp) between eQTL SNP location and an interacting bin, to consider that as eQTL bin
EQTLOffsetValue <- 5000
Outprefix <- 'pieQTL'
EQTLprefix <- 'Sample_EQTL'

##=== output directories - user can edit
outdir <- getwd()
eqtldir <- getwd()

##=== output files
AnnotatedLoopFile <- paste0(Outprefix, '_FitHiChIP_Annot_offset_', (OffsetValue / 1000), 'Kb.bed')
OutTextFile <- paste0(Outprefix, '_FitHiChIP_Annot_offset_', (OffsetValue / 1000), 'Kb.log')
OutPromoterLoopFile <- paste0(Outprefix, '_FitHiChIP_Prom_Loops_offset_', (OffsetValue / 1000), 'Kb.bed')
temp_ChrName_File <- paste0(Outprefix, '_temp_ChrName.bed')
temp_EQTLFile_CurrChr <- paste0(EQTLprefix, '_temp_EQTLFile_header.bed')

# P-Q direct connections (P: promoter, Q: EQTL)
Prom_Direct_EQTL_DumpFile <- paste0(EQTLprefix, '_Promoters_EQTL_Direct_PromSlack_', (OffsetValue / 1000), 'Kb_EQTLSlack_', (EQTLOffsetValue / 1000), 'Kb.bed')
Prom_Direct_EQTL_Dump_TSSFile <- paste0(EQTLprefix, '_Promoters_EQTL_Direct_PromSlack_', (OffsetValue / 1000), 'Kb_EQTLSlack_', (EQTLOffsetValue / 1000), 'Kb_TSS_Included.bed')

# P - P/E - Q indirect connections between them
Prom_Direct_1PorE_Indirect_EQTL_DumpFile <- paste0(EQTLprefix, '_Promoters_Direct_PorE_Indirect_EQTL_PromSlack_', (OffsetValue / 1000), 'Kb_EQTLSlack_', (EQTLOffsetValue / 1000), 'Kb.bed')
Prom_Direct_1PorE_Indirect_EQTL_Dump_TSSFile <- paste0(EQTLprefix, '_Promoters_Direct_PorE_Indirect_EQTL_PromSlack_', (OffsetValue / 1000), 'Kb_EQTLSlack_', (EQTLOffsetValue / 1000), 'Kb_TSS_Included.bed')

##==== stores direct pieQTLs + FitHiChIP loops
direct_pieQTL_loop_file <- 'Direct_pieQTL_FitHiChIP_Loops.bed')

##==== stores indirect pieQTLs + FitHiChIP loops
indirect_pieQTL_loop_file <- 'Indirect_pieQTL.bed'

##==== final list of pieQTLs (direct + indirect)
filtered_pieQTL_File <- 'dumped_pieQTL_FILT.txt'

##=== teporary output files
BinnedSeg_temp_File <- paste0(Outprefix, '_BinnedSeg_temp.bed')
BinnedSegFile <- paste0(Outprefix, '_BinnedSeg_Unique.bed')
AnnotLoopFile_CurrChr <- paste0(Outprefix, '_temp_AnnotLoopFile_CurrChr.bed')
EQTLFile_CurrChr <- paste0(Outprefix, '_temp_EQTLFile_CurrChr.bed')
GTFFile_CurrChr <- paste0(Outprefix, '_temp_GTFFile_CurrChr.bed')
GTFFile_GeneExprFile_CurrChr <- paste0(Outprefix, '_temp_GTFFile_GeneExprFile_CurrChr.bed')

#======================
# function to find the annotation of one interacting bin in FitHiChIP loops as promoters / enhancers / others
#======================
GetPromEnh_OneSegment <- function(Segdata, TSSdata, ChIPSeqData, Ov_Offset=5000) {

	OvTSS_Seg1 <- Overlap1D(Segdata, TSSdata, boundary=0, offset=Ov_Offset, uniqov=FALSE)
	PromIdx_Seg1 <- OvTSS_Seg1$A_AND_B
	OvPeak_Seg1 <- Overlap1D(Segdata, ChIPSeqData, boundary=0, offset=0, uniqov=FALSE)
	EnhIdx_Seg1 <- setdiff(OvPeak_Seg1$A_AND_B, PromIdx_Seg1)
	newList <- list(Prom = PromIdx_Seg1, Enh = EnhIdx_Seg1)
	return(newList)
}

#=========================
# dump loop annotation summary
#=========================
WriteLogAnnotLoop <- function(AnnotatedLoopData, OutTextFile, col1, col2) {

	P_P_Idx <- which((AnnotatedLoopData[,col1] == "P") & (AnnotatedLoopData[,col2] == "P"))
	P_E_Idx <- which(((AnnotatedLoopData[,col1] == "P") & (AnnotatedLoopData[,col2] == "E")) | ((AnnotatedLoopData[,col1] == "E") & (AnnotatedLoopData[,col2] == "P")))
	E_E_Idx <- which((AnnotatedLoopData[,col1] == "E") & (AnnotatedLoopData[,col2] == "E"))
	P_O_Idx <- which(((AnnotatedLoopData[,col1] == "P") & (AnnotatedLoopData[,col2] == "O")) | ((AnnotatedLoopData[,col1] == "O") & (AnnotatedLoopData[,col2] == "P")))
	E_O_Idx <- which(((AnnotatedLoopData[,col1] == "E") & (AnnotatedLoopData[,col2] == "O")) | ((AnnotatedLoopData[,col1] == "O") & (AnnotatedLoopData[,col2] == "E")))
	O_O_Idx <- which((AnnotatedLoopData[,col1] == "O") & (AnnotatedLoopData[,col2] == "O"))

	fp_out <- file(OutTextFile, "w")
	outtext <- paste0("\n Total number of loops: ", nrow(AnnotatedLoopData), "\n number of P-P loops: ", length(P_P_Idx), "\n number of P-E loops: ", length(P_E_Idx), "\n number of E-E loops: ", length(E_E_Idx), "\n number of P-O loops: ", length(P_O_Idx), "\n number of E-O loops: ", length(E_O_Idx), "\n number of O-O loops: ", length(O_O_Idx))
	writeLines(outtext, con=fp_out, sep="\n")
	close(fp_out)

}

#=========================
# function to annotate the interacting bins in FitHiChIP loops as promoters / enhancers / others
#=========================
Annotate_FitHiChIP_Loops_P_E <- function(InpLoopFile, InpTSSFile, InpChIPSeqFile, AnnotatedLoopFile, OutTextFile, Ov_Offset=5000) {

	FitHiChIPLoopData <- ReadContactData(InpLoopFile, headerInp=TRUE)
	InpGTFData <- data.table::fread(InpTSSFile, header=T)
	InpPeakData <- data.table::fread(InpChIPSeqFile, header=F)
	LabelVec_Seg1 <- rep("O", nrow(FitHiChIPLoopData))
	LabelVec_Seg2 <- rep("O", nrow(FitHiChIPLoopData))

	IdxSet1 <- GetPromEnh_OneSegment(FitHiChIPLoopData[,1:3], InpGTFData[,1:3], InpPeakData[,1:3], Ov_Offset)
	PromIdx_Seg1 <- IdxSet1$Prom
	EnhIdx_Seg1 <- IdxSet1$Enh
	LabelVec_Seg1[PromIdx_Seg1] <- "P"
	LabelVec_Seg1[EnhIdx_Seg1] <- "E"

	IdxSet2 <- GetPromEnh_OneSegment(FitHiChIPLoopData[,4:6], InpGTFData[,1:3], InpPeakData[,1:3], Ov_Offset)
	PromIdx_Seg2 <- IdxSet2$Prom
	EnhIdx_Seg2 <- IdxSet2$Enh
	LabelVec_Seg2[PromIdx_Seg2] <- "P"
	LabelVec_Seg2[EnhIdx_Seg2] <- "E"

	AnnotatedLoopData <- cbind.data.frame(FitHiChIPLoopData[,1:7], LabelVec_Seg1, LabelVec_Seg2, FitHiChIPLoopData[,8:ncol(FitHiChIPLoopData)])
	colnames(AnnotatedLoopData) <- c(colnames(FitHiChIPLoopData)[1:7], "Label1", "Label2", colnames(FitHiChIPLoopData)[8:ncol(FitHiChIPLoopData)])
	data.table::fwrite(AnnotatedLoopData, file = AnnotatedLoopFile, row.names = F, col.names = T, sep = "\t", quote=F, append=F)

	WriteLogAnnotLoop(AnnotatedLoopData, OutTextFile, 8, 9)

}

#===================================
# extract the promoter centric FitHiChIP loops
#===================================
Get_Promoter_Loops_FitHiChIP <- function(AnnotatedLoopFile, InpTSSFile, OutPromoterLoopFile, Ov_Offset=5000) {	

	AnnotatedLoopData <- UtilRPckg::ReadContactData(AnnotatedLoopFile, headerInp=TRUE)
	cat(sprintf("\n *** within function Get_Promoter_Loops_FitHiChIP -- read the input loops"))

	InpGTFData <- data.table::fread(InpTSSFile, header=T)
	cat(sprintf("\n *** within function Get_Promoter_Loops_FitHiChIP -- read the input GTF file"))

	ChrList <- sort(unique(as.vector(AnnotatedLoopData[,1])))

	# boolean flag for output file write
	outloop_flag <- FALSE

	for (chrIdx in 1:length(ChrList)) {
		currChr <- ChrList[chrIdx]
		cat(sprintf("\n\n ********* Processing current chromosome : %s ******** \n\n", currChr))

		AnnotatedLoopData_currChr <- AnnotatedLoopData[which(AnnotatedLoopData[,1] == currChr), ]
		if (nrow(AnnotatedLoopData_currChr) == 0) {
			cat(sprintf("\n No loop for the current chromosome - continue"))
			next
		}
		InpGTFData_currChr <- InpGTFData[which(InpGTFData[,1] == currChr), ]
		if (nrow(InpGTFData_currChr) == 0) {
			cat(sprintf("\n No TSS information for the current chromosome - continue"))
			next
		}

		#======================
		# now find the promoters and corresponding interacting segments
		#======================
		# find the overlap of AnnotatedLoopData_currChr (first interacting segment) with the gene (TSS) (slack = Ov_Offset)
		CurrOv1 <- Overlap1D(InpGTFData_currChr[,1:3], AnnotatedLoopData_currChr[,1:3], boundary=0, offset=Ov_Offset, uniqov=FALSE)
		cat(sprintf("\n -->> Computed overlap of GTF with the first interacting segment"))
		GTF_Idx_1 <- CurrOv1$A_AND_B
		Loop_Idx_1 <- CurrOv1$B_AND_A
		# cat(sprintf("\n -->> length GTF_Idx_1 : %s  length Loop_Idx_1 : %s ", length(GTF_Idx_1), length(Loop_Idx_1)))

		if ((length(GTF_Idx_1) > 0) & (length(Loop_Idx_1) > 0)) {	
			# construct a data frame by using the GTF gene information
			# and the second interacting segment of the overlapping loops
			# we also include the promoter bin
			MergedDataIdxVec <- c(1:7,9,16:21)
			OutDF1 <- cbind.data.frame(InpGTFData_currChr[GTF_Idx_1, c(1:6, 8)], AnnotatedLoopData_currChr[Loop_Idx_1, MergedDataIdxVec])
			OutDFNameVec <- c("TSSc", "TSSs", "TSSe", "TSSIdx", "Strand", "GeneID", "GeneName", "chrP", "startP", "endP", "chr", "start", "end", "cc", "label", "Coverage", "isPeak", "Bias", "Mapp", "GCContent", "RESites")		
			OutDF1 <- cbind.data.frame(OutDF1, AnnotatedLoopData_currChr[Loop_Idx_1, 22:27])
			OutDFNameVec <- c(OutDFNameVec, colnames(AnnotatedLoopData_currChr)[22:27])
			colnames(OutDF1) <- OutDFNameVec
		}

		# find the overlap of AnnotatedLoopData_currChr (second interacting segment) with the gene (TSS) (slack = Ov_Offset)
		CurrOv2 <- Overlap1D(InpGTFData_currChr[,1:3], AnnotatedLoopData_currChr[,4:6], boundary=0, offset=Ov_Offset, uniqov=FALSE)
		cat(sprintf("\n -->> Computed overlap of GTF with the second interacting segment"))
		GTF_Idx_2 <- CurrOv2$A_AND_B
		Loop_Idx_2 <- CurrOv2$B_AND_A
		# cat(sprintf("\n -->> length GTF_Idx_2 : %s  length Loop_Idx_2 : %s ", length(GTF_Idx_2), length(Loop_Idx_2)))

		if ((length(GTF_Idx_2) > 0) & (length(Loop_Idx_2) > 0)) {
			# construct a data frame by using the GTF gene information
			# and the second interacting segment of the overlapping loops
			# we also include the promoter bin
			MergedDataIdxVec <- c(4:6,1:3,7,8,10:15)
			OutDF2 <- cbind.data.frame(InpGTFData_currChr[GTF_Idx_2, c(1:6, 8)], AnnotatedLoopData_currChr[Loop_Idx_2, MergedDataIdxVec])
			OutDFNameVec <- c("TSSc", "TSSs", "TSSe", "TSSIdx", "Strand", "GeneID", "GeneName", "chrP", "startP", "endP", "chr", "start", "end", "cc", "label", "Coverage", "isPeak", "Bias", "Mapp", "GCContent", "RESites")
			OutDF2 <- cbind.data.frame(OutDF2, AnnotatedLoopData_currChr[Loop_Idx_2, 22:27])
			OutDFNameVec <- c(OutDFNameVec, colnames(AnnotatedLoopData_currChr)[22:27])
			colnames(OutDF2) <- OutDFNameVec
		}	

		# combine these two data frames
		if ((length(GTF_Idx_1) > 0) & (length(Loop_Idx_1) > 0) & (length(GTF_Idx_2) > 0) & (length(Loop_Idx_2) > 0)) {
			MasterSheetDF <- rbind.data.frame(OutDF1, OutDF2)	
		} else if ((length(GTF_Idx_1) > 0) & (length(Loop_Idx_1) > 0)) {
			MasterSheetDF <- OutDF1
		} else if ((length(GTF_Idx_2) > 0) & (length(Loop_Idx_2) > 0)) {
			MasterSheetDF <- OutDF2
		}

		# sort this master data frame with respect to the genes
		if (((length(GTF_Idx_1) > 0) & (length(Loop_Idx_1) > 0)) | ((length(GTF_Idx_2) > 0) & (length(Loop_Idx_2) > 0))) {
			# sort the data frame with respect to 
			# columns 2 (number), 12 (number) - chromosome intervals
			MasterSheetDF <- MasterSheetDF[ order( MasterSheetDF[,2], MasterSheetDF[,12] ), ]

			# write promoter centric interactions
			if (outloop_flag == FALSE) {
				write.table(MasterSheetDF, OutPromoterLoopFile, row.names=F, col.names=T, sep = "\t", quote=F, append=F)
				outloop_flag <- TRUE
			} else {
				write.table(MasterSheetDF, OutPromoterLoopFile, row.names=F, col.names=F, sep = "\t", quote=F, append=T)
			}
		}	# end condition 

	}	# end chromosome loop
}	

#==================================
# get the promoters which directly connect to one P/E and indirectly connect to one EQTL
# the P/E and EQTL segments are connected
#==================================
Get_Prom_Direct_1PorE_Indirect_EQTL <- function(Label_BinnedSeg_CurrChr, Conn_Idx_DF_Dist1, Conn_Idx_DF_Dist2, Prom_Nodes_Idx, Enh_Nodes_Idx, EQTL_Nodes_Idx) {

	# boolean variable 
	data_out <- FALSE

	# convert to data frames before processing
	Conn_Idx_DF_Dist1 <- as.data.frame(Conn_Idx_DF_Dist1)
	colnames(Conn_Idx_DF_Dist1) <- c("idx1", "idx2")
	Conn_Idx_DF_Dist2 <- as.data.frame(Conn_Idx_DF_Dist2)
	colnames(Conn_Idx_DF_Dist2) <- c("idx1", "idx2")

	# From the structure Conn_Idx_DF_Dist1, get the interacting segments 
	Dist1_Node1 <- Conn_Idx_DF_Dist1[, 1]
	Dist1_Node2 <- Conn_Idx_DF_Dist1[, 2]

	# From the structure Conn_Idx_DF_Dist2, get the interacting segments 
	Dist2_Node1 <- Conn_Idx_DF_Dist2[, 1]
	Dist2_Node2 <- Conn_Idx_DF_Dist2[, 2]

	cat(sprintf("\n\n\n *** Within function Get_Prom_Direct_1PorE_Indirect_EQTL *** \n\n\n"))

	# get the connections between P/E to EQTL segments (shortest distance = 1)
	# first column: denote by PorE
	# second column: denote by EQTL
	PorE_EQTL_Conn_DF <- Conn_Idx_DF_Dist1[which((Dist1_Node1 %in% union(Prom_Nodes_Idx, Enh_Nodes_Idx)) & (Dist1_Node2 %in% EQTL_Nodes_Idx)), ]
	colnames(PorE_EQTL_Conn_DF) <- c("PorE", "EQTL")
	cat(sprintf("\n nrow PorE_EQTL_Conn_DF : %s ", nrow(PorE_EQTL_Conn_DF)))

	# get prom_PorE connections
	# first column: denote by Prom
	# second column: denote by PorE
	prom_PorE_DF <- Conn_Idx_DF_Dist1[which((Dist1_Node1 %in% Prom_Nodes_Idx) & (Dist1_Node2 %in% union(Prom_Nodes_Idx, Enh_Nodes_Idx))), ]
	colnames(prom_PorE_DF) <- c("Prom", "PorE")
	cat(sprintf("\n nrow prom_PorE_DF : %s ", nrow(prom_PorE_DF)))

	# promoters indirectly connected with EQTL (shortest distance = 2)
	# first column: denote by Prom
	# second column: denote by EQTL
	prom_EQTL_Indirect_DF <- Conn_Idx_DF_Dist2[which((Dist2_Node1 %in% Prom_Nodes_Idx) & (Dist2_Node2 %in% EQTL_Nodes_Idx)), ]
	colnames(prom_EQTL_Indirect_DF) <- c("Prom", "EQTL")
	cat(sprintf("\n nrow prom_EQTL_Indirect_DF : %s ", nrow(prom_EQTL_Indirect_DF)))

	# nested inner join
	finalPromEnhDF <- (dplyr::inner_join(prom_PorE_DF, PorE_EQTL_Conn_DF) %>% dplyr::inner_join(prom_EQTL_Indirect_DF))

	if (nrow(finalPromEnhDF) == 0) {
		newList <- list(data_out = FALSE, finalDF = NULL)
	} else {
		prom_idx <- finalPromEnhDF$Prom 	#finalPromEnhDF[, 1]
		p_or_e_idx <- finalPromEnhDF$PorE 	#finalPromEnhDF[, 2]
		eqtl_idx <- finalPromEnhDF$EQTL 	#finalPromEnhDF[, 3]
		finalDF <- cbind.data.frame(Label_BinnedSeg_CurrChr[prom_idx, 1:3], Label_BinnedSeg_CurrChr[p_or_e_idx, 1:3], Label_BinnedSeg_CurrChr[eqtl_idx, 1:3])
		colnames(finalDF) <- c("PChr", "PStart", "PEnd", "Dist1PorEChr", "Dist1PorEStart", "Dist1PorEEnd", "Dist2EQTLChr", "Dist2EQTLBinStart", "Dist2EQTLBinEnd")
		newList <- list(data_out = TRUE, finalDF = finalDF)
	}

	return(newList)

}	# end function

#==============================
# this function returns the promoters 
# which are directly connected to one EQTL segment
#==============================
Get_Prom_Direct_EQTL <- function(Label_BinnedSeg_CurrChr, Conn_Idx_DF_Dist1, Prom_Nodes_Idx, EQTL_Nodes_Idx) {

	# boolean variable 
	data_out <- FALSE

	# convert to data frames before processing
	Conn_Idx_DF_Dist1 <- as.data.frame(Conn_Idx_DF_Dist1)
	colnames(Conn_Idx_DF_Dist1) <- c("idx1", "idx2")	

	cat(sprintf("\n\n\n *** Within function Get_Prom_Direct_EQTL *** \n\n\n"))

	# From the structure Conn_Idx_DF_Dist1, get the interacting segments 
	Dist1_Node1 <- Conn_Idx_DF_Dist1[, 1]
	Dist1_Node2 <- Conn_Idx_DF_Dist1[, 2]

	# promoters connected with EQTL
	# first column: denote by Prom
	# second column: denote by EQTL
	prom_EQTL_DF <- Conn_Idx_DF_Dist1[which((Dist1_Node1 %in% Prom_Nodes_Idx) & (Dist1_Node2 %in% EQTL_Nodes_Idx)), ]
	colnames(prom_EQTL_DF) <- c("Prom", "EQTL")
	cat(sprintf("\n nrow prom_EQTL_DF : %s ", nrow(prom_EQTL_DF)))	

	if (nrow(prom_EQTL_DF) == 0) {
		newList <- list(data_out = FALSE, finalDF = NULL)
	} else {
		prom_idx <- prom_EQTL_DF$Prom	
		eqtl_idx <- prom_EQTL_DF$EQTL
		finalDF <- cbind.data.frame(Label_BinnedSeg_CurrChr[prom_idx, 1:3], Label_BinnedSeg_CurrChr[eqtl_idx, 1:3])
		colnames(finalDF) <- c("PChr", "PStart", "PEnd", "EQTLChr", "EQTLBinStart", "EQTLBinEnd")
		newList <- list(data_out = TRUE, finalDF = finalDF)
	}
	
	return(newList)

}	# end function

##==================
# main code
##==================

AnnotatedLoopFile <- paste0(Outprefix, '_FitHiChIP_Annot_offset_', (OffsetValue / 1000), 'Kb.bed')
OutTextFile <- paste0(Outprefix, '_FitHiChIP_Annot_offset_', (OffsetValue / 1000), 'Kb.log')
OutPromoterLoopFile <- paste0(Outprefix, '_FitHiChIP_Prom_Loops_offset_', (OffsetValue / 1000), 'Kb.bed')

## annotate P-E loops
Annotate_FitHiChIP_Loops_P_E(LoopFile, TSSFile, ChIPSeqFile, AnnotatedLoopFile, OutTextFile, OffsetValue)	

## get promoter specific FitHiChIP loops
Get_Promoter_Loops_FitHiChIP(AnnotatedLoopFile, TSSFile, OutPromoterLoopFile, OffsetValue)

# get the bin size of participating loops
BinSize <- as.integer(system(paste0("awk \'{if (NR==2) {print ($3-$2)}}\' ", LoopFile), intern = TRUE))

# get the chromosomes participating in the interactions
temp_ChrName_File <- paste0(Outprefix, '_temp_ChrName.bed')
system(paste0("awk \'{if (NR>1) {print $1}}\' ", LoopFile, " | sort -k1,1 | uniq > ", temp_ChrName_File))
ChrNamesDF <- read.table(temp_ChrName_File, header=F)
chrnames <- ChrNamesDF[, 1]

# remove the temporary file
# system(paste("rm", temp_ChrName_File))

# get the EQTL header information
temp_EQTLFile_CurrChr <- paste0(EQTLprefix, '_temp_EQTLFile_header.bed')
system(paste("cat", EQTLFile, "| head -n 2 >", temp_EQTLFile_CurrChr))	
eqtl_header_data <- read.table(temp_EQTLFile_CurrChr, header=T, sep="\t", stringsAsFactors=F)
eqtl_header_info <- colnames(eqtl_header_data)
system(paste("rm", temp_EQTLFile_CurrChr))

#==========================
# P-Q direct connections (P: promoter, Q: EQTL)
Prom_Direct_EQTL_DumpFile <- paste0(EQTLprefix, '_Promoters_EQTL_Direct_PromSlack_', (OffsetValue / 1000), 'Kb_EQTLSlack_', (EQTLOffsetValue / 1000), 'Kb.bed')
Prom_Direct_EQTL_Dump_TSSFile <- paste0(EQTLprefix, '_Promoters_EQTL_Direct_PromSlack_', (OffsetValue / 1000), 'Kb_EQTLSlack_', (EQTLOffsetValue / 1000), 'Kb_TSS_Included.bed')

# boolean variable regarding this connection (chromosome wise processing)
bool_DF_Prom_Direct_EQTL <- FALSE

#==========================
# P - P/E - Q indirect connections between them
# here promoter and enhancer bins are defined according to the specified slack (OffsetValue)
# EQTL bins are also defined according to the specified slack (EQTLOffsetValue)
Prom_Direct_1PorE_Indirect_EQTL_DumpFile <- paste0(EQTLprefix, '_Promoters_Direct_PorE_Indirect_EQTL_PromSlack_', (OffsetValue / 1000), 'Kb_EQTLSlack_', (EQTLOffsetValue / 1000), 'Kb.bed')
Prom_Direct_1PorE_Indirect_EQTL_Dump_TSSFile <- paste0(EQTLprefix, '_Promoters_Direct_PorE_Indirect_EQTL_PromSlack_', (OffsetValue / 1000), 'Kb_EQTLSlack_', (EQTLOffsetValue / 1000), 'Kb_TSS_Included.bed')

# boolean variable regarding this connection (chromosome wise processing)
bool_DF_Prom_Direct_1PorE_Indirect_EQTL <- FALSE

#========================
# a few temporary files
# which can be deleted later
BinnedSeg_temp_File <- paste0(Outprefix, '_BinnedSeg_temp.bed')
BinnedSegFile <- paste0(Outprefix, '_BinnedSeg_Unique.bed')
AnnotLoopFile_CurrChr <- paste0(Outprefix, '_temp_AnnotLoopFile_CurrChr.bed')
EQTLFile_CurrChr <- paste0(Outprefix, '_temp_EQTLFile_CurrChr.bed')
GTFFile_CurrChr <- paste0(Outprefix, '_temp_GTFFile_CurrChr.bed')
GTFFile_GeneExprFile_CurrChr <- paste0(Outprefix, '_temp_GTFFile_GeneExprFile_CurrChr.bed')
#========================

#================================
# chromosome specific loop to process individual CIS interactions
#================================
for (i in (1:length(chrnames))) {

	cat(sprintf("\n\n\n\n ************** processing current chromosome : %s ", chrnames[i]))

	# extract the TSS information for the current chromosome
	# subject to the amount of slack mentioned in the input arguments
	system(paste0("awk \'((NR==1) || ($1==\"", chrnames[i], "\"))\' ", TSSFile, " > ", GTFFile_CurrChr))
	n <- GetNumLines(GTFFile_CurrChr)
	if (n <= 1) {
		next
	}	
	cat(sprintf("\n ====>>>> number of rows of the extracted GTF file for the current chromosome : %s ", n))

	# similarly extract the loops for the current chromosome
	# also include the header information
	system(paste0("awk \'((NR==1) || ((NR>1) && ($1==\"", chrnames[i], "\")))\' ", AnnotatedLoopFile, " > ", AnnotLoopFile_CurrChr))
	n <- GetNumLines(AnnotLoopFile_CurrChr)
	if (n <= 1) {
		next
	}	
	AnnotLoopData_CurrChr <- data.table::fread(AnnotLoopFile_CurrChr, header=T)	
	cat(sprintf("\n number of rows of AnnotLoopData_CurrChr (unique) : %s ", nrow(AnnotLoopData_CurrChr)))

	# process the annotated loops for the current chromosome
	# obtain the interacting bins (sorted) along with their annotations 
	# and assemble them in a data structure	
	system(paste0("awk \'{if ((NR>1) && ($1==\"", chrnames[i], "\")) {print $1\"\t\"$2\"\t\"$3\"\t\"$8}}\' ", AnnotLoopFile_CurrChr, " > ", BinnedSeg_temp_File))
	system(paste0("awk \'{if ((NR>1) && ($4==\"", chrnames[i], "\")) {print $4\"\t\"$5\"\t\"$6\"\t\"$9}}\' ", AnnotLoopFile_CurrChr, " >> ", BinnedSeg_temp_File))
	system(paste0("cat ", BinnedSeg_temp_File, " | sort -k1,1 -k2,2n | uniq > ", BinnedSegFile))

	# this structure stores the interacting bins
	Label_BinnedSeg_CurrChr <- data.table::fread(BinnedSegFile, header=F)	
	colnames(Label_BinnedSeg_CurrChr) <- c("chr", "start", "end", "label")
	cat(sprintf("\n number of rows of Label_BinnedSeg_CurrChr (unique) : %s ", nrow(Label_BinnedSeg_CurrChr)))
	if (nrow(Label_BinnedSeg_CurrChr) == 0) {
		next
	}	
	
	#================================================
	# list the enhancer and promoters for the current chromosome
	#================================================
	Enh_Nodes_Idx <- which(Label_BinnedSeg_CurrChr[,4] == "E")
	Prom_Nodes_Idx <- which(Label_BinnedSeg_CurrChr[,4] == "P")
	cat(sprintf("\n number of enhancer nodes for this chromosome : %s  \n number of promoter nodes for this chromosome : %s  subject to a slack of  %s ", length(Enh_Nodes_Idx), length(Prom_Nodes_Idx), OffsetValue))

	#================================================
	# ov1$B_AND_A provides the indices of the graph nodes which is interacting (first segment)
	ov1 <- Overlap1D(AnnotLoopData_CurrChr[, 1:3], Label_BinnedSeg_CurrChr, boundary=1, offset=0, uniqov=FALSE)

	# basically, ov2$B_AND_A provides the indices of the graph nodes which is interacting (second segment)
	ov2 <- Overlap1D(AnnotLoopData_CurrChr[, 4:6], Label_BinnedSeg_CurrChr, boundary=1, offset=0, uniqov=FALSE)

	# create a data frame which serves as the names of nodes
	# index of the "Label_BinnedSeg_CurrChr" is the name for nodes
	# it also stores the type of segments (i.e. promoters or enhancers)
	NodeNameDF <- data.frame(name=seq(1,nrow(Label_BinnedSeg_CurrChr)), type=Label_BinnedSeg_CurrChr[,4])

	# form the edges in a data frame 
	# connections are modeled with respect to the node indices of Label_BinnedSeg
	# we also store the contact count of respective loops
	EdgesDF <- data.frame(from=ov1$B_AND_A, to=ov2$B_AND_A, cc=AnnotLoopData_CurrChr[, 7], qval=AnnotLoopData_CurrChr[, ncol(AnnotLoopData_CurrChr)], dist=abs(AnnotLoopData_CurrChr[, 5] - AnnotLoopData_CurrChr[, 2]))
	if (nrow(EdgesDF) == 0) {
		next
	}

	# form the graph using the vertices and edges data frame 
	InpGraph <- graph_from_data_frame(EdgesDF, directed=FALSE, vertices=NodeNameDF)

	# compute the shortest path distances for all pairs of vertices
	# returns a matrix with N X N dimension
	ShortDist <- distances(InpGraph)

	# replace the "inf" entries in the shortest distance matrix as (-1)
	# check https://stackoverflow.com/questions/30990961/replace-inf-nan-and-na-values-with-zero-in-a-dataset-in-r
	is.na(ShortDist) <- sapply(ShortDist, is.infinite)
	ShortDist[is.na(ShortDist)] <- -1

	# first get the indices of connections (i.e. shortest distance = 1)
	# this function returns a data frame with two columns: row_idx and col_idx
	# rows indicate samples
	Conn_Idx_DF_Dist1 <- which( ShortDist == 1, arr.ind=T )
	cat(sprintf("\n short distance = 1 : nrow Conn_Idx_DF_Dist1 : %s ", nrow(Conn_Idx_DF_Dist1)))

	# similarly get the indices of shortest distance = 2
	# between promoters and / or enhancers
	Conn_Idx_DF_Dist2 <- which( ShortDist == 2, arr.ind=T )
	cat(sprintf("\n short distance = 2 : nrow Conn_Idx_DF_Dist2 : %s ", nrow(Conn_Idx_DF_Dist2)))

	# similarly get the indices of shortest distance = 2
	# between promoters and / or enhancers
	Conn_Idx_DF_Dist3 <- which( ShortDist == 3, arr.ind=T )
	cat(sprintf("\n short distance = 3 : nrow Conn_Idx_DF_Dist3 : %s ", nrow(Conn_Idx_DF_Dist3)))

	# extract the EQTL's for the current chromosome
	# EQTLChrCol: contains chromosome name
	# also include EQTLSNPCol (SNP information)
	# don't include the header
	system(paste0("awk \'{if ((NR>1) && ($", EQTLChrCol, "==\"", chrnames[i], "\")) {print $", EQTLChrCol, "\"\t\"$", EQTLSNPCol, "\"\t\"$", EQTLSNPCol, "\"\t\"$0}}\' ", EQTLFile, " > ", EQTLFile_CurrChr))
	n <- GetNumLines(EQTLFile_CurrChr)
	if (n == 0) {
		next
	}
	EQTLData_CurrChr <- data.table::fread(EQTLFile_CurrChr, header=F)

	# extract only the chromosome and SNP; include column names
	OnlyEQTLInfo_CurrChr <- EQTLData_CurrChr[, 1:3]
	colnames(OnlyEQTLInfo_CurrChr) <- c("EQTLChr", "EQTLSNP", "EQTLSNP2")

	# compute the EQTL containing bin (slack = EQTLOffsetValue - provided as input parameter)
	ov <- Overlap1D(OnlyEQTLInfo_CurrChr[, 1:3], Label_BinnedSeg_CurrChr, boundary=0, offset=EQTLOffsetValue, uniqov=F)

	if (length(ov$A_AND_B) > 0) {			
		# EQTL containing nodes (indices are numbered with respect to Label_BinnedSeg_CurrChr)
		EQTL_Nodes_Idx <- unique(ov$B_AND_A)
		cat(sprintf("\n number of EQTL nodes for this chromosome : %s  (subject to a slack of  %s  between EQTL and a bin) ", length(EQTL_Nodes_Idx), EQTLOffsetValue))
	} else {
		EQTL_Nodes_Idx <- c()
	}
	
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# get the P-EQTL direct connections
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~		
			
	if (length(EQTL_Nodes_Idx) > 0) {
		OutRes <- Get_Prom_Direct_EQTL(Label_BinnedSeg_CurrChr, Conn_Idx_DF_Dist1, Prom_Nodes_Idx, EQTL_Nodes_Idx)
		if (OutRes$data_out == TRUE) {
			OutDF <- OutRes$finalDF
			CN <- colnames(OutDF)

			# first compute the overlap between EQTL bins and EQTL locations
			# apply the EQTL slack also
			# place EQTL bins at the beginning
			tempfile_1 <- paste0(dirname(Prom_Direct_EQTL_Dump_TSSFile), '/tempfile_1.bed')
			write.table(cbind.data.frame(OutDF[, 4], (OutDF[, 5] - EQTLOffsetValue), (OutDF[, 6] + EQTLOffsetValue), OutDF[, 1:3]), tempfile_1, row.names=F, col.names=F, sep="\t", quote=F, append=F)

			# now compute the overlap with reference EQTL information file
			# use bedtools intersect
			# in the EQTL data, skip the first three fields
			# also reorder the data such that promoter bins are now placed at front
			# extend the promoter bin intervals to apply the promoter slack				
			tempfile_2 <- paste0(dirname(Prom_Direct_EQTL_Dump_TSSFile), '/tempfile_2.bed')
			system(paste0("bedtools intersect -a ", tempfile_1, " -b ", EQTLFile_CurrChr, " -wa -wb | awk -v FS=[\'\t\'] -v t=", EQTLOffsetValue, " -v p=", OffsetValue, " \'{if (NR>1) {printf \"\\n\"}; printf $4\"\t\"($5-p)\"\t\"($6+p)\"\t\"; printf $1\"\t\"($2+t)\"\t\"($3-t)\"\t\"; for (i=10;i<NF;i++){printf $i\"\t\"}; printf $NF}\' - > ", tempfile_2))

			# now apply the overlap with respect to TSS sites
			# using the promoter slack
			tempfile_3 <- paste0(dirname(Prom_Direct_EQTL_Dump_TSSFile), '/tempfile_3.bed')
			system(paste0("bedtools intersect -a ", GTFFile_CurrChr, " -b ", tempfile_2, " -wa -wb | awk -v FS=[\'\t\'] -v t=", OffsetValue, " \'{if (NR>1) {printf \"\\n\"}; for (i=1;i<=8;i++){printf $i\"\t\"}; printf $9\"\t\"($10+t)\"\t\"($11-t)\"\t\"; for (i=12;i<=NF;i++){printf $i \"\t\"}}\' - > ", tempfile_3))

			if (bool_DF_Prom_Direct_EQTL == FALSE) {
				# write the interacting bins
				write.table(OutDF, Prom_Direct_EQTL_DumpFile, row.names=F, col.names=T, sep="\t", quote=F, append=F)
				# write the header line of TSS information file
				fp_out <- file(Prom_Direct_EQTL_Dump_TSSFile, "w")
				writeLines(paste(c("TSSc", "TSSs", "TSSe", "TSSIdx", "Strand", "GeneID", "GeneType", "GeneName", CN, eqtl_header_info), collapse="\t"), con=fp_out)
				close(fp_out)
				# set the boolean flag
				bool_DF_Prom_Direct_EQTL <- TRUE					
			} else {
				write.table(OutDF, Prom_Direct_EQTL_DumpFile, row.names=F, col.names=F, sep="\t", quote=F, append=T)
			}
			# append the TSS file information
			system(paste("cat", tempfile_3, ">>", Prom_Direct_EQTL_Dump_TSSFile))
			# remove temporary files
			system(paste("rm", tempfile_1))
			system(paste("rm", tempfile_2))	
			system(paste("rm", tempfile_3))
		}				
	}

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~			
	# get the P-P/E (direct) and P---EQTL (indirect) connections
	# where P/E-EQTL are also connected
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~		

	if (length(EQTL_Nodes_Idx) > 0) {
		OutRes <- Get_Prom_Direct_1PorE_Indirect_EQTL(Label_BinnedSeg_CurrChr, Conn_Idx_DF_Dist1, Conn_Idx_DF_Dist2, Prom_Nodes_Idx, Enh_Nodes_Idx, EQTL_Nodes_Idx)

		if (OutRes$data_out == TRUE) {
			OutDF <- OutRes$finalDF
			CN <- colnames(OutDF)

			# first compute the overlap between EQTL bins and EQTL locations
			# apply the EQTL slack also
			# place EQTL bins at the beginning
			tempfile_1 <- paste0(dirname(Prom_Direct_1PorE_Indirect_EQTL_Dump_TSSFile), '/tempfile_1.bed')
			write.table(cbind.data.frame(OutDF[, 7], (OutDF[, 8] - EQTLOffsetValue), (OutDF[, 9] + EQTLOffsetValue), OutDF[, 1:6]), tempfile_1, row.names=F, col.names=F, sep="\t", quote=F, append=F)

			# now compute the overlap with reference EQTL information file
			# use bedtools intersect
			# in the EQTL data, skip the first three fields
			# also reorder the data such that promoter bins are now placed at front
			# extend the promoter bin intervals to apply the promoter slack				
			tempfile_2 <- paste0(dirname(Prom_Direct_1PorE_Indirect_EQTL_Dump_TSSFile), '/tempfile_2.bed')
			system(paste0("bedtools intersect -a ", tempfile_1, " -b ", EQTLFile_CurrChr, " -wa -wb | awk -v FS=[\'\t\'] -v t=", EQTLOffsetValue, " -v p=", OffsetValue, " \'{if (NR>1) {printf \"\\n\"}; printf $4\"\t\"($5-p)\"\t\"($6+p)\"\t\"; for (i=7;i<=9;i++){printf $i\"\t\"}; printf $1\"\t\"($2+t)\"\t\"($3-t)\"\t\"; for (i=13;i<NF;i++){printf $i\"\t\"}; printf $NF}\' - > ", tempfile_2))

			# now apply the overlap with respect to TSS sites
			# using the promoter slack
			tempfile_3 <- paste0(dirname(Prom_Direct_1PorE_Indirect_EQTL_Dump_TSSFile), '/tempfile_3.bed')
			system(paste0("bedtools intersect -a ", GTFFile_CurrChr, " -b ", tempfile_2, " -wa -wb | awk -v FS=[\'\t\'] -v t=", OffsetValue, " \'{if (NR>1) {printf \"\\n\"}; for (i=1;i<=8;i++){printf $i\"\t\"}; printf $9\"\t\"($10+t)\"\t\"($11-t)\"\t\"; for (i=12;i<=NF;i++){printf $i\"\t\"}}\' - > ", tempfile_3))

			if (bool_DF_Prom_Direct_1PorE_Indirect_EQTL == FALSE) {
				# write the interacting bins
				write.table(OutDF, Prom_Direct_1PorE_Indirect_EQTL_DumpFile, row.names=F, col.names=T, sep="\t", quote=F, append=F)
				# write the header line of TSS information file
				fp_out <- file(Prom_Direct_1PorE_Indirect_EQTL_Dump_TSSFile, "w")
				writeLines(paste(c("TSSc", "TSSs", "TSSe", "TSSIdx", "Strand", "GeneID", "GeneType", "GeneName", CN, eqtl_header_info), collapse="\t"), con=fp_out)
				close(fp_out)
				# set the boolean flag
				bool_DF_Prom_Direct_1PorE_Indirect_EQTL <- TRUE					
			} else {
				write.table(OutDF, Prom_Direct_1PorE_Indirect_EQTL_DumpFile, row.names=F, col.names=F, sep="\t", quote=F, append=T)
			}
			# append the TSS file information
			system(paste("cat", tempfile_3, ">>", Prom_Direct_1PorE_Indirect_EQTL_Dump_TSSFile))
			# remove temporary files
			system(paste("rm", tempfile_1))
			system(paste("rm", tempfile_2))	
			system(paste("rm", tempfile_3))
		}				
	}
	
	#*****************
	rm(Label_BinnedSeg_CurrChr)	
	rm(InpGraph)
	rm(ShortDist)
	rm(Conn_Idx_DF_Dist1)
	rm(Conn_Idx_DF_Dist2)
	rm(Conn_Idx_DF_Dist3)
	rm(NodeNameDF)
	rm(EdgesDF)
	rm(Prom_Nodes_Idx)
	rm(Enh_Nodes_Idx)
	rm(ov1)
	rm(ov2)
	rm(AnnotLoopData_CurrChr)

	rm(EQTL_Nodes_Idx)		
	rm(EQTLData_CurrChr)
	rm(OnlyEQTLInfo_CurrChr)

	# also call garbage collector routine
	gc()
	#*****************

#================================
}	# end chromosome loop
#================================

# remove temporary files
system(paste("rm", BinnedSeg_temp_File))
system(paste("rm", BinnedSegFile))
system(paste("rm", AnnotLoopFile_CurrChr))
system(paste("rm", EQTLFile_CurrChr))
system(paste("rm", GTFFile_CurrChr))
system(paste("rm", GTFFile_GeneExprFile_CurrChr))


##=======================
# finally extract the direct and indirect pieQTLs
##=======================

##======== extract direct pieQTLs - in the format of FitHiChIP interactions
system(paste0("awk -F\'[\t]\' \'function abs(x) {return x < 0 ? -x : x} {if ((NR==1) || (($6==$17) && (abs($2-$16)>10000))) {print $0}}\' ", Prom_Direct_EQTL_Dump_TSSFile, " > ", direct_pieQTL_loop_file))

##======== extract indirect pieQTLs - in the format of FitHiChIP interactions
system(paste0("awk -F\'[\t]\' \'function abs(x) {return x < 0 ? -x : x} {if ((NR==1) || (($6==$20) && (abs($2-$19)>10000))) {print $0}}\' ", Prom_Direct_1PorE_Indirect_EQTL_Dump_TSSFile, " > ", indirect_pieQTL_loop_file))

##======== now create a consolidated list containing direct and indirect pieQTLs
# fields: 'pieQTL.ID', 'Chromosome', 'pieQTL.Position', 'target_GeneID', 'target_GeneName', 'TSS', 'pvalue', 'FDR', 'beta', 'ref', 'alt', 'Interaction_Type'
dumped_pieQTL_File <- 'dumped_pieQTL.txt'
system(paste0("awk -F\'[\t]\' \'{if (NR>1) {print $18\"\t\"$1\"\t\"$16\"\t\"$6\"\t\"$8\"\t\"$2\"\t\"$19\"\t\"$20\"\t\"$21\"\t\"$22\"\t\"$23\"\tDirect_PieQTL\"}}\' ", direct_pieQTL_loop_file, " > ", dumped_pieQTL_File))
system(paste0("awk -F\'[\t]\' \'{if (NR>1) {print $21\"\t\"$1\"\t\"$19\"\t\"$6\"\t\"$8\"\t\"$2\"\t\"$22\"\t\"$23\"\t\"$24\"\t\"$25\"\t\"$26\"\tIndirect_PieQTL\"}}\' ", indirect_pieQTL_loop_file, " >> ", dumped_pieQTL_File))

##======== insert the column names and write back to the same file
orig_pieQTL_data <- data.table::fread(dumped_pieQTL_File, header=F)
colnames(orig_pieQTL_data) <- c('pieQTL.ID', 'Chromosome', 'pieQTL.Position', 'target_GeneID', 'target_GeneName', 'TSS', 'pvalue', 'FDR', 'beta', 'ref', 'alt', 'Interaction_Type')
cat(sprintf("\n\n *** dumped pieQTL (direct / indirect) - number of entries : %s ", nrow(orig_pieQTL_data)))
orig_pieQTL_data <- unique(orig_pieQTL_data)	# there can be duplicate entries - so filter them
cat(sprintf("\n\n *** dumped pieQTL (direct / indirect) after removing duplicates - number of entries : %s ", nrow(orig_pieQTL_data)))
write.table(orig_pieQTL_data, dumped_pieQTL_File, row.names=F, col.names=T, sep="\t", quote=F, append=F)

##====== important: discard the indirect pieQTLs of those genes which have at least one direct pieQTL 	
# list of genes having direct pieQTL
genes_direct_pieQTL <- as.vector(unique(orig_pieQTL_data[which(orig_pieQTL_data[, ncol(orig_pieQTL_data)] == "Direct_PieQTL"), 4]))
cat(sprintf("\n number of genes having direct pieQTL : %s ", length(genes_direct_pieQTL)))

# list of genes having indirect pieQTL
genes_indirect_pieQTL <- as.vector(unique(orig_pieQTL_data[which(orig_pieQTL_data[, ncol(orig_pieQTL_data)] == "Indirect_PieQTL"), 4]))
cat(sprintf("\n number of genes having indirect pieQTL : %s ", length(genes_indirect_pieQTL)))

# genes having both direct and indirect pieQTLs
genes_both_direct_indirect_pieQTL <- intersect(genes_direct_pieQTL, genes_indirect_pieQTL)
cat(sprintf("\n number of genes having both direct and indirect pieQTL : %s ", length(genes_both_direct_indirect_pieQTL)))

if (length(genes_both_direct_indirect_pieQTL) > 0) {
	rows_to_remove <- c()
	for (geneidx in 1:length(genes_both_direct_indirect_pieQTL)) {
		currgene <- genes_both_direct_indirect_pieQTL[geneidx]
		rows_to_remove <- c(rows_to_remove, which((orig_pieQTL_data[, 4] == currgene) & (orig_pieQTL_data[, ncol(orig_pieQTL_data)] == "Indirect_PieQTL")))
	}
	rows_to_remove <- unique(rows_to_remove)
	cat(sprintf("\n *** number of rows (indirect pieQTL - gene pair) to remove : %s ", length(rows_to_remove)))

	finalDF <- orig_pieQTL_data[-c(rows_to_remove), ]
	cat(sprintf("\n Modified pieQTL data entries : %s ", nrow(finalDF)))
} else {
	finalDF <- orig_pieQTL_data
	cat(sprintf("\n No filtering was required - final pieQTL data entries : %s ", nrow(finalDF)))
}
write.table(finalDF, filtered_pieQTL_File, row.names=F, col.names=T, sep="\t", quote=F, append=F)

