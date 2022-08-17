#!/usr/bin/env Rscript

source('Header.r')
source('Functions.r')

Color_Types <- c('deepskyblue', 'violet', 'red', 'green', 'gold')

##=== user needs to edit this
##=== list of donor specific FitHiChIP significant files (FDR < 0.01)
##=== each of these loops are assumed to be of 5 Kb bin size
##=== distance range = 10 Kb to 3 Mb
##=== and derived by peak-to-all interaction model with loose (P2P=0) background
CD4N_DonorList <- c('CD4N_R1.interactions_FitHiC_Q0.01.bed', 'CD4N_R2.interactions_FitHiC_Q0.01.bed', 'CD4N_R3.interactions_FitHiC_Q0.01.bed', 'CD4N_R4.interactions_FitHiC_Q0.01.bed', 'CD4N_R5.interactions_FitHiC_Q0.01.bed', 'CD4N_R6.interactions_FitHiC_Q0.01.bed')

CD8N_DonorList <- c('CD8N_R1.interactions_FitHiC_Q0.01.bed', 'CD8N_R2.interactions_FitHiC_Q0.01.bed', 'CD8N_R3.interactions_FitHiC_Q0.01.bed', 'CD8N_R4.interactions_FitHiC_Q0.01.bed', 'CD8N_R5.interactions_FitHiC_Q0.01.bed', 'CD8N_R6.interactions_FitHiC_Q0.01.bed')

Mono_DonorList <- c('Mono_R1.interactions_FitHiC_Q0.01.bed', 'Mono_R2.interactions_FitHiC_Q0.01.bed', 'Mono_R3.interactions_FitHiC_Q0.01.bed', 'Mono_R4.interactions_FitHiC_Q0.01.bed', 'Mono_R5.interactions_FitHiC_Q0.01.bed', 'Mono_R6.interactions_FitHiC_Q0.01.bed')

NB_DonorList <- c('NB_R1.interactions_FitHiC_Q0.01.bed', 'NB_R2.interactions_FitHiC_Q0.01.bed', 'NB_R3.interactions_FitHiC_Q0.01.bed', 'NB_R4.interactions_FitHiC_Q0.01.bed', 'NB_R5.interactions_FitHiC_Q0.01.bed', 'NB_R6.interactions_FitHiC_Q0.01.bed')

NK_DonorList <- c('NK_R1.interactions_FitHiC_Q0.01.bed', 'NK_R2.interactions_FitHiC_Q0.01.bed', 'NK_R3.interactions_FitHiC_Q0.01.bed', 'NK_R4.interactions_FitHiC_Q0.01.bed', 'NK_R5.interactions_FitHiC_Q0.01.bed', 'NK_R6.interactions_FitHiC_Q0.01.bed')

#============================
# this function generates various plots and statistics for PCA
#============================
Dump_Plot_PCA <- function(inpPCARes, outprefix, GroupVec, ColorVec, plotWidth=7, plotHeight=6) {

	# eigen values - three columns: 1) eigenvalues 2) % of variances 3) cumulative % of variance
	eig <- factoextra::get_eig(inpPCARes)
	# contribution of individuals (replicates)
	ind <- factoextra::get_pca_ind(inpPCARes)
	
	plotfile <- paste0(outprefix, '_viz_pca_ind_Plot_LATEST.pdf')
	factoextra::fviz_pca_ind(inpPCARes, pointsize = 4, col.ind = GroupVec, repel = TRUE, select.ind=rownames(inpPCARes)) + scale_color_manual(values=ColorVec) + scale_shape_manual(values=rep(1, length(ColorVec))) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	ggsave(plotfile, width=plotWidth, height=plotHeight)

	# used latest - without text label
	plotfile <- paste0(outprefix,'_viz_pca_ind_Plot_LATEST_2.pdf')
	factoextra::fviz_pca_ind(inpPCARes, pointsize = 4, geom=c("point"), col.ind = GroupVec, repel = TRUE, select.ind=rownames(inpPCARes)) + scale_color_manual(values=ColorVec) + scale_shape_manual(values=rep(1, length(ColorVec))) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	ggsave(plotfile, width=plotWidth-2, height=plotHeight-2)

}	# end function

#==========================================
# creation of feature values for PCA - top K loops from each donor, where K = 20000
#==========================================

K <- 20000

UnionLoopTempFile <- 'Union_Loops_Donors.bed'
UnionLoopFile <- 'Union_Loops_Donors_Sorted.bed'

for (CellIdx in (1:length(Cell_Types))) {
	cat(sprintf("\n ===>>>> Union loop creation -- Processing cell type : %s ", Cell_Types[CellIdx]))
	if (CellIdx == 1) {
		LoopDirReplVec <- Mono_DonorList
	} else if (CellIdx == 2) {
		LoopDirReplVec <- NK_DonorList
	} else if (CellIdx == 3) {
		LoopDirReplVec <- NB_DonorList
	} else if (CellIdx == 4) {
		LoopDirReplVec <- CD4N_DonorList
	} else {
		LoopDirReplVec <- CD8N_DonorList
	}	

	for (rc in (1:length(LoopDirReplVec))) {				
		cat(sprintf("\n ===>>>> Processing donor count : %s ", rc))
		ReplLoopFile <- LoopDirReplVec[rc]
		nf <- GetNumFields(ReplLoopFile)
		if ((CellIdx == 1) & (rc == 1)) {
			system(paste0("awk \'(NR>1)\' ", ReplLoopFile, " | sort -k", nf, ",", nf, "g - | head -n ", K, " | cut -f1-6 > ", UnionLoopTempFile))
		} else {
			system(paste0("awk \'(NR>1)\' ", ReplLoopFile, " | sort -k", nf, ",", nf, "g - | head -n ", K, " | cut -f1-6 >> ", UnionLoopTempFile))
		}
	}	# end replicate loop	
}	# end cell type loop

# now sort the generated loops and extract the unique ones
system(paste("sort -k1,1 -k2,2n -k5,5n", UnionLoopTempFile, " | uniq >", UnionLoopFile))

# final feature containing file (for all donors)
FinalFeatureFile <- 'Final_Features.bed'

# boolean: construction of donor specific feature data (chromosome wise)
bool_Donor_ChrWise <- FALSE

for (chr_idx in (1:length(ChrList_NameNum))) {
	chrName <- ChrList_NameNum[chr_idx]
	cat(sprintf("\n\n\n **** Processing chromosome : %s ", chrName))

	# extract the merged loops with respect to current chromosome
	InpTempFile <- paste0(dirname(UnionLoopFile), '/InpFile_Temp.bed')
	ExtractChrData(UnionLoopFile, chrName, InpTempFile, header=FALSE)
	cat(sprintf("\n === Extracted union loop data for the current chromosome "))
	# check if there is no loop for the current chromosome - then continue
	nreadCurr <- GetNumLines(InpTempFile)
	if (nreadCurr == 0) {
		next
	}
	# read the current chromosome specific loops
	CurrChrData <- data.table::fread(InpTempFile, header=F)
	CurrChrData <- CurrChrData[,1:6]

	CN <- c("chr1", "start1", "end1", "chr2", "start2", "end2")

	for (CellIdx in (1:length(Cell_Types))) {
		cat(sprintf("\n\n ===>>>> Processing cell type : %s ", Cell_Types[CellIdx]))
		# list of directories containing FitHiChIP loops
		if (CellIdx == 1) {
			LoopDirReplVec <- Mono_DonorList
		} else if (CellIdx == 2) {
			LoopDirReplVec <- NK_DonorList
		} else if (CellIdx == 3) {
			LoopDirReplVec <- NB_DonorList
		} else if (CellIdx == 4) {
			LoopDirReplVec <- CD4N_DonorList
		} else {
			LoopDirReplVec <- CD8N_DonorList
		}	

		for (rc in (1:length(LoopDirReplVec))) {				
			cat(sprintf("\n ===>>>> Processing replica count : %s ", rc))
			ReplLoopFile <- LoopDirReplVec[rc]

			CCVec <- rep(0, nrow(CurrChrData))
			qVec <- rep(-1, nrow(CurrChrData))
				
			ReplTempFile <- paste0(dirname(UnionLoopFile), '/ReplFile_Temp.bed')
			ExtractChrData(ReplLoopFile, chrName, ReplTempFile, header=TRUE)
			cat(sprintf("\n === Extracted replicate loop data for the current chromosome "))
			nreadRepl <- GetNumLines(ReplTempFile)
			if (nreadRepl == 0) {
				CurrChrData <- cbind.data.frame(CurrChrData, CCVec, qVec)
				next
			}
			ReplLoopData <- data.table::fread(ReplTempFile, header=F)
			ReplLoopData <- ReplLoopData[,c(1:7, ncol(ReplLoopData))]
			CurrOv <- OverlapLoop(CurrChrData[,1:6], ReplLoopData[,1:6], boundary=1, offset=0)
			# copy the contact count and q-value statistics
			CCVec[CurrOv$A_AND_B] <- ReplLoopData[CurrOv$B_AND_A, 7]
			qVec[CurrOv$A_AND_B] <- ReplLoopData[CurrOv$B_AND_A, 8]
			# update the data frame of the current chromosome
			CurrChrData <- cbind.data.frame(CurrChrData, CCVec, qVec)
			# update the column name vector
			CN <- c(CN, paste0(Cell_Types[CellIdx], '_', rc, '_CC'), paste0(Cell_Types[CellIdx], '_', rc, '_Qval'))

		}	# end replicate count loop
	}	# end cell index loop

	# now assign the column names of this chromosome specific data
	colnames(CurrChrData) <- CN
	if (bool_Donor_ChrWise == FALSE) {
		FinalDF <- CurrChrData
		bool_Donor_ChrWise <- TRUE
	} else {
		FinalDF <- rbind.data.frame(FinalDF, CurrChrData)
	}
}	# end chromosome loop

# write the final output
write.table(FinalDF, FinalFeatureFile, row.names = FALSE, col.names = TRUE, sep = "\t", quote=FALSE, append=FALSE)

#========================
# PCA feature vector
#========================
inpdata <- data.table::fread(FinalFeatureFile, header=T)
cln <- colnames(inpdata)

CC_ColList <- seq(7, ncol(inpdata), by=2)
Q_ColList <- seq(8, ncol(inpdata), by=2)
Lbls <- colnames(inpdata)[CC_ColList]
NewLbls <- gsub("_CC$", "", Lbls)
QValDF <- cbind.data.frame(inpdata[, Q_ColList])
colnames(QValDF) <- NewLbls

# transform the feature vector so that rows: individual replicates, and columns: different features
t_QValDF <- t(QValDF)
res.pca_1_QValDF <- prcomp(t_QValDF, scale = FALSE)
cat(sprintf("\n\n *** Dimension of res.pca_1_QValDF : %s  X  %s  ", nrow(res.pca_1_QValDF), ncol(res.pca_1_QValDF)))

# create color variable as a factor
groupvec <- as.factor(c(rep('1_Mono', 6), rep('2_NK', 6), rep('3_NB', 6), rep('4_CD4N', 6), rep('5_CD8N', 6)))

#========================
# plot PCA
#========================
Dump_Plot_PCA(res.pca_1_QValDF, paste0(dirname(FinalFeatureFile), '/PCA_1_QValDF_Reduced_Samples'), groupvec, Color_Types)
write.table(QValDF, paste0(dirname(FinalFeatureFile), '/PCA_1_QValDF_Reduced_Samples_Dumped_Data.bed'), row.names=F, col.names=T, sep="\t", append=F)

