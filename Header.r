#!/usr/bin/env Rscript

suppressMessages(library(GenomicRanges))
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(edgeR)
library(stats)
library(factoextra)
library(pheatmap)
library(dplyr)
library(parallel)
library(plotly)
library(igraph)
library(data.table)
library(ggpubr)
options(scipen = 10)
options(datatable.fread.datatable=FALSE)

HOMERExecPath <- '/home/sourya/packages/HOMER/bin/'
HOMERPeakAnnotExec <- '/home/sourya/packages/HOMER/bin/annotatePeaks.pl'
MotifFindExec <- '/home/sourya/packages/HOMER/bin/findMotifsGenome.pl'
deeptoolsPath <- '/home/sourya/packages/deepTools/deepTools2.0/bin/'

##=== we used hg19 reference genome, with gencode.v19 annotation
refgenome <- 'hg19'
hg19_fastafile <- 'hg19.fa'
# hg19_ucsc_annotationfile <- 'hg19.gtf'
hg19_annotfile <- 'gencode.v19.annotation.gtf'

##==== gene annotation files derived from the GTF input files
GTFType <- 2
InpGTFFile_ALL <- 'gEncode_Genes_Complete_NEW_TSS.gtf'
CustomGTF_AllGeneFile <- 'gEncode_Genes_Complete_type_GeneID_TranscriptID_NEW.gtf'
offset_val_ALLGene_GTF <- 2500

##==== cell types
Cell_Types <- c("Mono", "NK", "NB", "CD4N", "CD8N")
Cell_Types_DICEDB <- c("CM", "NK", "BN", "CD4N", "CD8N")

Colors <- c("#B22222", "#FFD700", "#32CD32", "#00BFFF", "#EE82EE")
# Colors <- c('red', 'yellow', 'green', 'deepskyblue', 'violet')	#c('deepskyblue', 'violet', 'orangered4', 'mediumseagreen', 'gold')	#c("red", "green", "blue", "orange", "violet")

##=== set of chromosomes for hg19
ChrList_NameNum <- c(paste("chr", seq(1,22), sep=""), "chrX", "chrY")



































#=====================================
# gene expression related files
#=====================================

# gene expression file for different cell types 
# this file is cell specific and condensed
# GeneExprFile <- '/mnt/BioAdHoc/Groups/vd-vijay/sourya/DATA/HiChIP/Merged_Vivek_Feb_March_2018/Reference_GTF_Imp/Gene_Expression_Cell_Specific/DERIVED_FILES/Cell_Type_Expression_Summary.bed'
GeneExprFile <- paste0(MAINBASEINPDIR, 'Data/Reference_GTF_Imp/Gene_Expression_All/Mean_TPM_5celltypes_Summary.txt')

# file list storing the EQTL information in cell specific manner
EQTL_FileList <- c(paste0(MAINBASEINPDIR, 'Data/Reference_GTF_Imp/EQTL/Mono/CM_eQTLs_FDR_0.05_SHORT.bed'), paste0(MAINBASEINPDIR, 'Data/Reference_GTF_Imp/EQTL/NK/NK_eQTLs_FDR_0.05_SHORT.bed'), paste0(MAINBASEINPDIR, 'Data/Reference_GTF_Imp/EQTL/NB/BN_eQTLs_FDR_0.05_SHORT.bed'), paste0(MAINBASEINPDIR, 'Data/Reference_GTF_Imp/EQTL/CD4N/CD4N_eQTLs_FDR_0.05_SHORT.bed'), paste0(MAINBASEINPDIR, 'Data/Reference_GTF_Imp/EQTL/CD8N/CD8N_eQTLs_FDR_0.05_SHORT.bed'))

# file list storing the EQTL information in cell specific manner
# here eQTLs overlap with reference ChIP-seq peaks for the corresponding cell type
EQTL_Ov_ChIP_Peak_FileList <- c(paste0(MAINBASEINPDIR, 'Data/Reference_GTF_Imp/EQTL/Mono/CM_eQTLs_FDR_0.05_SHORT_Overlap_ChIP_Peaks.bed'), paste0(MAINBASEINPDIR, 'Data/Reference_GTF_Imp/EQTL/NK/NK_eQTLs_FDR_0.05_SHORT_Overlap_ChIP_Peaks.bed'), paste0(MAINBASEINPDIR, 'Data/Reference_GTF_Imp/EQTL/NB/BN_eQTLs_FDR_0.05_SHORT_Overlap_ChIP_Peaks.bed'), paste0(MAINBASEINPDIR, 'Data/Reference_GTF_Imp/EQTL/CD4N/CD4N_eQTLs_FDR_0.05_SHORT_Overlap_ChIP_Peaks.bed'), paste0(MAINBASEINPDIR, 'Data/Reference_GTF_Imp/EQTL/CD8N/CD8N_eQTLs_FDR_0.05_SHORT_Overlap_ChIP_Peaks.bed'))

#=====================================
# directory base names containing FitHiChIP outputs
#=====================================

# input directory name bases
LoopDirList <- c('Merged_Monocyte', 'Merged_NK', 'Merged_NB', 'Merged_CD4_Naive', 'Merged_CD8_Naive')


#=====================================
# base directory containing PCHIC loops for different cell types
#=====================================

PCHIC_BaseDir <- '/mnt/BioAdHoc/Groups/vd-vijay/sourya/DATA/CaptureHiC/Javierre_2016/Detected_Interactions/PCHiC/Cell_Type_Specific_Loops/'

#=====================================
# reference ChIP-seq peak files
#=====================================

# Input reference ChIP-seq files
RefChIPSeqPeakFiles <- paste0(RESULTBASEDIR, '/ChIP_seq_Peaks/', Cell_Types, '/Peak_FDR0.01')

# files containing exclusive peaks for specific cell types
Excl_Peak_Files <- paste0(RESULTBASEDIR, "/Peak_Overlap_Stat/CompletePeaks/Cell_Specific_Peaks_Overlapping_FitHiChIP_Loops/OvExcl_", Cell_Types, ".bed")

#======================================================
# Miscellaneous parameters
#======================================================

# set of chromosomes which would be individually analyzed for testing overlap


# fold change threshold for differential analysis
FOLD_Change_Thr <- 3

# FDR threshold
FDR_Th_DESeq <- 0.01

# colors corresponding to the cell type
#======================================================
# distance range and other related parameters for FitHiChIP loops
#======================================================

BINSIZE_FITHICHIP_LOOP <- 5000
LOW_DIST_RANGE_FITHICHIP_LOOP <- 10000
HIGH_DIST_RANGE_FITHICHIP_LOOP <- 3000000

#======================================================
# replicate directories
# along with their respective BAM and peak files
# including the merged replicates as well
#======================================================

NB_Rep <- c('Merged_NB', 'VC1', 'VC101', 'VC106', 'VC201', 'VC206', 'VC6')

NB_Rep_DonorID <- c('Merged', '1829-RH-1', '1831-RH-1', '1800-RH-1', '1814-RH-1', '1816-RH-1', '1815-RH-1')

NB_Rep_BAMFiles <- c('/mnt/BioAdHoc/Groups/vd-vijay/sourya/genomes/REF_ChIP_Seq_Peaks/LJI_Custom/NB/Merge_Specific_Donors/merged_input_sorted.bam', '/mnt/NGSAnalyses/ChIP-Seq/Mapping/004291_R24_NB_AF/All_Uniquely_Mapped_BAM_Files/35_C_NaveB_3_C3_1829_RH1.bam', '/mnt/NGSAnalyses/ChIP-Seq/Mapping/004291_R24_NB_AF/All_Uniquely_Mapped_BAM_Files/5_A_NaveB_5_A5_1831_RH1.bam', '/mnt/NGSAnalyses/ChIP-Seq/Mapping/004291_R24_NB_AF/All_Uniquely_Mapped_BAM_Files/2_A_NaveB_2_A2_1800_RH1.bam', '/mnt/NGSAnalyses/ChIP-Seq/Mapping/004291_R24_NB_AF/All_Uniquely_Mapped_BAM_Files/17_B_NaveB_1_B1_1814_RH1.bam', '/mnt/NGSAnalyses/ChIP-Seq/Mapping/004291_R24_NB_AF/All_Uniquely_Mapped_BAM_Files/1_A_NaveB_1_A1_1816_RH1.bam', '/mnt/NGSAnalyses/ChIP-Seq/Mapping/004291_R24_NB_AF/All_Uniquely_Mapped_BAM_Files/21_B_NaveB_5_B5_1815_RH1.bam')

NB_Rep_PeakFiles <- c(paste0(RESULTBASEDIR, '/NB/Peak_FDR0.01'), paste0(MAINBASEINPDIR, 'Session_Details/B_Cell/PeaksFromDonors/NB_1829_RH1_peaks_Q0.01filt'), paste0(MAINBASEINPDIR, 'Session_Details/B_Cell/PeaksFromDonors/NB_1831_RH1_peaks_Q0.01filt'), paste0(MAINBASEINPDIR, 'Session_Details/B_Cell/PeaksFromDonors/NB_1800_RH1_peaks_Q0.01filt'), paste0(MAINBASEINPDIR, 'Session_Details/B_Cell/PeaksFromDonors/NB_1814_RH1_peaks_Q0.01filt'), paste0(MAINBASEINPDIR, 'Session_Details/B_Cell/PeaksFromDonors/NB_1816_RH1_peaks_Q0.01filt'), paste0(MAINBASEINPDIR, 'Session_Details/B_Cell/PeaksFromDonors/NB_1815_RH1_peaks_Q0.01filt'))


CD4N_Rep <- c('Merged_CD4_Naive', 'VC102', 'VC107', 'VC11', 'VC111', 'VC2', 'VC202', 'VC207', 'VC212', 'VC7')
CD4N_Rep_Reduced <- c('Merged_CD4_Naive', 'VC102', 'VC107', 'VC2', 'VC202', 'VC207', 'VC7')

CD4N_Rep_DonorID <- c('Merged', '1831-RH-1', '1800-RH-1', '1831-RH-1', '1815-RH-1', '1829-RH-1', '1814-RH-1', '1816-RH-1', '1829-RH-1', '1815-RH-1')

CD4N_Rep_BAMFiles <- c('/mnt/BioAdHoc/Groups/vd-vijay/sourya/genomes/REF_ChIP_Seq_Peaks/LJI_Custom/CD4N/Merge_Specific_Donors/merged_input_sorted.bam', '/mnt/NGSAnalyses/ChIP-Seq/Mapping/003951_Log7_DY_R24_CD4N_ChIP_Merged/All_Uniquely_Mapped_BAM_Files/A_CD4N_5_A5_1831_RH1_Final.bam', '/mnt/NGSAnalyses/ChIP-Seq/Mapping/003951_Log7_DY_R24_CD4N_ChIP_Merged/All_Uniquely_Mapped_BAM_Files/A_CD4N_2_A2_1800_RH1_Final.bam', '/mnt/NGSAnalyses/ChIP-Seq/Mapping/003951_Log7_DY_R24_CD4N_ChIP_Merged/All_Uniquely_Mapped_BAM_Files/A_CD4N_5_A5_1831_RH1_Final.bam', '/mnt/NGSAnalyses/ChIP-Seq/Mapping/003951_Log7_DY_R24_CD4N_ChIP_Merged/All_Uniquely_Mapped_BAM_Files/B_CD4N_5_B5_1815_RH1_Final.bam', '/mnt/NGSAnalyses/ChIP-Seq/Mapping/003951_Log7_DY_R24_CD4N_ChIP_Merged/All_Uniquely_Mapped_BAM_Files/C_CD4N_3_C3_1829_RH1_Final.bam', '/mnt/NGSAnalyses/ChIP-Seq/Mapping/003951_Log7_DY_R24_CD4N_ChIP_Merged/All_Uniquely_Mapped_BAM_Files/B_CD4N_1_B1_1814_RH1_Final.bam', '/mnt/NGSAnalyses/ChIP-Seq/Mapping/003951_Log7_DY_R24_CD4N_ChIP_Merged/All_Uniquely_Mapped_BAM_Files/A_CD4N_1_A1_1816_RH1_Final.bam', '/mnt/NGSAnalyses/ChIP-Seq/Mapping/003951_Log7_DY_R24_CD4N_ChIP_Merged/All_Uniquely_Mapped_BAM_Files/C_CD4N_3_C3_1829_RH1_Final.bam', '/mnt/NGSAnalyses/ChIP-Seq/Mapping/003951_Log7_DY_R24_CD4N_ChIP_Merged/All_Uniquely_Mapped_BAM_Files/B_CD4N_5_B5_1815_RH1_Final.bam')

CD4N_Rep_PeakFiles <- c(paste0(RESULTBASEDIR, '/CD4N/Peak_FDR0.01'), paste0(MAINBASEINPDIR, 'Session_Details/CD4_Naive/PeaksFromDonors/CD4N_1831_RH1_peaks_Q0.01filt'), paste0(MAINBASEINPDIR, 'Session_Details/CD4_Naive/PeaksFromDonors/CD4N_1800_RH1_peaks_Q0.01filt'), paste0(MAINBASEINPDIR, 'Session_Details/CD4_Naive/PeaksFromDonors/CD4N_1831_RH1_peaks_Q0.01filt'), paste0(MAINBASEINPDIR, 'Session_Details/CD4_Naive/PeaksFromDonors/CD4N_1815_RH1_peaks_Q0.01filt'), paste0(MAINBASEINPDIR, 'Session_Details/CD4_Naive/PeaksFromDonors/CD4N_1829_RH1_peaks_Q0.01filt'), paste0(MAINBASEINPDIR, 'Session_Details/CD4_Naive/PeaksFromDonors/CD4N_1814_RH1_peaks_Q0.01filt'), paste0(MAINBASEINPDIR, 'Session_Details/CD4_Naive/PeaksFromDonors/CD4N_1816_RH1_peaks_Q0.01filt'), paste0(MAINBASEINPDIR, 'Session_Details/CD4_Naive/PeaksFromDonors/CD4N_1829_RH1_peaks_Q0.01filt'), paste0(MAINBASEINPDIR, 'Session_Details/CD4_Naive/PeaksFromDonors/CD4N_1815_RH1_peaks_Q0.01filt'))



CD8N_Rep <- c('Merged_CD8_Naive', 'VC103', 'VC108', 'VC203', 'VC208', 'VC3', 'VC8')

CD8N_Rep_DonorID <- c('Merged', '1831-RH-1', '1800-RH-1', '1814-RH-1', '1816-RH-1', '1829-RH-1', '1815-RH-1')

CD8N_Rep_BAMFiles <- c('/mnt/BioAdHoc/Groups/vd-vijay/sourya/genomes/REF_ChIP_Seq_Peaks/LJI_Custom/CD8N/Merge_Specific_Donors/merged_input_sorted.bam', '/mnt/NGSAnalyses/ChIP-Seq/Mapping/004174_Log7_DY_R24_CD8N_ChIP_Merged/All_Uniquely_Mapped_BAM_Files/005_A_CD8N_5_A5_1831_RH1_Final.bam', '/mnt/NGSAnalyses/ChIP-Seq/Mapping/004174_Log7_DY_R24_CD8N_ChIP_Merged/All_Uniquely_Mapped_BAM_Files/002_A_CD8N_2_A2_1800_RH1_Final.bam', '/mnt/NGSAnalyses/ChIP-Seq/Mapping/004174_Log7_DY_R24_CD8N_ChIP_Merged/All_Uniquely_Mapped_BAM_Files/017_B_CD8N_1_B1_1814_RH1_Final.bam', '/mnt/NGSAnalyses/ChIP-Seq/Mapping/004174_Log7_DY_R24_CD8N_ChIP_Merged/All_Uniquely_Mapped_BAM_Files/001_A_CD8N_1_A1_1816_RH1_A_Final.bam', '/mnt/NGSAnalyses/ChIP-Seq/Mapping/004174_Log7_DY_R24_CD8N_ChIP_Merged/All_Uniquely_Mapped_BAM_Files/035_C_CD8N_3_C3_1829_RH1_Final.bam', '/mnt/NGSAnalyses/ChIP-Seq/Mapping/004174_Log7_DY_R24_CD8N_ChIP_Merged/All_Uniquely_Mapped_BAM_Files/021_B_CD8N_5_B5_1815_RH1_Final.bam')

CD8N_Rep_PeakFiles <- c(paste0(RESULTBASEDIR, '/CD8N/Peak_FDR0.01'), paste0(MAINBASEINPDIR, 'Session_Details/CD8_Naive/PeaksFromDonors/CD8N_1831_RH1_peaks_Q0.01filt'), paste0(MAINBASEINPDIR, 'Session_Details/CD8_Naive/PeaksFromDonors/CD8N_1800_RH1_peaks_Q0.01filt'), paste0(MAINBASEINPDIR, 'Session_Details/CD8_Naive/PeaksFromDonors/CD8N_1814_RH1_peaks_Q0.01filt'), paste0(MAINBASEINPDIR, 'Session_Details/CD8_Naive/PeaksFromDonors/CD8N_1816_RH1_peaks_Q0.01filt'), paste0(MAINBASEINPDIR, 'Session_Details/CD8_Naive/PeaksFromDonors/CD8N_1829_RH1_peaks_Q0.01filt'), paste0(MAINBASEINPDIR, 'Session_Details/CD8_Naive/PeaksFromDonors/CD8N_1815_RH1_peaks_Q0.01filt'))




Mono_Rep <- c('Merged_Monocyte', 'VC104', 'VC109', 'VC112', 'VC12', 'VC204', 'VC209', 'VC211', 'VC4', 'VC9')
Mono_Rep_Reduced <- c('Merged_Monocyte', 'VC104', 'VC109', 'VC204', 'VC209', 'VC4', 'VC9')

Mono_Rep_DonorID <- c('Merged', '1831-RH-1', '1800-RH-1', '1814-RH-1', '1816-RH-1', '1814-RH-1', '1816-RH-1', '1800-RH-1', '1829-RH-1', '1815-RH-1')

Mono_Rep_BAMFiles <- c('/mnt/BioAdHoc/Groups/vd-vijay/sourya/genomes/REF_ChIP_Seq_Peaks/LJI_Custom/Mono/Merge_Specific_Donors/merged_input_sorted.bam', '/mnt/NGSAnalyses/ChIP-Seq/Mapping/004020_R24_ClassMONO_RunsAF/All_Uniquely_Mapped_BAM_Files/A_ClassMCs_5_A5_1831_RH1.bam', '/mnt/NGSAnalyses/ChIP-Seq/Mapping/004020_R24_ClassMONO_RunsAF/All_Uniquely_Mapped_BAM_Files/A_ClassMCs_2_A2_1800_RH1.bam', '/mnt/NGSAnalyses/ChIP-Seq/Mapping/004020_R24_ClassMONO_RunsAF/All_Uniquely_Mapped_BAM_Files/B_ClassMCs_1_B1_1814_RH1.bam', '/mnt/NGSAnalyses/ChIP-Seq/Mapping/004020_R24_ClassMONO_RunsAF/All_Uniquely_Mapped_BAM_Files/A_ClassMCs_1_A1_1816_RH1.bam', '/mnt/NGSAnalyses/ChIP-Seq/Mapping/004020_R24_ClassMONO_RunsAF/All_Uniquely_Mapped_BAM_Files/B_ClassMCs_1_B1_1814_RH1.bam', '/mnt/NGSAnalyses/ChIP-Seq/Mapping/004020_R24_ClassMONO_RunsAF/All_Uniquely_Mapped_BAM_Files/A_ClassMCs_1_A1_1816_RH1.bam', '/mnt/NGSAnalyses/ChIP-Seq/Mapping/004020_R24_ClassMONO_RunsAF/All_Uniquely_Mapped_BAM_Files/A_ClassMCs_2_A2_1800_RH1.bam', '/mnt/NGSAnalyses/ChIP-Seq/Mapping/004020_R24_ClassMONO_RunsAF/All_Uniquely_Mapped_BAM_Files/C_ClassMCs_3_C3_1829_RH1.bam', '/mnt/NGSAnalyses/ChIP-Seq/Mapping/004020_R24_ClassMONO_RunsAF/All_Uniquely_Mapped_BAM_Files/B_ClassMCs_5_B5_1815_RH1.bam')

Mono_Rep_PeakFiles <- c(paste0(RESULTBASEDIR, '/Mono/Peak_FDR0.01'), paste0(MAINBASEINPDIR, 'Session_Details/Monocyte/PeaksFromDonors/Mono_1831_RH1_peaks_Q0.01filt'), paste0(MAINBASEINPDIR, 'Session_Details/Monocyte/PeaksFromDonors/Mono_1800_RH1_peaks_Q0.01filt'), paste0(MAINBASEINPDIR, 'Session_Details/Monocyte/PeaksFromDonors/Mono_1814_RH1_peaks_Q0.01filt'), paste0(MAINBASEINPDIR, 'Session_Details/Monocyte/PeaksFromDonors/Mono_1816_RH1_peaks_Q0.01filt'), paste0(MAINBASEINPDIR, 'Session_Details/Monocyte/PeaksFromDonors/Mono_1814_RH1_peaks_Q0.01filt'), paste0(MAINBASEINPDIR, 'Session_Details/Monocyte/PeaksFromDonors/Mono_1816_RH1_peaks_Q0.01filt'), paste0(MAINBASEINPDIR, 'Session_Details/Monocyte/PeaksFromDonors/Mono_1800_RH1_peaks_Q0.01filt'), paste0(MAINBASEINPDIR, 'Session_Details/Monocyte/PeaksFromDonors/Mono_1829_RH1_peaks_Q0.01filt'), paste0(MAINBASEINPDIR, 'Session_Details/Monocyte/PeaksFromDonors/Mono_1815_RH1_peaks_Q0.01filt'))


NK_Rep <- c('Merged_NK', 'VC10', 'VC105', 'VC110', 'VC205', 'VC210', 'VC5')

NK_Rep_DonorID <- c('Merged', '1815-RH-1', '1831-RH-1', '1800-RH-1', '1814-RH-1', '1816-RH-1', '1829-RH-1')

NK_Rep_BAMFiles <- c('/mnt/BioAdHoc/Groups/vd-vijay/sourya/genomes/REF_ChIP_Seq_Peaks/LJI_Custom/NK/Merge_Specific_Donors/merged_input_sorted.bam', '/mnt/NGSAnalyses/ChIP-Seq/Mapping/003939_Log7_DY_R24_NK_ChIP_Merged/All_Uniquely_Mapped_BAM_Files/021_B_NK_5_B5_1815_RH1_Final.bam', '/mnt/NGSAnalyses/ChIP-Seq/Mapping/003939_Log7_DY_R24_NK_ChIP_Merged/All_Uniquely_Mapped_BAM_Files/005_A_NK_5_A5_1831_RH1_Final.bam', '/mnt/NGSAnalyses/ChIP-Seq/Mapping/003939_Log7_DY_R24_NK_ChIP_Merged/All_Uniquely_Mapped_BAM_Files/002_A_NK_2_A2_1800_RH1_Final.bam', '/mnt/NGSAnalyses/ChIP-Seq/Mapping/003939_Log7_DY_R24_NK_ChIP_Merged/All_Uniquely_Mapped_BAM_Files/017_B_NK_1_B1_1814_RH1_Final.bam', '/mnt/NGSAnalyses/ChIP-Seq/Mapping/003939_Log7_DY_R24_NK_ChIP_Merged/All_Uniquely_Mapped_BAM_Files/001_A_NK_1_A1_1816_RH1_Final.bam', '/mnt/NGSAnalyses/ChIP-Seq/Mapping/003939_Log7_DY_R24_NK_ChIP_Merged/All_Uniquely_Mapped_BAM_Files/035_C_NK_3_C3_1829_RH1_Final.bam')

NK_Rep_PeakFiles <- c(paste0(RESULTBASEDIR, '/NK/Peak_FDR0.01'), paste0(MAINBASEINPDIR, 'Session_Details/NK/PeaksFromDonors/NK_1815_RH1_peaks_Q0.01filt'), paste0(MAINBASEINPDIR, 'Session_Details/NK/PeaksFromDonors/NK_1831_RH1_peaks_Q0.01filt'), paste0(MAINBASEINPDIR, 'Session_Details/NK/PeaksFromDonors/NK_1800_RH1_peaks_Q0.01filt'), paste0(MAINBASEINPDIR, 'Session_Details/NK/PeaksFromDonors/NK_1814_RH1_peaks_Q0.01filt'), paste0(MAINBASEINPDIR, 'Session_Details/NK/PeaksFromDonors/NK_1816_RH1_peaks_Q0.01filt'), paste0(MAINBASEINPDIR, 'Session_Details/NK/PeaksFromDonors/NK_1829_RH1_peaks_Q0.01filt'))

#*************************************
# end global variable declaration
#*************************************

# CD4N ChIP-seq bigwig files for different proteins or histone modifications of interest

# comment - sourya
# when only one replicate and corresponding bigwig files were considered
# CD4N_ChIP_BigWig_FileList <- c('/mnt/BioAdHoc/Groups/vd-vijay/sourya/genomes/REF_ChIP_Seq_Peaks/LJI_Custom/CD4N/Merge_Specific_Donors/merged_CD4_Naive_ChIPSeq.bw', paste0(RESULTBASEDIR, '/Blood_Sample_Reference_Peaks/CD4N_Primary/28977.Roadmap.SRS255301_101_8_pooled_leukopaks_Jan_20_2011_716.H3K27me3.signal.bigWig'), paste0(RESULTBASEDIR, '/Blood_Sample_Reference_Peaks/CD4N_Primary/29013.Roadmap.SRS255301_101_8_pooled_leukopaks_Jan_20_2011_716.H3K4me3.signal.bigWig'), paste0(RESULTBASEDIR, '/Blood_Sample_Reference_Peaks/CD4N_Primary/29048.Roadmap.SRS255304_100_7_pooled_leukopaks_Jan_7_2011_716.H3K9me3.signal.bigWig'), paste0(RESULTBASEDIR, '/Blood_Sample_Reference_Peaks/CD4N_Primary/29142.Roadmap.SRS255301_101_8_pooled_leukopaks_Jan_20_2011_716.H3K36me3.signal.bigWig'), paste0(RESULTBASEDIR, '/Blood_Sample_Reference_Peaks/CD4N_Primary/29285.Roadmap.SRS255301_101_8_pooled_leukopaks_Jan_20_2011_716.H3K4me1.signal.bigWig'))

CD4N_ChIP_BigWig_FileList <- c('/mnt/BioAdHoc/Groups/vd-vijay/sourya/genomes/REF_ChIP_Seq_Peaks/LJI_Custom/CD4N/Merge_Specific_Donors/merged_CD4_Naive_ChIPSeq.bw', paste0(RESULTBASEDIR, '/Blood_Sample_Reference_Peaks/CD4N_Primary/merged_ChIP_H3K27me3.bigWig'), paste0(RESULTBASEDIR, '/Blood_Sample_Reference_Peaks/CD4N_Primary/merged_ChIP_H3K4me3.bigWig'), paste0(RESULTBASEDIR, '/Blood_Sample_Reference_Peaks/CD4N_Primary/merged_ChIP_H3K9me3.bigWig'), paste0(RESULTBASEDIR, '/Blood_Sample_Reference_Peaks/CD4N_Primary/merged_ChIP_H3K36me3.bigWig'), paste0(RESULTBASEDIR, '/Blood_Sample_Reference_Peaks/CD4N_Primary/merged_ChIP_H3K4me1.bigWig'))


# CD8N ChIP-seq bigwig files for different proteins or histone modifications of interest

# comment - sourya
# when only one replicate and corresponding bigwig files were considered
# CD8N_ChIP_BigWig_FileList <- c('/mnt/BioAdHoc/Groups/vd-vijay/sourya/genomes/REF_ChIP_Seq_Peaks/LJI_Custom/CD8N/Merge_Specific_Donors/merged_CD8_Naive_ChIPSeq.bw', paste0(RESULTBASEDIR, '/Blood_Sample_Reference_Peaks/CD8N_Primary/29143.Roadmap.SRS255298_100_7_pooled_leukopaks_Jan_7_2011_540.H3K27me3.signal.bigWig'), paste0(RESULTBASEDIR, '/Blood_Sample_Reference_Peaks/CD8N_Primary/29181.Roadmap.SRS255297_101_8_pooled_leukopaks_Jan_20_2011_540.H3K4me3.signal.bigWig'), paste0(RESULTBASEDIR, '/Blood_Sample_Reference_Peaks/CD8N_Primary/29377.Roadmap.SRS255298_100_7_pooled_leukopaks_Jan_7_2011_540.H3K9me3.signal.bigWig'), paste0(RESULTBASEDIR, '/Blood_Sample_Reference_Peaks/CD8N_Primary/29313.Roadmap.SRS255298_100_7_pooled_leukopaks_Jan_7_2011_540.H3K36me3.signal.bigWig'), paste0(RESULTBASEDIR, '/Blood_Sample_Reference_Peaks/CD8N_Primary/29252.Roadmap.SRS255298_100_7_pooled_leukopaks_Jan_7_2011_540.H3K4me1.signal.bigWig'))

CD8N_ChIP_BigWig_FileList <- c('/mnt/BioAdHoc/Groups/vd-vijay/sourya/genomes/REF_ChIP_Seq_Peaks/LJI_Custom/CD8N/Merge_Specific_Donors/merged_CD8_Naive_ChIPSeq.bw', paste0(RESULTBASEDIR, '/Blood_Sample_Reference_Peaks/CD8N_Primary/merged_ChIP_H3K27me3.bigWig'), paste0(RESULTBASEDIR, '/Blood_Sample_Reference_Peaks/CD8N_Primary/merged_ChIP_H3K4me3.bigWig'), paste0(RESULTBASEDIR, '/Blood_Sample_Reference_Peaks/CD8N_Primary/merged_ChIP_H3K9me3.bigWig'), paste0(RESULTBASEDIR, '/Blood_Sample_Reference_Peaks/CD8N_Primary/merged_ChIP_H3K36me3.bigWig'), paste0(RESULTBASEDIR, '/Blood_Sample_Reference_Peaks/CD8N_Primary/merged_ChIP_H3K4me1.bigWig'))


# NB ChIP-seq bigwig files for different proteins or histone modifications of interest

# comment - sourya
# when only one replicate and corresponding bigwig files were considered
# NB_ChIP_BigWig_FileList <- c('/mnt/BioAdHoc/Groups/vd-vijay/sourya/genomes/REF_ChIP_Seq_Peaks/LJI_Custom/NB/Merge_Specific_Donors/merged_B_Cell_ChIPSeq_sorted.bw', paste0(RESULTBASEDIR, '/Blood_Sample_Reference_Peaks/CD38_Neg_NB/25059.Blueprint.ERS206573.H3K27me3.signal.bigWig'), paste0(RESULTBASEDIR, '/Blood_Sample_Reference_Peaks/CD38_Neg_NB/25173.Blueprint.ERS353040.H3K4me3.signal.bigWig'), paste0(RESULTBASEDIR, '/Blood_Sample_Reference_Peaks/CD38_Neg_NB/25062.Blueprint.ERS206573.H3K9me3.signal.bigWig'), paste0(RESULTBASEDIR, '/Blood_Sample_Reference_Peaks/CD38_Neg_NB/25171.Blueprint.ERS353040.H3K36me3.signal.bigWig'), paste0(RESULTBASEDIR, '/Blood_Sample_Reference_Peaks/CD38_Neg_NB/25172.Blueprint.ERS353040.H3K4me1.signal.bigWig'))

NB_ChIP_BigWig_FileList <- c('/mnt/BioAdHoc/Groups/vd-vijay/sourya/genomes/REF_ChIP_Seq_Peaks/LJI_Custom/NB/Merge_Specific_Donors/merged_B_Cell_ChIPSeq_sorted.bw', paste0(RESULTBASEDIR, '/Blood_Sample_Reference_Peaks/CD38_Neg_NB/merged_ChIP_H3K27me3.bigWig'), paste0(RESULTBASEDIR, '/Blood_Sample_Reference_Peaks/CD38_Neg_NB/merged_ChIP_H3K4me3.bigWig'), paste0(RESULTBASEDIR, '/Blood_Sample_Reference_Peaks/CD38_Neg_NB/merged_ChIP_H3K9me3.bigWig'), paste0(RESULTBASEDIR, '/Blood_Sample_Reference_Peaks/CD38_Neg_NB/merged_ChIP_H3K36me3.bigWig'), paste0(RESULTBASEDIR, '/Blood_Sample_Reference_Peaks/CD38_Neg_NB/merged_ChIP_H3K4me1.bigWig'))

# NK ChIP-seq bigwig files for different proteins or histone modifications of interest

# comment - sourya
# when only one replicate and corresponding bigwig files were considered
# NK_ChIP_BigWig_FileList <- c('/mnt/BioAdHoc/Groups/vd-vijay/sourya/genomes/REF_ChIP_Seq_Peaks/LJI_Custom/NK/Merge_Specific_Donors/merged_R24_NK_ChIPSeq.bw', paste0(RESULTBASEDIR, '/Blood_Sample_Reference_Peaks/CD56_NK/24984.Blueprint.ERS158628.H3K27me3.signal.bigWig'), paste0(RESULTBASEDIR, '/Blood_Sample_Reference_Peaks/CD56_NK/25226.Blueprint.ERS257333.H3K4me3.signal.bigWig'), paste0(RESULTBASEDIR, '/Blood_Sample_Reference_Peaks/CD56_NK/24985.Blueprint.ERS158628.H3K9me3.signal.bigWig'), paste0(RESULTBASEDIR, '/Blood_Sample_Reference_Peaks/CD56_NK/25038.Blueprint.ERS206510.H3K36me3.signal.bigWig'), paste0(RESULTBASEDIR, '/Blood_Sample_Reference_Peaks/CD56_NK/25039.Blueprint.ERS206510.H3K4me1.signal.bigWig'))

NK_ChIP_BigWig_FileList <- c('/mnt/BioAdHoc/Groups/vd-vijay/sourya/genomes/REF_ChIP_Seq_Peaks/LJI_Custom/NK/Merge_Specific_Donors/merged_R24_NK_ChIPSeq.bw', paste0(RESULTBASEDIR, '/Blood_Sample_Reference_Peaks/CD56_NK/merged_ChIP_H3K27me3.bigWig'), paste0(RESULTBASEDIR, '/Blood_Sample_Reference_Peaks/CD56_NK/merged_ChIP_H3K4me3.bigWig'), paste0(RESULTBASEDIR, '/Blood_Sample_Reference_Peaks/CD56_NK/merged_ChIP_H3K9me3.bigWig'), paste0(RESULTBASEDIR, '/Blood_Sample_Reference_Peaks/CD56_NK/merged_ChIP_H3K36me3.bigWig'), paste0(RESULTBASEDIR, '/Blood_Sample_Reference_Peaks/CD56_NK/merged_ChIP_H3K4me1.bigWig'))

