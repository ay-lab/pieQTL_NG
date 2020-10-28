#!/usr/bin/env Rscript

source('Header.r')
source('Functions.r')

##====================
# this script extracts the promoter interacting SNPs from the direct and indirect connections
##====================

###====== input parameters - users need to edit them

# promoter interacting SNP file - direct connection
PromIntSNP_DirectFile <- 'Promoters_Direct_PorE_Indirect_EQTL_PromSlack_5Kb_SNPSlack_5Kb_TSS_Included.bed'

# promoter interacting SNP file - indirect connection
PromIntSNP_InDirectFile <- 'Promoters_Direct_PorE_Indirect_EQTL_PromSlack_5Kb_SNPSlack_5Kb_TSS_Included.bed'

# output directory to store the results
curroutdir <- getwd()

##===== loop for all chromosomes
GeneID_SNP_File <- paste0(curroutdir, '/Input_Gene_SNP_pairs.txt')
if (file.exists(GeneID_SNP_File) == TRUE) {
	system(paste("rm", GeneID_SNP_File)
}

for (chridx in 1:length(ChrList_NameNum)) {
	chrname <- ChrList_NameNum[chridx]
	GeneID_SNP_File_Currchr <- paste0(curroutdir, '/Input_Gene_SNP_pairs_', chrname, '.txt')
	
	##=== extract the chromosomes, gene ID, and SNP position from the direct and indirect interactions
	system(paste0("awk -F\'[\t]\' \'{if ((NR>1) && ($1==\"", chrname, "\")) {print $1\"\t\"$6\"\t\"$16}}\' ", PromIntSNP_DirectFile, " > ", GeneID_SNP_File_Currchr))
	system(paste0("awk -F\'[\t]\' \'{if ((NR>1) && ($1==\"", chrname, "\")) {print $1\"\t\"$6\"\t\"$19}}\' ", PromIntSNP_InDirectFile, " >> ", GeneID_SNP_File_Currchr))

	system(paste0("sort -k2,2 -k4,4n ", GeneID_SNP_File_Currchr, " | uniq >> ", GeneID_SNP_File))

}	# end chromosome loop

