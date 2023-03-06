pieQTL_NG
-----------

Repository storing custom scripts to derive promoter interacting eQTLs (pieQTLs) from HiChIP data

Manuscript
-----------

Promoter-interacting expression quantitative trait loci are enriched for functional genetic variants

(accepted in Nature Genetics)

Authors: Vivek Chandra, Sourya Bhattacharyya, Benjamin Schmiedel, Ariel Madrigal, Cristian Gonzalez-Colin, Stephanie Fotsing, Austin Crinklaw, Greggory Seumois, Dr. Pejman Mohammadi, Dr. Mitchell Kronenberg, Dr. Bjoern Peters, Dr. Ferhat Ay and Dr. Pandurangan Vijayanand

La Jolla Institute For Immunology
La Jolla, CA 92037, USA

Prerequisites
--------------

(1) Download reference genome (hg19) annotation and fasta files using the "script_hg19_GTF.sh" script.

(2) Install the packages HOMER, deepTools and mention their paths in lines 20-23 of the script "Header.r"

(3) Install the following R (>=3.4.3) packages: GenomicRanges, RColorBrewer, gplots, ggplot2, edgeR, stats, factoextra, pheatmap, dplyr, parallel, plotly, igraph, data.table, ggpubr

(4) Install HiC-pro, take a look at the directory "HiC_pro_sample" to check the sample HiC-pro configuration script.

(5) Install FitHiChIP and check its documentation to run on HiChIP data. 

(6) Install bedtools (>= 2.26.0)


Modules
---------

(In each of these source codes, user needs to edit to put their custom files)

** HiChIP_Reproducibility.r

	Correlation of a pair of HiChIP contact maps, computed per chromosome.

** PCA.r	

	PCA of donors of different cell types, according to the FitHiChIP loops.

** pieQTL.r	

	Script to compute the pieQTLs. The output file "dumped_pieQTL.txt" contains the pieQTLs (direct + indirect).

** ultralong_pieQTL_Promoter_Interacting_SNP.r

	Script to extract the gene and SNP pairs to be used for ultra-long pieQTL analysis.

	1) First apply the script pieQTL.r but using the complete list of SNPs, instead of the eQTLs. The first and second column of the SNP file should contain the chromosome name, and the SNP coordinate.

	2) Suppose, the direct interactions between promoters and SNPs are stored in the file "Promoters_Direct_PorE_Indirect_EQTL_PromSlack_5Kb_SNPSlack_5Kb_TSS_Included.bed". Further, the indirect interactions between promoters and SNPs are stored in the file "Promoters_Direct_PorE_Indirect_EQTL_PromSlack_5Kb_SNPSlack_5Kb_TSS_Included.bed". 

	These two files will be provided as inputs to this script.

	The output file "Input_Gene_SNP_pairs.txt" contains the genes and SNPs which participate in either direct and indirect interactions.

	These gene-SNP pairs are then to be put in the matrixQTL framework.


Other utility scripts
--------------------------

Conditional eQTL analysis:

Check the link: https://github.com/ay-lab/Conditional_eQTL


Contact
----------

For any queries, please email:

Sourya Bhattacharyya (sourya@lji.org)
