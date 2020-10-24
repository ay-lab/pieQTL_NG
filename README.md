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





