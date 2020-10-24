#!/bin/bash

##=================
# script to download reference genome hg19 (GTF files and fasta files)
# gencode v19
##=================

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz

gunzip gencode.v19.annotation.gtf.gz
gunzip GRCh37.p13.genome.fa.gz

mv GRCh37.p13.genome.fa hg19.fa

## following fields are extracted
## Chr     Start   End     Strand  GeneID  GeneType        GeneName        TranscriptID
awk -F'[;\t ]' '{if (substr($1,1,1)!="#")  {print $1"\t"$4"\t"$5"\t"$7"\t"$10"\t"$16"\t"$22"\t"$13}}' gencode.v19.annotation.gtf > gEncode_Genes_Complete_type_GeneID_TranscriptID_NEW.gtf

## following fields are extracted
## Chr     Start   End     TSSIdx  Strand  GeneID  GeneType        GeneName
awk '{if ($4=="+") {print $1"\t"$2"\t"$2"\t"$4"\t"$5"\t"$6"\t"$7} else {print $1"\t"$3"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}}' gEncode_Genes_Complete_type_GeneID_TranscriptID_NEW.gtf | awk -F"\t" '!seen[$1, $4, $5, $6, $7]++' - | awk '{print $1"\t"$2"\t"$3"\t"(NR-1)"\t"$4"\t"$5"\t"$6"\t"$7}' - > gEncode_Genes_Complete_NEW_TSS.gtf


