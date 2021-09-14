# MHC-alleles-specificity-analysis
Studying correlation between MHC alleles specificity and protein similarities during internship at Huber group (EMBL).

- HLA_A_alignment_analysis.R is the script this whole analysis
- data.tar.gz consists of
  - A_prot.msf is alignment of HLA-A alleles downloaded from IMGT database
  - allele_specif_cor_matr.csv - EL-scores correlations matrix
  - allele_specif_cor_matr_BA.csv - BA-scores correlations matrix
  - similarity_hlaa_wo_dupl.csv - an example of similarities estimated from the alignment
  - script_mhcI.sh - script for netMHCpan 4.1 running for all MHCI alleles
