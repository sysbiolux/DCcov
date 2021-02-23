### Goal: Download the raw expression data and calculate the differentially exxoressed mrtabolic pathways

# Create a new directory to downlaoid the read count data in
mkdir data
mkdir data/time_series_study
mkdir data/severity_study

## Begin downloading the FPKM / RPKM and metadata

# The time-series study
# Download manually the FPKM for the polyA of GSE148729: https://filetransfer.mdc-berlin.de/?u=CVXckugR&p=MACT6Xw9

# SupplementaryData2_Calu3_polyA_series1_fpkm.tsv	
# SupplementaryData5_H1299_polyA_fpkm.tsv
# SupplementaryData3_Calu3_polyA_series2_fpkm.tsv

# The severity study
#wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147507/suppl/GSE147507_RawReadCounts_Human.tsv.gz -P data/severity_study/
#gunzip data/severity_study/GSE147507_RawReadCounts_Human.tsv.gz

# PCA analysis on the severity
mkdir Figs
#Rscript PCA_analysis.R
# DEG analysis
mkdir ./data/severity_study/DEGs
Rscript DEG_analysis.R

# Calculate differentially expressed reactions (DER) for the severity study
matlab -nodisplay -nodesktop -r Differentially_Expressed_Reactions_Table.m
# Visualize DER\nouput in Figs
jupyter nbconvert --to notebook --inplace --execute Visualize_Differentially_Expressed_Pathways.ipynb


