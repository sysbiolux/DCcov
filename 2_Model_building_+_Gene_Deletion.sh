## Extract the metadata from the FPKM and read count files

# The time-series study
# Extract metadata and\nDividing the lung FPKM according to ecah phenotype to a seperate csv file for model reconstruction 
jupyter nbconvert --to notebook --inplace --execute FPKM_Dividing_Timeseries_Study.ipynb

# Extract gene lengths for RPKM for the severity study
Rscript Extract_gene_lengths.R

# Extract metadata and calculate RPKM for model building for the severity study
jupyter nbconvert --to notebook --inplace --execute  Calculate_RPKM.ipynb


mkdir models
mkdir models/timeseries models/severity
mkdir KO_data
mkdir KO_data/timeseries KO_data/severity

# Model building and single and double gene deletion
## on the time series study
matlab -nodisplay -nodesktop -r Model_building_timeseries_Recon3D.m
matlab -nodisplay -nodesktop -r Model_building_timeseries_Recon2.m

## Model building and single and double gene deletion on the time series study
matlab -nodisplay -nodesktop -r Model_building_severity_Recon3D.m
matlab -nodisplay -nodesktop -r Model_building_severity_Recon2.m

# Summerize all KO results from different conditions and models
matlab -nodisplay -nodesktop -r Summarize_KO_Results_severity.m
matlab -nodisplay -nodesktop -r Summarize_KO_Results_timeseries.m

# Pathway analysis of the SKO and DKO using KEGG
Enrichment_Analysis_KEGG.ipynb

# Metabolic Pathway Analysis of the SKO genes
matlab -nodisplay -nodesktop -r  SKO_Metabolic_Pathway_Analysis.m
jupyter nbconvert --to notebook --inplace --execute Visualize_SKO_Pathways.ipynb

## Drug repurpusing ##
mkdir drugbank

## download DrukBank .xml  manually

## convert .xml to csv and integrate MedDRA side effect database
jupyter nbconvert --to notebook --inplace --execute  DrugBank-MedDRA_Integration.ipynb
# Calculate essentiality and safety and apply drug repurusing using DrugBank 
jupyter nbconvert --to notebook --inplace --execute  Drug_Repurposing_DrugBank.ipynb
