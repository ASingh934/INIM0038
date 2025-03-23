**DataPrep.R**
This R script contains the code which prepares the gene expression data for processing. The gene_data is transposed to ensure that the row names match the patient IDs, and the column names match the gene names. The non-protein coding genes are filtered out with a CSV file containing the gene names and their corresponding biotypes. Patient grouping is then performed - there are 4 categories of patients: Oxford, UCL, UHB and Controls.

**PCA.R**
PCA is performed on gene_data. 
