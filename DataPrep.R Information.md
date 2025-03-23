**DataPrep.R**

This R script contains the code which prepares the gene expression data for processing. The _**gene_data**_ is transposed to ensure that the row names match the patient IDs, and the column names match the gene names. The non-protein coding genes are filtered out with a CSV file containing the gene names and their corresponding biotypes. Patient grouping is then performed - there are 4 categories of patients: Oxford, UCL, UHB and Controls. The code merges the first 3 groups into a singular BioAID group.

It is assumed that all BioAID patients have ARIs, whilst all controls are healthy. 
The clinical_data contains the final confirmatory diagnosis of all Oxford and UCL patients (not the Controls or UHB patients)
