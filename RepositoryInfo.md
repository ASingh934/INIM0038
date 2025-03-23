**DataPrep.R**

This R script contains the code which prepares the gene expression data for processing. The **gene_data** is transposed to ensure that the row names match the patient IDs, and the column names match the gene names. The non-protein coding genes are filtered out with a CSV file containing the gene names and their corresponding biotypes. Patient grouping is then performed - there are 4 categories of patients: Oxford, UCL, UHB and Controls. The code merges the first 3 groups into a singular BioAID group.

It is assumed that all BioAID patients have ARIs, whilst all controls are healthy. 
The clinical_data contains the final confirmatory diagnosis of all Oxford and UCL patients (not the Controls or UHB patients)

**PCA.R**

PCA is performed on **gene_data**. For each PC, the Standard Deviation and proportion of total variance is displayed. For PC1, the top 50 most weighted genes are displayed to show which genes are most effective at clustering the patients. Then, 2 PCA plots are drawn.

**PCA plot 1**: All patients are plotted here, grouped according to whether they are BioAID or controls. We can see from visual inspection that the controls cluster in the Right.

![image](https://github.com/user-attachments/assets/314b2404-061b-4ead-b66f-267eedd52f9c)

**PCA plot 2**: For the Oxford & UCL patients whose clinical labels are available, they are replotted on the same PCA plot, and grouped by whether or not they have an ARI (and if so, then the pathogen which is causing it). From this, we can see no clear clustering of any group in particular, and we conclude that PCA plots are not an effective method of endotyping ARI patients.

![image](https://github.com/user-attachments/assets/11d5d1af-ee41-4d3b-a153-d01b4f8f8de5)

A scree plot is shown in the code for completeness.

