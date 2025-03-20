library(ggbiplot)
library(dplyr)
library(ggrepel)
library(ggpubr)
library(ggfortify)
library(tidyr)
library(igraph)
library(ggplot2)
library(gridExtra)
library(igraph)
library(dplyr)
library(pheatmap)
library(scales)

gene_data <- read.csv("C:/Users/ABSin/Downloads/BIOAID_tpm_PC0.001_log2_genesymbol_dedup (1).csv")
clinical_data <- readRDS("C:/Users/ABSin/Downloads/BioAID_resp_master_v2_March25.rds")

################################################################################

# Transpose data so patients are observations (rows) and genes are variables (cols)
gene_data <- as.data.frame(t(gene_data))

###############################################################################

# Update the column headings to be the gene names
new_headers <- gene_data[1,] 
gene_data <- gene_data[-1,] 
colnames(gene_data) <- new_headers 

# Convert columns to numeric properly (ensuring negative values are kept)
gene_data[] <- lapply(gene_data, function(x) as.numeric(as.character(x)))

##############################################################################

# Filter out non-protein coding genes

# Import the csv file containing gene symbols and their biotype
annotated_data <- read.csv("C:/Users/ABSin/Downloads/Annotatedexample_tpm_PC0.001_log2.csv")

# Filter to retain only protein-coding genes
valid_genes <- annotated_data %>%
  filter(gene_biotype == "protein_coding") %>%
  select(external_gene_name) %>%
  pull()

# Keep only columns in gene_data that match the valid gene names
gene_data <- gene_data %>% select(all_of(intersect(names(gene_data), valid_genes)))
