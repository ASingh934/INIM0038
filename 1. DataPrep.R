library(ggbiplot)
library(dplyr)
library(ggrepel)
library(ggpubr)
library(ggfortify)
library(tidyr)
library(Matrix)
library(igraph)
library(ggplot2)
library(gridExtra)
library(igraph)
library(dplyr)
library(pheatmap)
library(scales)
library(stringr)
library(tibble)

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

################################################################################
# Sort patients into groups: Merge Oxford, UCL, and UHB into BioAID
##############################################################################
sample_ids <- rownames(gene_data)
group_labels <- case_when(
  grepl("^(OX|ox)", sample_ids) ~ "BioAID",
  grepl("^UP", sample_ids) ~ "BioAID",
  grepl("^(667|X)", sample_ids) ~ "BioAID",
  grepl("^WH0", sample_ids) ~ "Controls"
)

# Convert to factor with correct order
group_labels <- factor(group_labels, levels = c("BioAID", "Controls"))

###############################################################################
# Determine which patients in gene_data have clinical labels i.e. are in clinical_data
#################################################################################

# Create a new column in clinical_data: matched_id
clinical_data <- clinical_data %>%
  mutate(matched_id = case_when(
    # For UP entries: use the UIN directly, except UP3447 becomes "UP3447" 
    grepl("^UP", UIN) ~ if_else(UIN == "UP3447", "UP3447", UIN),
    # For OX entries (case-insensitive): remove special characters so "OX-10411" becomes "OX10411"
    grepl("^(OX|ox)", UIN) ~ gsub("[^A-Za-z0-9]", "", UIN),
    TRUE ~ NA_character_
  ))

# --- Step 2: Process gene_data row names ---
# Create a data frame of gene_data row names
gene_ids <- data.frame(original_id = rownames(gene_data), stringsAsFactors = FALSE)

# Clean gene_data row names similar to clinical_data
gene_ids <- gene_ids %>%
  mutate(cleaned_id = case_when(
    # For OX entries, remove any non-alphanumeric characters (e.g. "OX.10411" becomes "OX10411")
    grepl("^(OX|ox)", original_id) ~ gsub("[^A-Za-z0-9]", "", original_id),
    # For UP entries: if the row name is "UP.3447", remove the dot to match clinical "UP3447"
    original_id == "UP.3447" ~ "UP3447",
    # Otherwise, use the row name as is (e.g., "UP3071")
    TRUE ~ original_id
  ))

# --- Step 3: Identify matching patients ---
# Filter clinical_data to keep only rows with valid matched_id (UP or OX type)
clinical_match <- clinical_data %>% filter(!is.na(matched_id))

# Find the intersection of cleaned gene IDs and clinical matched IDs
common_ids <- intersect(gene_ids$cleaned_id, clinical_match$matched_id)

# Get the original gene_data row names corresponding to these common IDs
matched_patient_ids <- gene_ids$original_id[gene_ids$cleaned_id %in% common_ids]

# Subset gene_data for matched patients
matched_patients <- gene_data[matched_patient_ids, ]
cat("matched_patients dataframe has been created with", nrow(matched_patients), "patients.\n")

###########################################################################################
# matched_patient_ids contains all the patient IDs (row names) in gene_data who are in clinical_data
#############################################################################################
