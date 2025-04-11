optimum_cluster <- read.csv("C:/Users/ABSin/Downloads/4. OptimumLeidenCluster.csv")


###############################################################################

# Sensitivity & Specificity

##############################################################################

# Step 1: Identify the community with the most viral patients
viral_counts <- optimum_cluster %>%
  filter(micro_diagnosis == "Viral") %>%
  count(community, name = "viral_count")

dominant_viral_comm <- viral_counts %>%
  arrange(desc(viral_count)) %>%
  slice(1) %>%
  pull(community)

cat("Dominant Viral Community:", dominant_viral_comm, "\n")

# Step 2: Create predicted and true labels
optimum_cluster <- optimum_cluster %>%
  mutate(
    predicted_label = if_else(community == dominant_viral_comm, "Viral", "Non-Viral"),
    true_label = if_else(micro_diagnosis == "Viral", "Viral", "Non-Viral")
  )

# Step 3: Build confusion matrix
conf_matrix <- table(Predicted = optimum_cluster$predicted_label,
                     Actual = optimum_cluster$true_label)

print(conf_matrix)

# Step 4: Calculate Sensitivity and Specificity
TP <- conf_matrix["Viral", "Viral"]
FP <- conf_matrix["Viral", "Non-Viral"]
TN <- conf_matrix["Non-Viral", "Non-Viral"]
FN <- conf_matrix["Non-Viral", "Viral"]

sensitivity <- TP / (TP + FN)
specificity <- TN / (TN + FP)

cat("Sensitivity:", round(sensitivity, 3), "\n")
cat("Specificity:", round(specificity, 3), "\n")

################################################################################

# Limma 

################################################################################

# --- 1. Load Data ---
# Assuming 'gene_data' (patients as rows, genes as columns, numeric) is still in your environment
# Load the clustering results
optimum_cluster <- read.csv("C:/Users/ABSin/Downloads/4. OptimumLeidenCluster.csv")

# --- 2. Match Samples between Gene Data and Cluster Data ---

# Create cleaned IDs for gene_data row names (reuse logic from previous steps)
gene_ids_df <- data.frame(original_id = rownames(gene_data), stringsAsFactors = FALSE) %>%
  mutate(cleaned_id = case_when(
    grepl("^(OX|ox)", original_id) ~ gsub("[^A-Za-z0-9]", "", original_id),
    original_id == "UP.3447" ~ "UP3447", # Specific case if needed based on your data
    TRUE ~ original_id
  ))

# Create cleaned IDs for the cluster data patient_id
optimum_cluster <- optimum_cluster %>%
  mutate(cleaned_id = case_when(
    grepl("^(OX|ox)", patient_id) ~ gsub("[^A-Za-z0-9]", "", patient_id),
    patient_id == "UP3447" ~ "UP3447", # Match potential clinical data format
    TRUE ~ patient_id
  ))

# Find common cleaned IDs
common_cleaned_ids <- intersect(gene_ids_df$cleaned_id, optimum_cluster$cleaned_id)
cat("Found", length(common_cleaned_ids), "common samples between gene expression data and cluster data.\n")

# Filter cluster data to common samples and get original gene_data IDs
cluster_info_matched <- optimum_cluster %>%
  filter(cleaned_id %in% common_cleaned_ids) %>%
  left_join(gene_ids_df, by = "cleaned_id") %>%
  select(original_id, community, micro_diagnosis, cleaned_id) # Keep original_id for filtering gene_data

# Filter gene_data to common samples
# Ensure the order matches the cluster_info_matched dataframe
gene_data_matched_cluster <- gene_data[cluster_info_matched$original_id, ]

# --- 3. Prepare Data for Limma ---

# Transpose gene data: genes as rows, samples as columns
expr_matrix <- t(gene_data_matched_cluster)
# Ensure it's numeric
expr_matrix <- apply(expr_matrix, 2, as.numeric)
rownames(expr_matrix) <- colnames(gene_data_matched_cluster) # Gene symbols
colnames(expr_matrix) <- rownames(gene_data_matched_cluster) # Original sample IDs

# Create the design matrix based on 'community'
# Ensure community is treated as a factor
cluster_info_matched$community <- factor(paste0("C", cluster_info_matched$community)) # Prepend "C" for clarity
design <- model.matrix(~0 + community, data = cluster_info_matched)
colnames(design) <- levels(cluster_info_matched$community) # Clean column names

# Check dimensions
print(dim(expr_matrix))
print(dim(design))
if (ncol(expr_matrix) != nrow(design)) {
  stop("Mismatch between expression matrix columns and design matrix rows.")
}

# --- 4. Define Contrasts ---
cont.matrix <- makeContrasts(
  C1vsRest = C1 - (C2 + C3 + C4) / 3,
  C2vsRest = C2 - (C1 + C3 + C4) / 3,
  C3vsRest = C3 - (C1 + C2 + C4) / 3,
  C4vsRest = C4 - (C1 + C2 + C3) / 3,
  levels = design
)

print("Contrast Matrix:")
print(cont.matrix)

# --- 5. Run Limma ---
fit <- lmFit(expr_matrix, design)
fit_contrasts <- contrasts.fit(fit, contrasts = cont.matrix)
fit_ebayes <- eBayes(fit_contrasts)

# --- 6. Extract DEGs for each contrast (Revised Storage) ---

# Define thresholds
adj_p_threshold <- 0.05
logfc_threshold <- 1.0

# Initialize the main list to store all details
deg_details <- list()
# Initialize list to store summary counts for the plot
deg_summary <- list()

# Loop through each contrast
contrast_names <- colnames(cont.matrix)

for (contrast in contrast_names) {
  cat("\nProcessing contrast:", contrast, "\n")
  
  # Extract results using topTable
  res_table <- topTable(fit_ebayes, coef = contrast, adjust.method = "BH", number = Inf)
  
  # Filter based on adjusted p-value and logFC
  up_genes <- res_table %>%
    rownames_to_column("Gene") %>%
    filter(adj.P.Val < adj_p_threshold & logFC > logfc_threshold) %>%
    pull(Gene)
  
  down_genes <- res_table %>%
    rownames_to_column("Gene") %>%
    filter(adj.P.Val < adj_p_threshold & logFC < -logfc_threshold) %>%
    pull(Gene)
  
  # Store the gene lists in the nested structure
  deg_details[[contrast]] <- list(
    upregulated = up_genes,
    downregulated = down_genes
  )
  
  # Store the counts for the summary plot
  deg_summary[[contrast]] <- data.frame(
    Contrast = contrast,
    Upregulated = length(up_genes),
    Downregulated = length(down_genes)
  )
  
  cat("  Upregulated:", length(up_genes), "\n")
  cat("  Downregulated:", length(down_genes), "\n")
}

# Combine summary counts into a single data frame
summary_df <- bind_rows(deg_summary)

print("\nSummary of DEGs:")
print(summary_df)

# --- 7. Visualize DEG Counts ---

# Reshape data for ggplot
plot_data <- summary_df %>%
  pivot_longer(cols = c("Upregulated", "Downregulated"), names_to = "Direction", values_to = "Count") %>%
  mutate(Direction = factor(Direction, levels = c("Upregulated", "Downregulated")))

# Create the bar plot
deg_barplot <- ggplot(plot_data, aes(x = Contrast, y = Count, fill = Direction)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = Count),
            position = position_dodge(width = 0.8),
            vjust = -0.3,
            size = 3) +
  scale_fill_manual(values = c("Upregulated" = "firebrick", "Downregulated" = "steelblue")) +
  labs(
    title = "Differentially Expressed Genes per Contrast",
    subtitle = paste("Adjusted P-value <", adj_p_threshold, " & |LogFC| >", logfc_threshold),
    x = "Contrast Comparison",
    y = "Number of Genes",
    fill = "Regulation"
  ) +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print the plot
print(deg_barplot)

############################### Lists of DEGs ##################################
print(deg_details$C1vsRest$upregulated)
# deg_details$C1vsRest$upregulated - Genes upregulated in C1 compared to C2, C3, and C4
# deg_details$C1vsRest$downregulated - Genes downregulated in C1 compared to the rest

# deg_details$C2vsRest$upregulated - Genes upregulated in C2 compared to C1, C3, and C4
# deg_details$C2vsRest$downregulated - Genes downregulated in C2 compared to the rest

# deg_details$C3vsRest$upregulated - Genes upregulated in C3 compared to C1, C2, and C4
# deg_details$C3vsRest$downregulated - Genes downregulated in C3 compared to the rest

# deg_details$C4vsRest$upregulated - Genes upregulated in C4 compared to C1, C2, and C3
# deg_details$C4vsRest$downregulated - Genes downregulated in C4 compared to the rest

###############################################################################

# XGR

###############################################################################

library(viridis)

# Save each DEG list as a separate .rds file
for (contrast in names(deg_details)) {
  for (direction in c("upregulated", "downregulated")) {
    gene_vector <- deg_details[[contrast]][[direction]]
    
    if (length(gene_vector) == 0) next
    
    file_name <- paste0(contrast, "_", direction, ".rds")
    saveRDS(gene_vector, file = file_name)
    cat("Saved:", file_name, "\n")
  }
}

# Helper: Read RDS → Enrich → Plot (top 10 terms max)
plot_from_rds <- function(file_path) {
  gene_vector <- readRDS(file_path)
  parts <- strsplit(basename(file_path), "_|\\.rds")[[1]]
  contrast_name <- parts[1]
  direction <- parts[2]
  
  cat("Processing:", contrast_name, "-", direction, "\n")
  
  if (length(gene_vector) == 0) {
    cat("No genes in file:", file_path, "\n")
    return(NULL)
  }
  
  # Enrichment
  eTerm <- xEnricherGenes(data = gene_vector, ontology = "MsigdbC2REACTOME", background = "protein coding")
  res <- xEnrichViewer(eTerm, sortBy = "adjp", details = TRUE)
  
  # Take top 10 only (if they exist)
  res <- res %>% arrange(adjp) %>% slice_head(n = 10)
  
  if (nrow(res) == 0 || all(is.na(res$fc))) {
    cat("No enrichment to show for", contrast_name, direction, "\n")
    return(NULL)
  }
  
  # Format adjusted p-values for labels
  res$formatted_adjp <- sapply(res$adjp, function(p) {
    sci <- format(p, scientific = TRUE)
    parts <- strsplit(sci, "e")[[1]]
    base <- as.numeric(parts[1])
    exponent <- as.integer(parts[2])
    bquote(.(base) %*% 10^.(exponent))
  })
  
  # Plot top 10
  p <- ggplot(res, aes(x = fc, y = reorder(name, fc))) +
    geom_point(aes(size = -log10(adjp), color = -log10(adjp))) +
    geom_text(aes(label = formatted_adjp), hjust = -0.2, size = 3, parse = TRUE) +
    scale_color_viridis_c() +
    guides(size = "none") +
    labs(
      title = paste(contrast_name, "-", direction, "Genes (Top 10 Pathways)"),
      x = "Fold Enrichment",
      y = NULL,
      color = expression(-log[10]~"adj p-value")
    ) +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(size = 10),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.margin = margin(10, 20, 10, 10)
    ) +
    coord_cartesian(xlim = c(0, max(res$fc, na.rm = TRUE) + 2))
  
  return(p)
}

# Get list of saved .rds DEG files
deg_files <- list.files(pattern = "C[1-4]vsRest_(upregulated|downregulated)\\.rds$")

# Loop and plot each (top 10 terms only)
for (file in deg_files) {
  plot_obj <- plot_from_rds(file)
  if (!is.null(plot_obj)) print(plot_obj)
}
