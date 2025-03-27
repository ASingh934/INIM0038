##############################################################################
# Perform PCA
##############################################################################
pca_result <- prcomp(gene_data, center = TRUE, scale. = FALSE)

###############################################################################
# Summary of PCA - prints the SD & proportion of variance of each PC
##############################################################################
print(summary(pca_result)) 

###############################################################################
# Express PC1 as a linear combination of its dimensions (top 50 weighted genes)
################################################################################
pc1_loadings <- pca_result$rotation[,1]
sorted_indices <- order(abs(pc1_loadings), decreasing = TRUE)
top_50_indices <- sorted_indices[1:50]
top_50_loadings <- pc1_loadings[top_50_indices]
pc1_equation <- paste0(round(top_50_loadings, 4), " * ", names(top_50_loadings), 
                       collapse = " + ")

cat("\nPC1 is expressed as:\n")
cat("PC1 =", pc1_equation, "\n")

###############################################################################
# First PCA plot with patients grouped according to whether they are BioAID or Control
# All patients are plotted here
################################################################################

# Visualise PCA using a PCA plot with only 2 colors
plot_data <- data.frame(group_labels = group_labels)
rownames(plot_data) <- rownames(gene_data)

pca_plot <- autoplot(
  pca_result,
  data = plot_data,
  colour = "group_labels",
  size = 3,
  shape = 16
) +
  scale_color_manual(values = c(
    "BioAID"   = "red", 
    "Controls" = "blue"
  )) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  labs(title = "PCA of Gene Expression Data (BioAID vs Controls)", color = "Group") +
  theme_pubr()

# Print the updated PCA plot
print(pca_plot)

###############################################################################
# Second PCA plot with patients grouped according to clinical labels
# IN this plot, patients whose clinical data aren't available are not plotted
###############################################################################

# --- Step 1: Prepare plotting information ---
# Merge gene_ids with clinical_match to get micro_diagnosis info
plot_info <- gene_ids %>%
  filter(original_id %in% matched_patient_ids) %>%
  left_join(clinical_match, by = c("cleaned_id" = "matched_id"))

# Ensure the row names of plot_info match those in gene_data
rownames(plot_info) <- plot_info$original_id

# --- Step 2: Subset PCA results to matched patients ---
# Assume pca_result was computed on the full gene_data earlier:
# pca_result <- prcomp(gene_data, center = TRUE, scale. = FALSE)

# Subset the PCA scores to only include matched patients
pca_result_matched <- pca_result
pca_result_matched$x <- pca_result$x[matched_patient_ids, ]

# --- Step 3: Create the PCA plot using autoplot ---
# Define color mapping based on micro_diagnosis:
color_mapping <- c("None" = "blue", 
                   "Bacterial" = "red", 
                   "Viral" = "green", 
                   "Bacterial & Viral" = "black")

pca_plot_matched <- autoplot(pca_result_matched, data = plot_info,
                             colour = "micro_diagnosis", size = 3, shape = 16) +
  scale_color_manual(values = color_mapping) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  labs(title = "PCA Plot for Matched Patients by Clinical micro_diagnosis",
       color = "micro_diagnosis") +
  theme_pubr()

# Display the PCA plot
print(pca_plot_matched)
################################################################################
# Produce a scree plot 
###############################################################################

# Scree plot to visualise the variance captured by the first 50 PCs
explained_variance <- pca_result$sdev^2 / sum(pca_result$sdev^2)
scree_data <- data.frame(PC = 1:
                           length(explained_variance), Variance = explained_variance)
scree_data_50 <- scree_data[1:50, ]

scree_plot <- ggplot(scree_data_50, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  labs(title = "Scree Plot (First 50 PCs)", x = "Principal Components", y = "Variance Explained") +
  theme_minimal()

print(scree_plot)
