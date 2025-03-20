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

# Calculate percentage variance explained for PC1 and PC2
explained_variance <- pca_result$sdev^2 / sum(pca_result$sdev^2)
pc1_percent <- round(explained_variance[1] * 100, 2)
pc2_percent <- round(explained_variance[2] * 100, 2)

# Create the PCA plot using only patients with clinical data,
# coloring points based on their micro_diagnosis and updating axis labels accordingly
p <- ggplot(pca_data_clinical, aes(x = PC1, y = PC2, color = micro_diagnosis)) +
  geom_point(size = 3, shape = 16) +
  scale_color_manual(
    name = "Micro Diagnosis",
    values = c(
      "Bacterial"         = "red",
      "Viral"             = "blue",
      "None"              = "green",
      "Bacterial & Viral" = "purple"
    )
  ) +
  labs(
    title = "PCA of Gene Expression Data with Clinical Grouping",
    x = paste0("PC1 (", pc1_percent, "%)"),
    y = paste0("PC2 (", pc2_percent, "%)")
  ) +
  theme_pubr()

# Print the updated PCA plot with clinical grouping
print(p)

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
