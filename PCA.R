library(ggplot2)
library(ggbiplot)
library(dplyr)
library(ggrepel)
library(ggpubr)
library(ggfortify)
library(tidyr)

gene_data <- read.csv("C:/Users/ABSin/Downloads/BIOAID_tpm_PC0.001_log2_genesymbol_dedup (1).csv")

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

###############################################################################

# Prepare gene_data for PCA
PCA_ready_data <- gene_data

###############################################################################

# Perform PCA
pca_result <- prcomp(PCA_ready_data, center = TRUE, scale. = FALSE)

###############################################################################

# Summary of PCA - prints the SD & proportion of variance of each PC
print(summary(pca_result)) 

###############################################################################

# Express PC1 as a linear combination of its dimensions (top 50 weighted genes)
pc1_loadings <- pca_result$rotation[,1]
sorted_indices <- order(abs(pc1_loadings), decreasing = TRUE)
top_50_indices <- sorted_indices[1:50]
top_50_loadings <- pc1_loadings[top_50_indices]
pc1_equation <- paste0(round(top_50_loadings, 4), " * ", names(top_50_loadings), collapse = " + ")

cat("\nPC1 is expressed as:\n")
cat("PC1 =", pc1_equation, "\n")

################################################################################

# Sort patients into groups based on sample IDs
sample_ids <- rownames(PCA_ready_data)
group_labels <- case_when(
  grepl("^(OX|ox)", sample_ids) ~ "Oxford",
  grepl("^UP", sample_ids) ~ "UCL",
  grepl("^(667|X)", sample_ids) ~ "UHB",
  grepl("^WH0", sample_ids) ~ "Controls"
)
group_labels <- factor(group_labels, levels = c("Oxford", "UCL", "UHB", "Controls"))

################################################################################

# Visualise PCA using a PCA plot
plot_data <- data.frame(group_labels = group_labels)
rownames(plot_data) <- rownames(PCA_ready_data)

pca_plot <- autoplot(
  pca_result,
  data = plot_data,
  colour = "group_labels",
  size = 1,
  shape = 16
) +
  scale_color_manual(values = c(
    "Oxford"   = "blue", 
    "UCL"      = "red", 
    "UHB"      = "green", 
    "Controls" = "purple"
  )) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  labs(title = "PCA of Gene Expression Data", color = "Group") +
  theme_pubr()

print(pca_plot)

###############################################################################

# Scree plot to visualise the variance captured by the first 50 PCs
explained_variance <- pca_result$sdev^2 / sum(pca_result$sdev^2)
scree_data <- data.frame(PC = 1:length(explained_variance), Variance = explained_variance)
scree_data_50 <- scree_data[1:50, ]

scree_plot <- ggplot(scree_data_50, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  labs(title = "Scree Plot (First 50 PCs)", x = "Principal Components", y = "Variance Explained") +
  theme_minimal()

print(scree_plot)

#############################################################################

# Plot distribution of MX1, IFI27 and IFI44L expression by group
genes_interest <- c("MX1", "IFI44L", "IFI27")
gene_subset <- PCA_ready_data[, genes_interest, drop = FALSE]
gene_subset$Group <- group_labels
gene_subset$SampleID <- rownames(PCA_ready_data)

gene_long <- gene_subset %>%
  pivot_longer(
    cols = all_of(genes_interest),
    names_to = "Gene",
    values_to = "Expression"
  )

gene_expression_plot <- ggplot(gene_long, aes(x = Group, y = Expression, color = Group)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  facet_wrap(~ Gene, scales = "free_y") +
  scale_color_manual(values = c(
    "Oxford"   = "blue", 
    "UCL"      = "red", 
    "UHB"      = "green", 
    "Controls" = "purple"
  )) +
  labs(
    title = "Expression Distribution of MX1, IFI44L, and IFI27",
    x = "Group",
    y = "Expression (log2 TPM)"   # Added units to y-axis label
  ) +
  theme_pubr() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(gene_expression_plot)

summary_stats <- gene_long %>%
  group_by(Gene, Group) %>%
  summarise(
    Count  = n(),
    Median = median(Expression, na.rm = TRUE),
    Q1     = quantile(Expression, 0.25, na.rm = TRUE),
    Q3     = quantile(Expression, 0.75, na.rm = TRUE),
    IQR    = IQR(Expression, na.rm = TRUE)
  ) %>%
  ungroup()

print(summary_stats)

#############################################################################

# Plot distribution of MX1, IFI27 and IFI44L expression by control/infected

# Collapse the four groups into two: "Controls" and "Intervention"
gene_long2 <- gene_long %>%
  mutate(Group2 = if_else(Group == "Controls", "Controls", "Intervention"))

gene_expression_plot_2 <- ggplot(gene_long2, aes(x = Group2, y = Expression, color = Group2)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  facet_wrap(~ Gene, scales = "free_y") +
  scale_color_manual(values = c(
    "Controls"     = "purple", 
    "Intervention" = "blue"
  )) +
  labs(
    title = "Expression Distribution of MX1, IFI44L, and IFI27 (2 Groups)",
    x = "Group",
    y = "Expression (log2 TPM)"   # Added units to y-axis label
  ) +
  theme_pubr() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(gene_expression_plot_2)

summary_stats2 <- gene_long2 %>%
  group_by(Gene, Group2) %>%
  summarise(
    Count  = n(),
    Median = median(Expression, na.rm = TRUE),
    Q1     = quantile(Expression, 0.25, na.rm = TRUE),
    Q3     = quantile(Expression, 0.75, na.rm = TRUE),
    IQR    = IQR(Expression, na.rm = TRUE)
  ) %>%
  ungroup()

print(summary_stats2)

###############################################################################

# Convert Expression (log2 TPM) into biomarker z scores using Controls as reference for each gene.
gene_long <- gene_long %>%
  group_by(Gene) %>%
  mutate(Expression_z = (Expression - mean(Expression[Group == "Controls"], na.rm = TRUE)) /
           sd(Expression[Group == "Controls"], na.rm = TRUE)) %>%
  ungroup()

# Create the plot using the z-scored values
gene_expression_plot <- ggplot(gene_long, aes(x = Group, y = Expression_z, color = Group)) +
  # Boxplot to display the overall distribution per group
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  # Jittered points to display individual sample values
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  # Create separate facets for each gene
  facet_wrap(~ Gene, scales = "free_y") +
  # Apply the same manual colour scheme as your PCA plot
  scale_color_manual(values = c(
    "Oxford"   = "blue", 
    "UCL"      = "red", 
    "UHB"      = "green", 
    "Controls" = "purple"
  )) +
  labs(
    title = "Expression Distribution of MX1, IFI44L, and IFI27 (Biomarker Z Scores)",
    x = "Group",
    y = "Biomarker Z Score"
  ) +
  theme_pubr() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Display the box plots
print(gene_expression_plot)

#############################################################################
# Plot distribution of MX1, IFI27 and IFI44L expression by control / infected

# Collapse the four groups into two: "Controls" and "Intervention"
# If the sample is labeled as "Controls", it remains "Controls"; otherwise, it becomes "Intervention".
gene_long2 <- gene_long %>%
  mutate(Group2 = if_else(Group == "Controls", "Controls", "Intervention"))

# Create the box plot for the two groups using the z-scored values
gene_expression_plot_2 <- ggplot(gene_long2, aes(x = Group2, y = Expression_z, color = Group2)) +
  # Boxplot to display the overall distribution per new group
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  # Jittered points to display individual sample values
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  # Facet by gene so each gene has its own panel
  facet_wrap(~ Gene, scales = "free_y") +
  # Apply a new manual colour scheme for the two groups
  scale_color_manual(values = c(
    "Controls"     = "purple", 
    "Intervention" = "blue"
  )) +
  labs(
    title = "Expression Distribution of MX1, IFI44L, and IFI27 (2 Groups; Biomarker Z Scores)",
    x = "Group",
    y = "Biomarker Z Score"
  ) +
  theme_pubr() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Display the new box plots
print(gene_expression_plot_2)
###############################################################################
# Create a new grouping variable (if not already done)
gene_long2 <- gene_long %>%
  mutate(Group2 = if_else(Group == "Controls", "Controls", "Intervention"))

# Histogram of Gene Expression Z-Scores by Group
histogram_plot_z <- ggplot(gene_long2, aes(x = Expression_z, fill = Group2)) +
  geom_histogram(alpha = 0.6, bins = 30, position = "identity") +
  facet_wrap(~ Gene, scales = "free") +
  labs(
    title = "Histogram of Gene Expression Z-Scores (2 Groups)",
    x = "Z-Score",
    y = "Count"
  ) +
  theme_pubr() +
  theme(legend.position = "right")

###############################################################################
# Controls with |z|> 1 have now been removed

# -------------------------------
# Step 1: Reassign groups and prepare for outlier removal & z-score recalculation
# -------------------------------
gene_long_recalc <- gene_long %>%
  # Create a new grouping variable: "Controls" remain, all others become "Intervention"
  mutate(Group2 = if_else(Group == "Controls", "Controls", "Intervention")) %>%
  group_by(Gene) %>%
  # Calculate temporary control-based statistics and temporary z-scores:
  mutate(
    temp_mean = mean(Expression[Group == "Controls"], na.rm = TRUE),
    temp_sd   = sd(Expression[Group == "Controls"], na.rm = TRUE),
    temp_z    = (Expression - temp_mean) / temp_sd
  ) %>%
  # Remove control samples that are outliers (|temp_z| > 1)
  filter(!(Group == "Controls" & abs(temp_z) > 1)) %>%
  # Recalculate control mean and SD using only the filtered controls
  mutate(
    new_control_mean = mean(Expression[Group == "Controls"], na.rm = TRUE),
    new_control_sd   = sd(Expression[Group == "Controls"], na.rm = TRUE)
  ) %>%
  # Recalculate the z-score for all samples based on the new control reference
  mutate(Expression_z = (Expression - new_control_mean) / new_control_sd) %>%
  ungroup()

# -------------------------------
# Step 2: Plot the results using the recalculated z-scores
# -------------------------------

# Box Plot of Recalculated Z-Scores by Group (Controls vs. Intervention)
box_plot_z <- ggplot(gene_long_recalc, aes(x = Group2, y = Expression_z, color = Group2)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  facet_wrap(~ Gene, scales = "free_y") +
  scale_color_manual(values = c("Controls" = "purple", "Intervention" = "blue")) +
  labs(
    title = "Expression Distribution (Biomarker Z Scores) after Removing Outlier Controls",
    x = "Group",
    y = "Biomarker Z Score"
  ) +
  theme_pubr() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
print(box_plot_z)

# Histogram of Recalculated Z-Scores by Group
histogram_plot_z <- ggplot(gene_long_recalc, aes(x = Expression_z, fill = Group2)) +
  geom_histogram(alpha = 0.6, bins = 30, position = "identity") +
  facet_wrap(~ Gene, scales = "free") +
  scale_fill_manual(values = c("Controls" = "purple", "Intervention" = "blue")) +
  labs(
    title = "Histogram of Gene Expression Z-Scores (Filtered Controls)",
    x = "Z-Score",
    y = "Count"
  ) +
  theme_pubr() +
  theme(legend.position = "right")
print(histogram_plot_z)

# Density Plot of Recalculated Z-Scores by Group
density_plot_z <- ggplot(gene_long_recalc, aes(x = Expression_z, fill = Group2, color = Group2)) +
  geom_density(alpha = 0.4) +
  facet_wrap(~ Gene, scales = "free") +
  scale_fill_manual(values = c("Controls" = "purple", "Intervention" = "blue")) +
  scale_color_manual(values = c("Controls" = "purple", "Intervention" = "blue")) +
  labs(
    title = "Density Plot of Gene Expression Z-Scores (Filtered Controls)",
    x = "Z-Score",
    y = "Density"
  ) +
  theme_pubr() +
  theme(legend.position = "right")
print(density_plot_z)

print(histogram_plot_z)

# Density Plot of Gene Expression Z-Scores by Group
density_plot_z <- ggplot(gene_long2, aes(x = Expression_z, fill = Group2, color = Group2)) +
  geom_density(alpha = 0.4) +
  facet_wrap(~ Gene, scales = "free") +
  labs(
    title = "Density Plot of Gene Expression Z-Scores (2 Groups)",
    x = "Z-Score",
    y = "Density"
  ) +
  theme_pubr() +
  theme(legend.position = "right")

print(density_plot_z)


############################################################################
# Louvain clustering

# Load required packages
library(igraph)
library(ggplot2)
library(gridExtra)

# ------------------------------
# Assumptions:
#   - PCA_ready_data: gene expression matrix (rows: patients, columns: genes)
#   - group_labels: factor vector (levels: "Oxford", "UCL", "UHB", "Controls")
# ------------------------------

# Compute the correlation among patients
patient_corr <- cor(t(PCA_ready_data), use = "pairwise.complete.obs")

# Define the thresholds to test
thresholds <- c(0.5, 0.7, 0.9, 0.91, 0.93, 0.95)

# Prepare a data frame to store results
results <- data.frame(threshold = thresholds,
                      num_communities = NA,
                      modularity = NA)

# Loop over each threshold to perform Louvain clustering
for(i in seq_along(thresholds)){
  th <- thresholds[i]
  
  # Create an adjacency matrix: set correlations below the threshold to 0
  adj <- patient_corr
  adj[adj < th] <- 0
  
  # Build the graph from the thresholded matrix
  g <- graph_from_adjacency_matrix(adj, mode = "undirected", weighted = TRUE, diag = FALSE)
  
  # Perform Louvain community detection
  lc <- cluster_louvain(g)
  
  # Record the number of communities and the modularity
  results$num_communities[i] <- length(lc)
  results$modularity[i] <- modularity(lc, weights = E(g)$weight)
}

# Print the results
print(results)

# Create two visual plots using ggplot2:
# Plot 1: Number of Communities vs. Correlation Threshold
p1 <- ggplot(results, aes(x = threshold, y = num_communities)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "blue", size = 3) +
  labs(title = "Number of Communities vs. Correlation Threshold",
       x = "Correlation Threshold", y = "Number of Communities") +
  theme_minimal()

# Plot 2: Modularity vs. Correlation Threshold
p2 <- ggplot(results, aes(x = threshold, y = modularity)) +
  geom_line(color = "darkgreen", size = 1) +
  geom_point(color = "darkgreen", size = 3) +
  labs(title = "Modularity vs. Correlation Threshold",
       x = "Correlation Threshold", y = "Modularity") +
  theme_minimal()

# Arrange the two plots side-by-side
grid.arrange(p1, p2, ncol = 2)

##########################################
# Louvain clustering performed for 0.91 - 0.95

# 1. Compute the correlation matrix among patients
patient_corr <- cor(t(PCA_ready_data), use = "pairwise.complete.obs")

# 2. Define thresholds from 0.91 to 0.95 (increments of 0.01)
thresholds <- seq(0.91, 0.95, by = 0.01)

# 3. Prepare a data frame to store the number of communities & modularity
results <- data.frame(threshold = thresholds,
                      num_communities = NA,
                      modularity = NA)

# 4. Combine Oxford, UCL, and UHB into one group "BioAID"; keep "Controls" separate
combined_groups <- ifelse(as.character(group_labels) == "Controls", 
                          "Controls", 
                          "BioAID")

# 5. Perform Louvain clustering at each threshold
louvain_tables <- list()  # to store the "Community vs. Combined Groups" tables

for(i in seq_along(thresholds)) {
  th <- thresholds[i]
  
  # Create an adjacency matrix for this threshold
  adjacency <- patient_corr
  adjacency[adjacency < th] <- 0
  
  # Build the graph
  g <- graph_from_adjacency_matrix(adjacency, mode = "undirected", 
                                   weighted = TRUE, diag = FALSE)
  
  # Perform Louvain community detection
  louvain_result <- cluster_louvain(g)
  
  # Record the number of communities & modularity
  results$num_communities[i] <- length(louvain_result)
  results$modularity[i] <- modularity(louvain_result, weights = E(g)$weight)
  
  # Compare each community to the combined groups (BioAID vs. Controls)
  comm_vs_group <- table(
    Community = membership(louvain_result),
    Group = combined_groups
  )
  
  # Store this table so we can look at it later
  louvain_tables[[as.character(th)]] <- comm_vs_group
  
  # Print the table for clarity
  cat("\n=================================================\n")
  cat("Threshold =", th, "\n")
  cat("Number of Communities:", length(louvain_result), "\n")
  cat("Modularity:", results$modularity[i], "\n")
  cat("\nCommunity vs. BioAID/Controls:\n")
  print(comm_vs_group)
}

# 6. Visualize the Number of Communities and Modularity across thresholds
p1 <- ggplot(results, aes(x = threshold, y = num_communities)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "blue", size = 3) +
  labs(title = "Number of Communities vs. Correlation Threshold",
       x = "Correlation Threshold", y = "Number of Communities") +
  theme_minimal()

p2 <- ggplot(results, aes(x = threshold, y = modularity)) +
  geom_line(color = "darkgreen", size = 1) +
  geom_point(color = "darkgreen", size = 3) +
  labs(title = "Modularity vs. Correlation Threshold",
       x = "Correlation Threshold", y = "Modularity") +
  theme_minimal()

# Arrange both plots side by side
grid.arrange(p1, p2, ncol = 2)

###########################################################################
# Deep dive analysis of 0.95 Louvain cluster
library(igraph)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(pheatmap)

# 1. Extract membership from Louvain clustering
community_membership <- membership(louvain_result)

# 2. Find the 3 largest communities
comm_sizes <- table(community_membership)
largest_3 <- names(sort(comm_sizes, decreasing = TRUE))[1:3]

# 3. Subset your data and define a factor
samples_in_largest_3 <- names(community_membership)[community_membership %in% largest_3]
data_largest_3 <- PCA_ready_data[samples_in_largest_3, ]
community_factor <- factor(community_membership[samples_in_largest_3], 
                           levels = largest_3)

# 4. ANOVA to find top discriminating genes
pvals <- apply(
  data_largest_3, 
  2, 
  function(gene_expr) {
    fit <- aov(gene_expr ~ community_factor)
    summary(fit)[[1]][["Pr(>F)"]][1]
  }
)
pvals_adj <- p.adjust(pvals, method = "BH")
sorted_genes <- names(sort(pvals_adj, decreasing = FALSE))
top_genes <- head(sorted_genes, 50)

# 5. Heatmap of top 50 genes
top_genes_expr <- data_largest_3[, top_genes]
top_genes_z <- t(scale(t(top_genes_expr)))  # row-wise z-scaling
annotation_row <- data.frame(Community = community_factor)
rownames(annotation_row) <- rownames(top_genes_z)

pheatmap(
  top_genes_z,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_row = annotation_row,
  show_rownames = FALSE,
  main = "Heatmap of Top 50 Genes Distinguishing the 3 Largest Communities"
)

# 6. PCA on top genes
pca_top <- prcomp(top_genes_expr, center = TRUE, scale. = TRUE)
pca_top_df <- data.frame(
  PC1 = pca_top$x[,1],
  PC2 = pca_top$x[,2],
  Community = community_factor
)

ggplot(pca_top_df, aes(x = PC1, y = PC2, color = Community)) +
  geom_point(size = 3) +
  # geom_text_repel(aes(label = rownames(pca_top_df)), size = 3) + 
  labs(
    title = "PCA (Top Discriminating Genes) for 3 Largest Louvain Communities",
    x = "PC1",
    y = "PC2"
  ) +
  theme_pubr()

##########################################################
# K-Means Clustering Pipeline (Trial & Error: k = 4 to 8)
##########################################################
# 3. Remove Constant Columns (zero variance) to avoid NAs during scaling
variances <- apply(PCA_ready_data, 2, var, na.rm = TRUE)
PCA_ready_data_filtered <- PCA_ready_data[, variances > 0]

# Optionally, remove rows with NA values (if any)
PCA_ready_data_filtered <- na.omit(PCA_ready_data_filtered)

# 4. Scale Data (recommended for K-means)
scaled_data <- scale(PCA_ready_data_filtered)

# Double-check for issues:
if (sum(is.na(scaled_data)) > 0 || sum(is.infinite(scaled_data)) > 0) {
  stop("Scaled data contains NA/Inf values. Please check your input data.")
}

# 5. Try K-Means from k = 2 to k = 8 (Trial & Error)
set.seed(123)  # for reproducibility
k_values <- 2:8
wss_values <- numeric(length(k_values))  # store total within-cluster sum of squares

for (i in seq_along(k_values)) {
  k <- k_values[i]
  km_out <- kmeans(scaled_data, centers = k, nstart = 25)
  wss_values[i] <- km_out$tot.withinss
}

# 6. Plot Elbow Curve
wss_df <- data.frame(k = k_values, wss = wss_values)
elbow_plot <- ggplot(wss_df, aes(x = k, y = wss)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "blue", size = 3) +
  labs(
    title = "Elbow Plot for K-Means Clustering (k = 2 to 8)",
    x = "Number of Clusters (k)",
    y = "Total Within-Cluster Sum of Squares"
  ) +
  theme_minimal()
print(elbow_plot)

##############################
# Deep-dive analysis of K = 3 - 8

# 1. Load Required Packages
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(dplyr)
library(pheatmap)
library(scales)

# 2. Assume 'PCA_ready_data' is your gene expression matrix 
#    (Rows = samples, Columns = genes)
#    Also assume 'group_labels' is a vector of sample group labels 
#    (with values "Oxford", "UCL", "UHB", or "Controls") that aligns with the rownames of PCA_ready_data.

# 3. Remove Constant Columns and Rows with NAs
variances <- apply(PCA_ready_data, 2, var, na.rm = TRUE)
PCA_ready_data_filtered <- PCA_ready_data[, variances > 0]
PCA_ready_data_filtered <- na.omit(PCA_ready_data_filtered)

# 4. Scale Data (recommended for K-means)
scaled_data <- scale(PCA_ready_data_filtered)

# Double-check scaled data for issues:
if (sum(is.na(scaled_data)) > 0 || sum(is.infinite(scaled_data)) > 0) {
  stop("Scaled data contains NA/Inf values. Please check your input data.")
}

# 5. Define the Combined Group:
#    BioAID is an amalgamation of Oxford, UCL, and UHB. Controls remain separate.
combined_group <- ifelse(group_labels == "Controls", "Controls", "BioAID")
combined_group <- factor(combined_group, levels = c("Controls", "BioAID"))

# 6. Run K-Means for k = 3 to 8 and store results for comparison
plot_data_list <- list()

for (k in 3:8) {
  set.seed(123)  # for reproducibility
  km_out <- kmeans(scaled_data, centers = k, nstart = 25)
  clusters <- km_out$cluster
  
  # Create and print a contingency table: Cluster vs. Combined Group
  cont_tab <- table(Cluster = clusters, Group = combined_group)
  cat("\n---------------------\n")
  cat("K =", k, "\n")
  print(cont_tab)
  
  # Build a data frame with sample-level info for plotting
  df_plot <- data.frame(
    Sample = rownames(scaled_data),
    Cluster = factor(clusters),
    Group = combined_group,
    k = paste0("k=", k)
  )
  
  plot_data_list[[as.character(k)]] <- df_plot
}

# Combine the data frames for all k values
all_plot_data <- do.call(rbind, plot_data_list)

# 7. Create a Faceted Bar Plot (stacked percentages) to compare cluster composition
bar_plot <- ggplot(all_plot_data, aes(x = Cluster, fill = Group)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = percent) +
  labs(
    title = "Cluster Composition (Controls vs. BioAID) for k = 3 to 8",
    x = "K-Means Cluster",
    y = "Percentage of Samples",
    fill = "Group"
  ) +
  facet_wrap(~ k, nrow = 2) +
  theme_minimal()

print(bar_plot)



