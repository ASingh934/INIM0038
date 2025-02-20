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
