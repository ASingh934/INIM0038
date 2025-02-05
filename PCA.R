file.rename("Untitled1.R", "PCA.R")

library(ggplot2)
library(ggbiplot)
library(dplyr)
library(ggrepel)
library(ggpubr)

gene_data <- read.csv("C:/Users/ABSin/Downloads/BIOAID_tpm_PC0.001_log2_genesymbol_dedup (1).csv")

################################################################################

# Transposes gene data so that patients are observations and genes are variables
gene_data <- as.data.frame( t(gene_data) )

###############################################################################

# The column headings need to be changed to the gene names

# The column headings (gene names) are stored
new_headers <- gene_data[1,] #

# Removes the first row (gene names) from this dataframe
gene_data <- gene_data[-1,] 

# The column headings are replaced with new_headers
colnames(gene_data) <- new_headers 

gene_data[] <- lapply(gene_data, as.numeric) # Treats str characters as numeric

##############################################################################

# Filter out non-protein coding genes

# Imports the csv file containing the gene symbols and if they are protein-coding
annotated_data <- read.csv("C:/Users/ABSin/Downloads/Annotatedexample_tpm_PC0.001_log2.csv")

# Filter annotation data to retain only protein-coding genes
valid_genes <- annotated_data %>%
  filter(gene_biotype == "protein_coding") %>%
  select(external_gene_name) %>%
  pull()

# Keep only columns in gene_data that match the valid gene names
gene_data <- gene_data %>% select(all_of(intersect(names(gene_data), valid_genes)))

###############################################################################

# Prepares the gene_data for PCA
PCA_ready_data <- gene_data

###############################################################################

# Perform PCA
pca_result <- prcomp(PCA_ready_data, center = TRUE, scale. = FALSE)

###############################################################################

# Summary of PCA - prints the SD & proportion of variance of each PC
print(summary(pca_result)) 

###############################################################################

# Express PC1 as a linear combination of its dimensions, with only the 50
# most weighted dimensions displayed

# Extract PC1 loadings (eigenvector coefficients for PC1)
pc1_loadings <- pca_result$rotation[,1]  # First principal component

# Order the coefficients by absolute magnitude in descending order
sorted_indices <- order(abs(pc1_loadings), decreasing = TRUE)
sorted_pc1_loadings <- pc1_loadings[sorted_indices]

# Select the top 50 most weighted genes
top_50_indices <- sorted_indices[1:50]  # Get indices of top 50 most weighted genes
top_50_loadings <- pc1_loadings[top_50_indices]

# Format the equation with the largest coefficients first
pc1_equation <- paste0(round(top_50_loadings, 4), " * ", names(top_50_loadings), collapse = " + ")

# Print the PC1 equation (Top 50 genes only)
cat("\nPC1 is expressed as:\n")
cat("PC1 =", pc1_equation, "\n")

################################################################################

# Sort the patients into groups based on whether they are controls or from OX, UHB or UCL.

# Extract sample IDs from row names
sample_ids <- rownames(PCA_ready_data)

# Create a grouping variable based on sample IDs
group_labels <- case_when(
  grepl("^(OX|ox)", sample_ids) ~ "Oxford",
  grepl("^UP", sample_ids) ~ "UCL",
  grepl("^(667|X)", sample_ids) ~ "UHB",
  grepl("^WH0", sample_ids) ~ "Controls"
)

# Convert to factor for consistent ordering
group_labels <- factor(group_labels, levels = c("Oxford", "UCL", "UHB", "Controls"))

################################################################################

# Visualise PCA using a PCA plot

pca_plot <- ggbiplot(
  pca_result,
  obs.scale = 1, 
  var.scale = 1,
  ellipse = FALSE, 
  circle = FALSE,
  var.axes = FALSE,  
  labels = NULL
) +
  geom_point(aes(color = group_labels), size = 1, shape = 16, stroke = 0) +
  scale_color_manual(values = c("Oxford" = "blue", 
                                "UCL" = "red", 
                                "UHB" = "green", 
                                "Controls" = "purple")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  labs(title = "PCA of Gene Expression Data", color = "Group") +
  theme_minimal()


print(pca_plot)

###############################################################################

# Scree plot to visualize the proportion of total variance captured by the first 50 PCs

# Variance explained by each principal component
explained_variance <- pca_result$sdev^2 / sum(pca_result$sdev^2)

# Create a data frame for all PCs
scree_data <- data.frame(
  PC       = 1:length(explained_variance), 
  Variance = explained_variance
)

# Keep only the first 50 PCs
scree_data_50 <- scree_data[1:50, ]

scree_plot <- ggplot(scree_data_50, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  labs(title = "Scree Plot (First 50 PCs)", 
       x = "Principal Components", 
       y = "Variance Explained") +
  theme_minimal()

print(scree_plot)






