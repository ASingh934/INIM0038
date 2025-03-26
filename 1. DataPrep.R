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

##########
# 1. Fixed correlation threshold and Pearson's correlation
cor_method <- "pearson"
threshold <- 0.95

# Compute the correlation matrix (genes are in columns)
patient_corr <- cor(t(gene_data), use = "pairwise.complete.obs", method = cor_method)

# Threshold the matrix: set correlations below threshold to 0
adj <- patient_corr
adj[adj < threshold] <- 0

# Build the undirected, weighted graph
g <- graph_from_adjacency_matrix(adj, mode = "undirected", weighted = TRUE, diag = FALSE)

# 2. Define the vectorized custom modularity function with resolution parameter gamma
compute_modularity_custom <- function(g, membership, gamma) {
  W <- as.matrix(as_adjacency_matrix(g, attr = "weight"))
  m <- sum(W) / 2
  s <- rowSums(W)
  
  # Expected edge weights with resolution parameter gamma
  X <- W - gamma * (outer(s, s)) / (2 * m)
  
  # Membership indicator matrix (one column per community)
  communities <- sort(unique(membership))
  M <- sapply(communities, function(c) as.numeric(membership == c))
  
  # Compute modularity: sum of within-community contributions divided by (2*m)
  Q <- sum(diag(t(M) %*% X %*% M)) / (2 * m)
  return(Q)
}

# 3. Vary the resolution parameter using cluster_leiden() for community detection
# Define a sequence of resolution limits (gamma)
resolutions <- seq(0.5, 2, by = 0.1)
results <- data.frame(resolution = resolutions,
                      num_communities = NA,
                      modularity = NA)

for (i in seq_along(resolutions)) {
  gamma <- resolutions[i]
  
  # Perform community detection using the Leiden algorithm with the given resolution
  lc <- cluster_leiden(g, resolution_parameter = gamma)
  membership <- membership(lc)
  
  results$num_communities[i] <- length(unique(membership))
  results$modularity[i] <- compute_modularity_custom(g, membership, gamma)
}

# 4. Plotting

# Plot: Number of Communities vs. Resolution Limit (γ)
p_num_comm <- ggplot(results, aes(x = resolution, y = num_communities)) +
  geom_line() +
  geom_point() +
  labs(title = "Number of Communities vs. Resolution Limit",
       x = "Resolution Limit (γ)",
       y = "Number of Communities") +
  theme_minimal()

# Plot: Modularity vs. Resolution Limit (γ)
p_modularity <- ggplot(results, aes(x = resolution, y = modularity)) +
  geom_line() +
  geom_point() +
  labs(title = "Modularity vs. Resolution Limit",
       x = "Resolution Limit (γ)",
       y = "Modularity") +
  theme_minimal()

# Display the two plots side-by-side
grid.arrange(p_num_comm, p_modularity, ncol = 2)
