# Perform PCA
pca_result <- prcomp(gene_data, center = TRUE, scale. = FALSE)

###############################################################################

# Summary of PCA - prints the SD & proportion of variance of each PC
print(summary(pca_result)) 

###############################################################################

# Express PC1 as a linear combination of its dimensions (top 50 weighted genes)
pc1_loadings <- pca_result$rotation[,1]
sorted_indices <- order(abs(pc1_loadings), decreasing = TRUE)
top_50_indices <- sorted_indices[1:50]
top_50_loadings <- pc1_loadings[top_50_indices]
pc1_equation <- paste0(round(top_50_loadings, 4), " * ", names(top_50_loadings), 
                       collapse = " + ")

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
