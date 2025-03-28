# Ensure rownames in matched_patients align with plot_info
matched_expression <- matched_patients[rownames(plot_info), ]

# Extract micro_diagnosis for each patient
micro_diag_map <- plot_info %>%
  select(original_id, micro_diagnosis) %>%
  deframe()

# Function to parse micro_diagnosis
parse_diagnosis <- function(x) {
  if (is.na(x)) return(character(0))
  x <- gsub(" ", "", x)
  unlist(strsplit(x, "&"))
}

# Main loop over PCC and resolution thresholds
for (pcc_threshold in c(0.91, 0.92, 0.93, 0.94, 0.95)) {
  for (resolution in seq(0.1, 1.0, 0.1)) {
    
    # Compute correlation matrix
    corr_matrix <- cor(t(matched_expression), method = "pearson")
    corr_matrix[is.na(corr_matrix)] <- 0  # Replace NAs
    
    # Threshold the correlation matrix to create adjacency matrix
    adj_matrix <- corr_matrix
    adj_matrix[abs(adj_matrix) < pcc_threshold] <- 0
    diag(adj_matrix) <- 0
    
    # Create graph from adjacency matrix
    g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)
    
    # Apply Louvain clustering with multilevel.community + resolution scaling
    cluster_result <- cluster_louvain(g)
    
    # Get communities
    community_list <- split(names(membership(cluster_result)), membership(cluster_result))
    
    # Filter out communities with <5 patients
    community_list <- community_list[sapply(community_list, length) > 4]
    
    
    for (i in seq_along(community_list)) {
      community <- community_list[[i]]
      diagnoses <- unlist(lapply(micro_diag_map[community], parse_diagnosis))
      diagnosis_counts <- sort(table(diagnoses), decreasing = TRUE)
    }
  }
}

# Plot 1: % of nodes captured in communities > 4
plot_percent <- ggplot(purity_summary, aes(x = Resolution, y = Percent_Captured, group = PCC, color = as.factor(PCC))) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  labs(
    title = "% of Nodes in Communities > 4",
    x = "Resolution",
    y = "% of Nodes Captured",
    color = "PCC Threshold"
  ) +
  theme_minimal()

# Plot 2: Weighted Purity of Louvain Clusters
plot_purity <- ggplot(purity_summary, aes(x = Resolution, y = Weighted_Purity, group = PCC, color = as.factor(PCC))) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  labs(
    title = "Weighted Purity of Louvain Clusters",
    x = "Resolution",
    y = "Weighted Purity",
    color = "PCC Threshold"
  ) +
  theme_minimal()

# Print side by side
grid.arrange(plot_percent, plot_purity, ncol = 2)
