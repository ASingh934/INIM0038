# ---------------------------------------------------------------
# NEW CODE: Louvain clustering over varying PCC thresholds & resolution multipliers
# ---------------------------------------------------------------
set.seed(123)  # For reproducibility

# Use the matched_patients and plot_info data frames already created in your code.
# matched_patients: gene expression data for patients present in both datasets.
# plot_info: clinical data (with rownames as the original patient IDs) containing micro_diagnosis.

n_patients <- nrow(matched_patients)
sample_size <- floor(0.8 * n_patients)

# Randomly sample 80% of the matched patients
sampled_ids <- sample(rownames(matched_patients), size = sample_size)
sampled_data <- matched_patients[sampled_ids, ]

# Define the PCC thresholds and resolution multipliers
pcc_thresholds <- c(0.85, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99)
resolution_multipliers <- seq(0, 1, by = 0.1)

# Loop over each combination of threshold and resolution multiplier
for (thresh in pcc_thresholds) {
  for (res in resolution_multipliers) {
    cat("------------------------------------------------------\n")
    cat("PCC Threshold:", thresh, "| Resolution multiplier:", res, "\n")
    
    # Compute the Pearson correlation matrix between patients.
    # Since rows are patients and columns are genes, we compute correlation on the transpose.
    corr_matrix <- cor(t(sampled_data), method = "pearson")
    
    # Threshold the matrix: set values below the threshold to 0.
    corr_matrix[corr_matrix < thresh] <- 0
    
    # Remove self-loops by zeroing the diagonal.
    diag(corr_matrix) <- 0
    
    # Apply the resolution multiplier to simulate a resolution effect on the edge weights.
    weighted_matrix <- corr_matrix * res
    
    # Create an undirected weighted graph using the weighted matrix.
    # Only nonzero edges will be included; vertex names come from the rownames of sampled_data.
    g <- graph.adjacency(weighted_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)
    
    # Remove isolated nodes (if any)
    g <- delete.vertices(g, which(degree(g) == 0))
    
    # If the graph is empty for this combination, skip to the next iteration.
    if (vcount(g) == 0) {
      cat("Graph is empty for this combination. Skipping...\n")
      next
    }
    
    # Perform Louvain clustering using the edge weights.
    louvain_result <- cluster_louvain(g, weights = E(g)$weight)
    total_communities <- length(louvain_result)
    cat("Total communities detected:", total_communities, "\n")
    
    # Get the sizes of each community.
    comm_sizes <- sizes(louvain_result)
    
    # Identify the three largest communities by size.
    top3_ids <- names(sort(comm_sizes, decreasing = TRUE))[1:3]
    
    # For each of the top 3 communities, print out the community size and the micro_diagnosis breakdown.
    for (i in seq_along(top3_ids)) {
      comm_id <- top3_ids[i]
      
      # membership() returns a named vector where names are patient IDs (from matched_patients) 
      # and values are community IDs.
      members <- names(membership(louvain_result))[membership(louvain_result) == as.numeric(comm_id)]
      
      # Use plot_info to extract clinical micro_diagnosis for the members.
      # (Since plot_info rownames were set to the original IDs, they should match the vertex names.)
      diag_info <- plot_info[members, "micro_diagnosis", drop = FALSE]
      
      # Create a table for the micro_diagnosis breakdown (Viral, Bacterial, None, etc.)
      diag_table <- table(diag_info$micro_diagnosis, useNA = "ifany")
      
      cat("\nTop community", i, "(Community ID:", comm_id, 
          "- Size:", comm_sizes[comm_id], ")\n")
      print(diag_table)
    }
    
    cat("------------------------------------------------------\n\n")
  }
}
