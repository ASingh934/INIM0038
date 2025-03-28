# Define PCC thresholds and resolution limits
pcc_thresholds <- c(0.6, 0.7, 0.8, 0.85, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96)
resolution_limits <- seq(0.1, 2.0, by = 0.1)

# Prepare patient micro_diagnosis lookup table
diagnosis_lookup <- clinical_match %>%
  select(matched_id, micro_diagnosis) %>%
  distinct()

# Match to cleaned IDs
id_lookup <- gene_ids %>%
  filter(original_id %in% matched_patient_ids) %>%
  select(original_id, cleaned_id)

diagnosis_lookup <- left_join(id_lookup, diagnosis_lookup, by = c("cleaned_id" = "matched_id"))

# Standardise again to ensure gene expression data is scaled between -1 to 1
matched_patients_scaled <- matched_patients

# Create a list to collect all summary tables
all_summary_tables <- list()

# Main clustering loop
for (pcc in pcc_thresholds) {
  for (res in resolution_limits) {
    
    # Step 1: Build adjacency graph using PCC threshold
    cor_matrix <- cor(t(matched_patients_scaled), method = "pearson")
    adj_matrix <- cor_matrix
    adj_matrix[abs(adj_matrix) < pcc] <- 0
    diag(adj_matrix) <- 0  # remove self-loops
    
    # Convert to igraph object
    graph <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)
    
    # Step 2: Perform Leiden clustering
    set.seed(123)  # Ensures consistent clustering results
    leiden_clusters <- leiden_find_partition(graph, partition_type = "RBConfiguration", resolution_parameter = res)
    
    # Assign community labels
    V(graph)$community <- leiden_clusters$membership
    comm_df <- data.frame(
      patient_id = names(V(graph)),
      community = V(graph)$community
    )
    
    # Merge with diagnosis
    comm_df <- left_join(comm_df, diagnosis_lookup, by = c("patient_id" = "original_id"))
    
    # Exclude communities smaller than 20
    comm_counts <- comm_df %>% count(community)
    valid_communities <- comm_counts %>% filter(n >= 20) %>% pull(community)
    comm_df <- comm_df %>% filter(community %in% valid_communities)
    
    # Summary table for this configuration
    summary_table <- comm_df %>%
      rowwise() %>%
      mutate(micro_diagnosis_split = str_split(micro_diagnosis, " & ")) %>%
      unnest(micro_diagnosis_split) %>%
      group_by(community, micro_diagnosis_split) %>%
      summarise(count = n(), .groups = "drop") %>%
      pivot_wider(names_from = micro_diagnosis_split, values_from = count, values_fill = 0) %>%
      left_join(comm_counts, by = "community") %>%
      rename(Patient_Count = n) %>%
      arrange(desc(Patient_Count))
    
    # Append to collection if valid
    if (nrow(summary_table) > 0) {
      summary_table$PCC_Threshold <- pcc
      summary_table$Resolution_Limit <- res
      all_summary_tables[[length(all_summary_tables) + 1]] <- summary_table
    }
  }
}

# Combine all summary tables into one data frame
if (length(all_summary_tables) > 0) {
  combined_summary <- do.call(rbind, all_summary_tables)
  
  # Reorder columns for readability
  combined_summary <- combined_summary %>%
    relocate(PCC_Threshold, Resolution_Limit, community, Patient_Count)
} else {
  combined_summary <- NULL
}
