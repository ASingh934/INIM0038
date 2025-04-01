###############################################################################

# Leiden clustering

#############################################################################

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

# Standardize again to ensure gene expression data is scaled between -1 to 1
matched_patients_scaled <- matched_patients

# 🔒 Ensure reproducible row order for PCC and graph
matched_patients_scaled <- matched_patients_scaled[order(rownames(matched_patients_scaled)), ]

# Create lists to collect all summary tables and all patient assignments
all_summary_tables <- list()
all_assignments <- list()

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
    set.seed(123)  # Ensures consistent R-side results
    reticulate::py_run_string("import random; import numpy as np; random.seed(123); np.random.seed(123)")
    leiden_clusters <- leiden_find_partition(graph, partition_type = "RBConfiguration", resolution_parameter = res)
    
    # Assign community labels to graph nodes
    V(graph)$community <- leiden_clusters$membership
    
    # Create a dataframe of patient assignments (patient IDs are row names in gene_data)
    comm_df <- data.frame(
      patient_id = names(V(graph)),
      community = V(graph)$community,
      stringsAsFactors = FALSE
    )
    
    # Save the patient assignments in a list using a key based on the parameters
    key <- paste0("pcc_", pcc, "_res_", res)
    all_assignments[[key]] <- comm_df
    
    # Merge with diagnosis for summary generation
    comm_df_summary <- left_join(comm_df, diagnosis_lookup, by = c("patient_id" = "original_id"))
    
    # Exclude communities smaller than 20 patients
    comm_counts <- comm_df_summary %>% count(community)
    valid_communities <- comm_counts %>% filter(n >= 20) %>% pull(community)
    comm_df_summary <- comm_df_summary %>% filter(community %in% valid_communities)
    
    # Create summary table for this configuration
    summary_table <- comm_df_summary %>%
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

# Combine all summary tables into one data-frame
if (length(all_summary_tables) > 0) {
  combined_summary <- dplyr::bind_rows(all_summary_tables) %>%
    relocate(PCC_Threshold, Resolution_Limit, community, Patient_Count)
} else {
  combined_summary <- NULL
}

# Print this combined_summary to get a list of all the Leiden clusters
# with a diagnosis breakdown of each community

###############################################################################

# Plotting code 

#############################################################################

# Total number of patients
total_patients <- 362

# 1️⃣ Percent Captured Plot ------------------------------------------
# Compute percentage of total patients captured for each PCC/res combination
percent_df <- combined_summary %>%
  group_by(PCC_Threshold, Resolution_Limit) %>%
  summarise(Patients_Captured = sum(Patient_Count), .groups = "drop") %>%
  mutate(Percentage_Captured = round(100 * Patients_Captured / total_patients, 2))

plot1 <- ggplot(percent_df, aes(x = Resolution_Limit, y = Percentage_Captured, color = factor(PCC_Threshold))) +
  geom_line(size = 1) +
  geom_point(size = 1.5) +
  labs(title = "Percentage of Patients Captured",
       x = "Resolution Parameter",
       y = "% of Patients (≥20 size communities)",
       color = "PCC Threshold") +
  theme_minimal()

# 2️⃣ Weighted Purity Plot ------------------------------------------

# Calculate purity per community
purity_df <- combined_summary %>%
  pivot_longer(cols = -c(PCC_Threshold, Resolution_Limit, community, Patient_Count),
               names_to = "Diagnosis", values_to = "Count") %>%
  group_by(PCC_Threshold, Resolution_Limit, community) %>%
  mutate(Purity = Count / sum(Count)) %>%
  slice_max(order_by = Purity, n = 1, with_ties = FALSE) %>%
  ungroup()

# Compute WP for each PCC+res
wp_df <- purity_df %>%
  group_by(PCC_Threshold, Resolution_Limit) %>%
  summarise(Weighted_Purity = sum(Purity * Patient_Count) / sum(Patient_Count), .groups = "drop")

plot2 <- ggplot(wp_df, aes(x = Resolution_Limit, y = Weighted_Purity, color = factor(PCC_Threshold))) +
  geom_line(size = 1) +
  geom_point(size = 1.5) +
  labs(title = "Weighted Purity of Communities",
       x = "Resolution Parameter",
       y = "Weighted Purity",
       color = "PCC Threshold") +
  theme_minimal()

# 3️⃣ Kappa Concordant statistic Plot (Viral) ----------------------------
# Create a dataframe to store Kappa values
kappa_results <- data.frame()

# Loop through saved assignments
for (key in names(all_assignments)) {
  
  comm_df <- all_assignments[[key]]
  
  # Extract PCC and Resolution from key
  key_parts <- strsplit(key, "_")[[1]]
  pcc <- as.numeric(key_parts[2])
  res <- as.numeric(key_parts[4])
  
  # Join with diagnosis info
  labeled_df <- left_join(comm_df, diagnosis_lookup, by = c("patient_id" = "original_id"))
  
  # Skip if too few labeled patients
  if (nrow(labeled_df) < 20) next
  
  # Count community sizes
  comm_counts <- labeled_df %>% count(community)
  valid_communities <- comm_counts %>% filter(n >= 20) %>% pull(community)
  
  if (length(valid_communities) == 0) next
  
  # Filter to valid communities only
  labeled_df <- labeled_df %>% filter(community %in% valid_communities)
  
  # Determine number of viral patients in each community
  labeled_df <- labeled_df %>%
    mutate(is_viral = grepl("Viral", micro_diagnosis, ignore.case = TRUE))
  
  viral_counts <- labeled_df %>%
    group_by(community) %>%
    summarise(n_viral = sum(is_viral), .groups = "drop")
  
  # Viral community = community with most viral patients
  viral_community <- viral_counts %>%
    slice_max(n_viral, with_ties = FALSE) %>%
    pull(community)
  
  # Predicted: in viral_community => Viral, else Non-Viral
  labeled_df <- labeled_df %>%
    mutate(
      predicted_label = ifelse(community == viral_community, "Viral", "Non-Viral"),
      true_label = ifelse(is_viral, "Viral", "Non-Viral")
    )
  
  # Prepare table for Kappa
  if (length(unique(labeled_df$predicted_label)) > 1 && length(unique(labeled_df$true_label)) > 1) {
    kappa_score <- kappa2(labeled_df[, c("predicted_label", "true_label")])$value
  } else {
    kappa_score <- NA  # Undefined if one label only
  }
  
  # Store result
  kappa_results <- rbind(kappa_results, data.frame(
    PCC_Threshold = pcc,
    Resolution_Limit = res,
    Kappa = kappa_score
  ))
}

# Plot Kappa results
plot3 <- ggplot(kappa_results, aes(x = Resolution_Limit, y = Kappa, color = factor(PCC_Threshold))) +
  geom_line(size = 1) +
  geom_point(size = 1.5) +
  labs(
    title = "Cohen's Kappa Concordance for Leiden Clusters",
    x = "Resolution Parameter",
    y = "Kappa Statistic",
    color = "PCC Threshold"
  ) +
  theme_minimal()

# ---------------------------------------------------------------------------

# Print all three plots
print(plot1)
print(plot2)
print(plot3)

# ---- Print Optimum Leiden Cluster ----------------------------------------------
# Combine all evaluation metrics
eval_df <- kappa_results %>%
  inner_join(wp_df, by = c("PCC_Threshold", "Resolution_Limit")) %>%
  inner_join(percent_df, by = c("PCC_Threshold", "Resolution_Limit")) %>%
  filter(Percentage_Captured >= 75)

# Find max Kappa
max_kappa <- max(eval_df$Kappa, na.rm = TRUE)

# Keep rows within 0.05 of max Kappa
top_kappa_rows <- eval_df %>%
  filter(Kappa >= max_kappa - 0.05)

# Select best row by highest weighted purity
best_config <- top_kappa_rows %>%
  slice_max(order_by = Weighted_Purity, n = 1, with_ties = FALSE)

# Print best configuration
cat("\n✅ Optimum Leiden Clustering Configuration:\n")
cat("PCC Threshold:", best_config$PCC_Threshold, "\n")
cat("Resolution Limit:", best_config$Resolution_Limit, "\n")
cat("Kappa:", round(best_config$Kappa, 4), "\n")
cat("Weighted Purity:", round(best_config$Weighted_Purity, 4), "\n")
cat("Patient Capture (%):", round(best_config$Percentage_Captured, 2), "\n")
