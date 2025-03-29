###############################################################################

# Leiden 

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

# 3️⃣ Print both plots
gridExtra::grid.arrange()

# 4️⃣ Identify the best combo based on both metrics
# Normalize both metrics to [0, 1], then average them
combined_scores <- percent_df %>%
  inner_join(wp_df, by = c("PCC_Threshold", "Resolution_Limit")) %>%
  mutate(Norm_Percent = (Percentage_Captured - min(Percentage_Captured)) / (max(Percentage_Captured) - min(Percentage_Captured)),
         Norm_Purity  = (Weighted_Purity - min(Weighted_Purity)) / (max(Weighted_Purity) - min(Weighted_Purity)),
         Combined_Score = (Norm_Percent + Norm_Purity) / 2)

best_row <- combined_scores %>% slice_max(order_by = Combined_Score, n = 1)

cat("\n🔥 Best combination (based on combined WP and patient capture):\n")
cat("PCC Threshold:", best_row$PCC_Threshold, "\n")
cat("Resolution Limit:", best_row$Resolution_Limit, "\n")
cat("Weighted Purity:", round(best_row$Weighted_Purity, 4), "\n")
cat("Patient Capture (%):", round(best_row$Percentage_Captured, 2), "\n")
