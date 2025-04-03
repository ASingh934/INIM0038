################################################################################

# Perform K-means clustering from K = 1 to K = 15

################################################################################

# Remove any genes with zero variance (recheck just to be safe)
matched_patients_filtered <- matched_patients[, apply(matched_patients, 2, function(x) sd(x) > 0)]

# Scale data (mean 0, SD 1 per gene)
scaled_data <- scale(matched_patients_filtered)

# Check for NaNs in scaled version
anyNA(scaled_data)

k_values <- 1:15
kmeans_results <- list()
wss_values <- numeric(length(k_values))

for (k in k_values) {
  set.seed(123)
  kmeans_fit <- kmeans(scaled_data, centers = k, nstart = 25)
  kmeans_results[[as.character(k)]] <- kmeans_fit
  wss_values[k] <- kmeans_fit$tot.withinss
}

################################################################################

# Elbow plot

################################################################################

# Rescale WSS values to millions
elbow_df <- data.frame(K = k_values, WSS = wss_values / 1e6)

# Plot with adjusted axis label and values
elbow_plot <- ggplot(elbow_df, aes(x = K, y = WSS)) +
  geom_line(color = "steelblue", size = 1.2) +
  geom_point(size = 2, color = "red") +
  labs(
    title = "Elbow Plot: K-Means Clustering",
    x = "Number of Clusters (K)",
    y = expression("Total Within-Cluster Sum of Squares (×10 "^6*")")
  ) +
  theme_minimal()

print(elbow_plot)

################################################################################

# Cluster composition summary table

################################################################################

# Step 1: Create diagnosis_lookup
diagnosis_lookup <- gene_ids %>%
  filter(original_id %in% matched_patient_ids) %>%
  select(original_id, cleaned_id) %>%
  left_join(
    clinical_match %>% select(matched_id, micro_diagnosis),
    by = c("cleaned_id" = "matched_id")
  ) %>%
  select(original_id, micro_diagnosis)

# Step 2: Generate summaries
kmeans_summary_tables <- list()

for (k in names(kmeans_results)) {
  cluster_assignments <- kmeans_results[[k]]$cluster
  
  df <- data.frame(
    patient_id = names(cluster_assignments),
    cluster = paste0("C", cluster_assignments),
    stringsAsFactors = FALSE
  )
  
  # Join with diagnosis
  df <- left_join(df, diagnosis_lookup, by = c("patient_id" = "original_id"))
  
  # Filter out small clusters
  cluster_sizes <- df %>% count(cluster, name = "Patient_Count")
  valid_clusters <- cluster_sizes %>% filter(Patient_Count >= 20) %>% pull(cluster)
  df <- df %>% filter(cluster %in% valid_clusters)
  
  # Expand "Bacterial & Viral"
  df_expanded <- df %>%
    filter(micro_diagnosis == "Bacterial & Viral") %>%
    mutate(micro_diagnosis = "Bacterial") %>%
    bind_rows(
      df %>% filter(micro_diagnosis == "Bacterial & Viral") %>%
        mutate(micro_diagnosis = "Viral")
    )
  
  df_simple <- df %>% filter(micro_diagnosis != "Bacterial & Viral")
  df_final <- bind_rows(df_simple, df_expanded)
  
  # Create summary
  summary_df <- df_final %>%
    group_by(cluster, micro_diagnosis) %>%
    summarise(count = n(), .groups = "drop") %>%
    pivot_wider(names_from = micro_diagnosis, values_from = count, values_fill = 0)
  
  # Add total patient counts
  patient_counts <- df %>% count(cluster, name = "Patient_Count")
  
  summary_df <- left_join(summary_df, patient_counts, by = "cluster") %>%
    mutate(K = as.integer(k)) %>%
    select(K, cluster, Patient_Count, everything()) %>%
    arrange(K, desc(Patient_Count))
  
  kmeans_summary_tables[[k]] <- summary_df
}

# Combine all K cluster summaries
combined_kmeans_summary <- bind_rows(kmeans_summary_tables)

# To view table, please uncomment line below
print(combined_kmeans_summary, n=81)

################################################################################

# Plotting code 

###############################################################################

#----------------------- Retention rate --------------------------------------

# Compute % of retained patients in large clusters for each K
retention_df <- combined_kmeans_summary %>%
  group_by(K) %>%
  summarise(Patients_Captured = sum(Patient_Count)) %>%
  mutate(Percentage_Captured = round(100 * Patients_Captured / length(matched_patient_ids), 2))

# Plot
retention_plot <- ggplot(retention_df, aes(x = K, y = Percentage_Captured)) +
  geom_line(color = "darkgreen", size = 1.2) +
  geom_point(color = "black", size = 2) +
  labs(
    title = "Retention in Large Clusters (≥20 patients)",
    x = "K",
    y = "% of Patients Captured"
  ) +
  theme_minimal()

print(retention_plot)

#-------------------- Kappa Concordance (Viral) -----------------------------
kappa_results_kmeans <- data.frame()

for (k in names(kmeans_results)) {
  df <- combined_kmeans_summary %>% filter(K == as.integer(k))
  cluster_assignments <- kmeans_results[[k]]$cluster
  
  assignment_df <- data.frame(
    patient_id = names(cluster_assignments),
    cluster = paste0("C", cluster_assignments),
    stringsAsFactors = FALSE
  )
  
  assignment_df <- left_join(assignment_df, diagnosis_lookup, by = c("patient_id" = "original_id"))
  
  # Filter for large clusters only
  cluster_sizes <- assignment_df %>% count(cluster)
  valid_clusters <- cluster_sizes %>% filter(n >= 20) %>% pull(cluster)
  assignment_df <- assignment_df %>% filter(cluster %in% valid_clusters)
  
  if (nrow(assignment_df) < 20) next
  
  # Identify Viral cluster (with most viral patients)
  assignment_df <- assignment_df %>%
    mutate(is_viral = grepl("Viral", micro_diagnosis, ignore.case = TRUE))
  
  viral_counts <- assignment_df %>%
    group_by(cluster) %>%
    summarise(n_viral = sum(is_viral), .groups = "drop")
  
  viral_cluster <- viral_counts %>%
    slice_max(n_viral, with_ties = FALSE) %>%
    pull(cluster)
  
  assignment_df <- assignment_df %>%
    mutate(
      predicted_label = ifelse(cluster == viral_cluster, "Viral", "Non-Viral"),
      true_label = ifelse(is_viral, "Viral", "Non-Viral")
    )
  
  # Calculate Kappa
  if (length(unique(assignment_df$predicted_label)) > 1 &&
      length(unique(assignment_df$true_label)) > 1) {
    kappa_score <- kappa2(assignment_df[, c("predicted_label", "true_label")])$value
  } else {
    kappa_score <- NA
  }
  
  kappa_results_kmeans <- rbind(kappa_results_kmeans, data.frame(
    K = as.integer(k),
    Kappa = kappa_score
  ))
}

# Plot
kappa_plot <- ggplot(kappa_results_kmeans, aes(x = K, y = Kappa)) +
  geom_line(color = "purple", size = 1.2) +
  geom_point(size = 2, color = "black") +
  labs(
    title = "Cohen's Kappa Concordance (Viral vs. Non-Viral)",
    x = "K",
    y = "Kappa Statistic"
  ) +
  theme_minimal()

print(kappa_plot)

#------------------------ Weighted Purity -------------------------------------
# Calculate weighted purity for each K
wp_df_kmeans <- combined_kmeans_summary %>%
  pivot_longer(cols = -c(K, cluster, Patient_Count), names_to = "Diagnosis", values_to = "Count") %>%
  group_by(K, cluster) %>%
  mutate(Purity = Count / sum(Count)) %>%
  slice_max(order_by = Purity, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  group_by(K) %>%
  summarise(Weighted_Purity = sum(Purity * Patient_Count) / sum(Patient_Count), .groups = "drop")

# Plot
wp_plot <- ggplot(wp_df_kmeans, aes(x = K, y = Weighted_Purity)) +
  geom_line(color = "darkorange", size = 1.2) +
  geom_point(size = 2, color = "black") +
  labs(
    title = "Weighted Purity of K-means Clusters",
    x = "K",
    y = "Weighted Purity"
  ) +
  theme_minimal()

print(wp_plot)

##############################################################################

# Optimum K-value 

#############################################################################

eval_kmeans_df <- kappa_results_kmeans %>%
  inner_join(wp_df_kmeans, by = "K") %>%
  inner_join(retention_df, by = "K") %>%
  filter(Percentage_Captured >= 80)  # Optional filter for large retention

# Select configuration with highest Kappa, breaking ties with WP
max_kappa <- max(eval_kmeans_df$Kappa, na.rm = TRUE)
best_k <- eval_kmeans_df %>%
  filter(Kappa >= max_kappa - 0.05) %>%
  slice_max(order_by = Weighted_Purity, n = 1)

cat("\n🎯 Best K-means Clustering Configuration:\n")
print(best_k)

################################################################################

#  Optimum K-cluster Information - Stored in final_cluster_df

###############################################################################
# -- Replace this with your actual optimum K value if you've identified it programmatically
optimum_k <- best_k$K[1]  # Assuming you used the earlier selection code

# -- Get the kmeans results for that K
opt_cluster_assignments <- kmeans_results[[as.character(optimum_k)]]$cluster

# -- Create a basic assignment data frame
opt_df <- data.frame(
  patient_id = names(opt_cluster_assignments),
  cluster = opt_cluster_assignments,
  stringsAsFactors = FALSE
)

# -- Join with diagnosis
opt_df <- left_join(opt_df, diagnosis_lookup, by = c("patient_id" = "original_id"))

# -- Filter out clusters with <20 patients
cluster_sizes <- opt_df %>% count(cluster, name = "Patient_Count")
valid_clusters <- cluster_sizes %>% filter(Patient_Count >= 20) %>% pull(cluster)

opt_df <- opt_df %>% filter(cluster %in% valid_clusters)

# -- Recode "Bacterial & Viral" to "Viral"
opt_df <- opt_df %>%
  mutate(micro_diagnosis = if_else(micro_diagnosis == "Bacterial & Viral", "Viral", micro_diagnosis))

# -- Optional: format cluster as label (e.g., "C1") or keep as number
# opt_df <- opt_df %>% mutate(cluster = paste0("C", cluster))

# -- Final result with desired columns
final_cluster_df <- opt_df %>% select(patient_id, cluster, micro_diagnosis)

# -- Print or inspect
# print(final_cluster_df, n=30)

