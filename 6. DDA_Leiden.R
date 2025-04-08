optimum_cluster <- read.csv("C:/Users/ABSin/Downloads/4. OptimumLeidenCluster.csv")


###############################################################################

# Sensitivity & Specificity

##############################################################################

# Step 1: Identify the community with the most viral patients
viral_counts <- optimum_cluster %>%
  filter(micro_diagnosis == "Viral") %>%
  count(community, name = "viral_count")

dominant_viral_comm <- viral_counts %>%
  arrange(desc(viral_count)) %>%
  slice(1) %>%
  pull(community)

cat("Dominant Viral Community:", dominant_viral_comm, "\n")

# Step 2: Create predicted and true labels
optimum_cluster <- optimum_cluster %>%
  mutate(
    predicted_label = if_else(community == dominant_viral_comm, "Viral", "Non-Viral"),
    true_label = if_else(micro_diagnosis == "Viral", "Viral", "Non-Viral")
  )

# Step 3: Build confusion matrix
conf_matrix <- table(Predicted = optimum_cluster$predicted_label,
                     Actual = optimum_cluster$true_label)

print(conf_matrix)

# Step 4: Calculate Sensitivity and Specificity
TP <- conf_matrix["Viral", "Viral"]
FP <- conf_matrix["Viral", "Non-Viral"]
TN <- conf_matrix["Non-Viral", "Non-Viral"]
FN <- conf_matrix["Non-Viral", "Viral"]

sensitivity <- TP / (TP + FN)
specificity <- TN / (TN + FP)

cat("Sensitivity:", round(sensitivity, 3), "\n")
cat("Specificity:", round(specificity, 3), "\n")

################################################################################

# Top 50 Discriminating Genes Analysis for Optimum Leiden Clusters

###############################################################################

# Ensure matched_patients is in memory and rownames are patient IDs
# matched_patients <- your scaled log2 TPM matrix [already in memory]

# Match sample order between expression matrix and cluster labels
common_ids <- intersect(rownames(matched_patients), optimum_cluster$patient_id)
matched_expr <- matched_patients[common_ids, ]
cluster_info <- optimum_cluster %>%
  filter(patient_id %in% common_ids) %>%
  arrange(match(patient_id, rownames(matched_expr)))

# Sanity check
stopifnot(all(cluster_info$patient_id == rownames(matched_expr)))

# Convert community to factor
cluster_factor <- factor(cluster_info$community)

# Compute p-values via one-way ANOVA
pvals <- apply(matched_expr, 2, function(gene_expr) {
  fit <- aov(gene_expr ~ cluster_factor)
  summary(fit)[[1]][["Pr(>F)"]][1]
})

# Adjust p-values and select top 50
pvals_adj <- p.adjust(pvals, method = "BH")
top_genes <- names(sort(pvals_adj))[1:50]

cat("Top 50 genes:\n", paste(top_genes, collapse = ", "), "\n")

top_expr <- matched_expr[, top_genes]
top_expr_z <- t(scale(t(top_expr)))  # row-wise Z-score

annotation_row <- data.frame(Community = cluster_factor)
rownames(annotation_row) <- rownames(top_expr_z)

pheatmap(
  top_expr_z,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_row = annotation_row,
  show_rownames = FALSE,
  main = "Top 50 Discriminating Genes Across Leiden Clusters"
)

###############################################################################

# Limma 

##############################################################################
# Step 1: Align clusters and expression matrix
common_ids <- intersect(rownames(matched_patients), optimum_cluster$patient_id)

matched_expr <- matched_patients[common_ids, ]
cluster_info <- optimum_cluster %>%
  filter(patient_id %in% common_ids) %>%
  arrange(match(patient_id, rownames(matched_expr)))

stopifnot(all(cluster_info$patient_id == rownames(matched_expr)))
cluster_factor <- factor(cluster_info$community)

# Step 2: Build design matrix (no intercept)
design <- model.matrix(~ 0 + cluster_factor)
colnames(design) <- paste0("Cluster", levels(cluster_factor))

# Step 3: Fit the model
fit <- lmFit(t(matched_expr), design)

# Step 4: Define pairwise contrasts
contrast_matrix <- makeContrasts(
  Cluster1vs2 = Cluster1 - Cluster2,
  Cluster1vs3 = Cluster1 - Cluster3,
  Cluster1vs4 = Cluster1 - Cluster4,
  Cluster2vs3 = Cluster2 - Cluster3,
  Cluster2vs4 = Cluster2 - Cluster4,
  Cluster3vs4 = Cluster3 - Cluster4,
  levels = design
)

# Step 5: Apply contrasts and compute statistics
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Step 6: Extract DEGs for each contrast (FDR + logFC filtering)
fdr_threshold <- 0.05
logfc_threshold <- 1

results_list <- lapply(colnames(contrast_matrix), function(contrast_name) {
  res <- topTable(fit2, coef = contrast_name, adjust = "BH", number = Inf)
  res <- res %>% filter(adj.P.Val < fdr_threshold & abs(logFC) > logfc_threshold)
  res$contrast <- contrast_name
  res$Gene <- rownames(res)
  return(res)
})

# Step 7: Combine all results
deg_results_filtered <- bind_rows(results_list)

# Step 8: Count DEGs per contrast (for a basic bar plot)
deg_counts_filtered <- deg_results_filtered %>%
  group_by(contrast) %>%
  summarise(FDR_filtered_DEGs = n())
print(deg_counts_filtered)

# You might keep your original bar plot if needed:
ggplot(deg_counts_filtered, aes(x = contrast, y = FDR_filtered_DEGs)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Number of DEGs for each Pairwise Comparison",
       x = "Contrast",
       y = "Number of DEGs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


###############################################################################
# UpSet Plot with Inclusive vs. Exclusive Intersection Counts
###############################################################################

# Check if there are any DEGs to plot
if (nrow(deg_results_filtered) == 0) {
  warning("No DEGs found in 'deg_results_filtered'. Cannot create UpSet plot.")
} else {
  
  # --- Step 1: Create the list for UpSet analysis ---
  # Group by contrast: for each contrast, we get the unique gene names.
  deg_list_for_upset <- lapply(
    split(deg_results_filtered$Gene, deg_results_filtered$contrast),
    unique
  )
  # Force a stable, reproducible ordering (alphabetical order of contrast names)
  deg_list_for_upset <- deg_list_for_upset[sort(names(deg_list_for_upset))]
  
  # Convert the list to a binary membership matrix (rows = genes, columns = contrasts)
  upset_data <- fromList(deg_list_for_upset)
  
  # --- Step 2: Compute inclusive & exclusive intersection sizes ---
  # Define a function that, for each combination of contrasts, computes:
  # - inclusive_size: genes present in all of those contrasts (even if in additional sets)
  # - exclusive_size: genes present exactly in those contrasts (and no others)
  calculate_inclusive_exclusive <- function(upset_mat) {
    set_names <- colnames(upset_mat)
    
    all_combos <- unlist(
      lapply(seq_along(set_names), function(k)
        utils::combn(set_names, k, simplify = FALSE)
      ),
      recursive = FALSE
    )
    
    results <- lapply(all_combos, function(sets) {
      # Inclusive: count genes that are present in every one of the selected sets.
      inclusive_count <- sum(rowSums(upset_mat[, sets, drop = FALSE]) == length(sets))
      # Exclusive: count genes that are present in the selected sets and in no others.
      exclusive_count <- sum(
        rowSums(upset_mat[, sets, drop = FALSE]) == length(sets) &
          rowSums(upset_mat) == length(sets)
      )
      
      data.frame(
        combination    = paste(sets, collapse = " & "),
        sets_count     = length(sets),
        inclusive_size = inclusive_count,
        exclusive_size = exclusive_count,
        stringsAsFactors = FALSE
      )
    })
    
    results_df <- do.call(rbind, results)
    # Sort stably: first by descending inclusive_size, then lexicographically.
    results_df <- results_df %>% arrange(desc(inclusive_size), combination)
    return(results_df)
  }
  
  inc_exc_df <- calculate_inclusive_exclusive(upset_data)
  
  # --- Step 3: Select a reproducible subset (e.g., top 40 intersections) ---
  # (Adjust 'n' if needed.)
  top_inc_exc_df <- inc_exc_df %>% 
    slice_max(inclusive_size, n = 40) %>% 
    mutate(comb_id = factor(row_number(), levels = 1:n()))
  
  # Reshape the data for a grouped bar plot (one row per combination per intersection type)
  df_melted <- top_inc_exc_df %>%
    pivot_longer(
      cols = c("inclusive_size", "exclusive_size"),
      names_to = "intersection_type",
      values_to = "size"
    ) %>%
    mutate(intersection_type = factor(intersection_type,
                                      levels = c("inclusive_size", "exclusive_size"),
                                      labels = c("Inclusive", "Exclusive")))
  
  # --- Step 4: Create the Bar Plot (with no x-axis text labels) ---
  # Here we assign custom colors (e.g. Inclusive in orange, Exclusive in steelblue)
  bar_plot <- ggplot(df_melted, aes(x = comb_id, y = size, fill = intersection_type)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
    scale_fill_manual(values = c("Inclusive" = "orange", "Exclusive" = "steelblue")) +
    theme_minimal(base_size = 13) +
    labs(y = "Number of DEGs", fill = "Intersection Type",
         title = "Overlap of Differentially Expressed Genes Across Contrasts") +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),  # Remove x-axis text labels
      axis.ticks.x = element_blank(),
      plot.margin = margin(b = 0)
    )
  
  # --- Step 5: Create the Membership (Dot Matrix) Panel ---
  # This panel mimics the lower part of a traditional UpSet plot.
  all_sets <- colnames(upset_data)
  
  # Create a data frame mapping each combination to each contrast
  membership_df <- expand.grid(
    comb_id = 1:nrow(top_inc_exc_df),
    contrast = all_sets,
    stringsAsFactors = FALSE
  )
  # Retrieve the combination string for each comb_id
  membership_df$combination <- top_inc_exc_df$combination[membership_df$comb_id]
  
  # Mark membership: for each row, a 1 if the contrast is present in the combination, else 0.
  membership_df$membership <- mapply(function(comb, cont) {
    # Use word boundaries to avoid partial matching
    if (grepl(paste0("\\b", cont, "\\b"), comb)) 1 else 0
  }, membership_df$combination, membership_df$contrast)
  
  # Build the dot matrix plot (reversing the order of contrasts as in a standard UpSet plot)
  membership_plot <- ggplot() +
    geom_point(data = membership_df %>% filter(membership == 1),
               aes(x = factor(comb_id), y = factor(contrast, levels = rev(all_sets))),
               size = 4) +
    geom_line(data = membership_df %>% filter(membership == 1) %>%
                arrange(comb_id, factor(contrast, levels = rev(all_sets))),
              aes(x = as.numeric(factor(comb_id)),
                  y = as.numeric(factor(contrast, levels = rev(all_sets))),
                  group = comb_id)) +
    scale_y_discrete(name = "Contrasts") +
    xlab("") +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(t = 0)
    )
  
  # --- Step 6: Combine the Bar Plot and Membership Plot ---
  # Using the patchwork package for a stacked layout.
  final_plot <- bar_plot / membership_plot + plot_layout(heights = c(3, 1))
  
  print(final_plot)
  
  cat("\nUpSet-style plot (with both Inclusive and Exclusive intersections) generated.\n")
  
} # End of check for non-empty deg_results_filtered
