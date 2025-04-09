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

# Limma 

################################################################################

# Threshold
adj_p_threshold <- 0.05

# Store DEG counts and lists
deg_counts <- data.frame()
deg_details <- list()

for (i in seq_along(contrast_names)) {
  contrast_name <- contrast_names[i]
  
  # Get top table
  results_table <- topTable(fit_ebayes, coef = i, number = Inf, sort.by = "P") %>%
    rownames_to_column("Gene") %>%
    filter(adj.P.Val < adj_p_threshold)
  
  # Upregulated: logFC > 0
  upregulated <- results_table %>%
    filter(logFC > 1) %>%
    arrange(adj.P.Val) %>%
    pull(Gene)
  
  # Downregulated: logFC < 0
  downregulated <- results_table %>%
    filter(logFC < -1) %>%
    arrange(adj.P.Val) %>%
    pull(Gene)
  
  # Store DEG counts
  deg_counts <- rbind(deg_counts, data.frame(
    Contrast = contrast_name,
    Direction = "Upregulated",
    Count = length(upregulated)
  ))
  
  deg_counts <- rbind(deg_counts, data.frame(
    Contrast = contrast_name,
    Direction = "Downregulated",
    Count = length(downregulated)
  ))
  
  # Store DEG lists
  deg_details[[contrast_name]] <- list(
    upregulated = upregulated,
    downregulated = downregulated
  )
  
  cat("Contrast:", contrast_name, "- Up:", length(upregulated), "Down:", length(downregulated), "\n")
}

# --- Bar Plot of DEG Counts ---
ggplot(deg_counts, aes(x = Contrast, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  labs(
    title = "Differentially Expressed Genes by Contrast",
    y = "Number of DEGs", x = "Contrast"
  ) +
  theme_minimal(base_size = 12) +
  scale_fill_manual(values = c("Upregulated" = "steelblue", "Downregulated" = "firebrick")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

############################### Lists of DEGs ##################################

# deg_details$C1vsRest$upregulated - Genes upregulated in C1 compared to C2, C3, and C4
# deg_details$C1vsRest$downregulated - Genes downregulated in C1 compared to the rest

# deg_details$C2vsRest$upregulated - Genes upregulated in C2 compared to C1, C3, and C4
# deg_details$C2vsRest$downregulated - Genes downregulated in C2 compared to the rest

# deg_details$C3vsRest$upregulated - Genes upregulated in C3 compared to C1, C2, and C4
# deg_details$C3vsRest$downregulated - Genes downregulated in C3 compared to the rest

# deg_details$C4vsRest$upregulated - Genes upregulated in C4 compared to C1, C2, and C3
# deg_details$C4vsRest$downregulated - Genes downregulated in C4 compared to the rest
