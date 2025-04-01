####################################################################################

# BOX PLOTS

####################################################################################
# Step 1: Select the relevant genes and reshape to long format
##############################
# Subset gene_data to include only IFI27, MX1, and IFI44L
gene_subset <- gene_data[, c("IFI27", "MX1", "IFI44L")]

# Add a SampleID column from row names (without altering gene_data)
gene_subset$SampleID <- rownames(gene_data)

# Convert from wide to long format
library(tidyr)
gene_long <- pivot_longer(gene_subset,
                          cols = c("IFI27", "MX1", "IFI44L"),
                          names_to = "Gene",
                          values_to = "Expression")

##############################
# Step 2: Create a temporary Group variable on the fly
##############################
# If the SampleID starts with "WH0", label as "Controls"; otherwise, label as "BioAID"
gene_long <- gene_long %>%
  mutate(Group = if_else(grepl("^WH0", SampleID), "Controls", "BioAID"))

##############################
# Step 3: Compute the z‑scores for each gene using only Controls as the reference
##############################
# Function to filter Controls and compute z-scores per gene
filter_and_zscore <- function(df) {
  df %>%
    group_by(Gene) %>%
    group_modify(~ {
      gene_df <- .
      control_df <- gene_df %>% filter(Group == "Controls")
      
      # Calculate IQR and threshold
      Q1 <- quantile(control_df$Expression, 0.25, na.rm = TRUE)
      Q3 <- quantile(control_df$Expression, 0.75, na.rm = TRUE)
      IQR_val <- Q3 - Q1
      upper_thresh <- Q3 + 1.5 * IQR_val
      
      # Filter non-outlier controls
      filtered_controls <- control_df %>% filter(Expression <= upper_thresh)
      
      # Compute mean and sd from filtered Controls
      mean_val <- mean(filtered_controls$Expression, na.rm = TRUE)
      sd_val <- sd(filtered_controls$Expression, na.rm = TRUE)
      
      # Compute z-score for everyone in this gene group
      gene_df <- gene_df %>%
        mutate(Expression_z = (Expression - mean_val) / sd_val)
      
      # Exclude outlier Controls from the output
      gene_df %>%
        filter(!(Group == "Controls" & Expression > upper_thresh))
    }) %>%
    ungroup()
}

# Apply the filtering and z-score calculation
gene_long <- filter_and_zscore(gene_long)
##############################
# Step 4: Create a facetted box plot for IFI27, MX1, and IFI44L
##############################
box_plot <- ggplot(gene_long, aes(x = Group, y = Expression_z, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +   # Box plots without drawing outliers
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +   # Add individual sample points
  facet_wrap(~ Gene, scales = "free_y") +            # Facet by gene
  scale_fill_manual(values = c("Controls" = "royalblue", "BioAID" = "tomato")) +
  labs(
    title = "Expression Z-Scores: IFI27, MX1, and IFI44L",
    x = "Group",
    y = "Biomarker Z Score"
  ) +
  theme_pubr() +
  theme(legend.position = "none")

print(box_plot)

################################################################################
# Testing the hypothesis
###############################################################################

# --- Assume the gene_long data frame is already created as in your earlier code ---
# gene_long contains: SampleID, Gene, Expression, Expression_z and a Group column 
# (created as: if_else(grepl("^WH0", SampleID), "Controls", "BioAID"))

# Filter: For BioAID, keep only patients with clinical data (using matched_patient_ids)
gene_long_matched <- gene_long %>%
  filter(Group == "Controls" | (Group == "BioAID" & SampleID %in% matched_patient_ids)) %>%
  # For BioAID patients, join the micro_diagnosis from plot_info
  left_join(plot_info %>% select(original_id, micro_diagnosis),
            by = c("SampleID" = "original_id")) %>%
  # Create a new column to use for color mapping:
  # Controls remain "Controls" and BioAID patients get labeled with their micro_diagnosis.
  mutate(fill_group = if_else(Group == "Controls", "Controls", micro_diagnosis))

# Define the new color mapping:
# Controls remain a fixed color, while BioAID patients are colored by their micro_diagnosis
color_mapping_new <- c("Controls" = "royalblue",
                       "None" = "blue", 
                       "Bacterial" = "red", 
                       "Viral" = "green", 
                       "Bacterial & Viral" = "black")

# Now, create the box plot.
# We plot Expression_z by Group (Controls vs. BioAID), and overlay jittered points that are 
# colored according to the new fill_group variable.
box_plot_matched <- ggplot(gene_long_matched, aes(x = Group, y = Expression_z)) +
  geom_boxplot(aes(fill = Group), outlier.shape = NA, alpha = 0.7) +  # aggregated box plot per group
  geom_jitter(aes(color = fill_group), width = 0.2, size = 2, alpha = 0.8) +  # individual points
  facet_wrap(~ Gene, scales = "free_y") +  # one facet per gene
  scale_fill_manual(values = c("Controls" = "royalblue", "BioAID" = "tomato")) +
  scale_color_manual(values = color_mapping_new) +
  labs(title = "Expression Z-Scores for ISGs by Group and micro_diagnosis",
       x = "Group", y = "Biomarker Z Score") +
  theme_pubr() +
  theme(legend.position = "right")

# Display the plot
print(box_plot_matched)

####################################################################################

# HISTOGRAM AND DENSITY PLOTS

####################################################################################

# First, filter the micro-diagnosis data to remove the "Bacterial & Viral" group
gene_long_micro_filtered <- gene_long_matched %>%
  filter(fill_group != "Bacterial & Viral")

# Update the color mapping to remove the "Bacterial & Viral" key
color_mapping_micro <- c("Controls" = "royalblue",
                         "None" = "blue", 
                         "Bacterial" = "red", 
                         "Viral" = "green")

# ---------------------------
# Micro_diagnosis Plots (with "Bacterial & Viral" removed)
# ---------------------------

# Histogram: micro_diagnosis (with Controls unchanged)
histogram_micro <- ggplot(gene_long_micro_filtered, aes(x = Expression_z, fill = fill_group)) +
  geom_histogram(position = "dodge", bins = 15, color = "black", alpha = 0.7) +
  facet_wrap(~ Gene, scales = "free_y") +
  scale_fill_manual(values = color_mapping_micro) +
  labs(title = "Histogram of ISG Expression Z-Scores by micro_diagnosis",
       x = "Expression Z-Score",
       y = "Count") +
  theme_pubr()

print(histogram_micro)

# Density Plot: micro_diagnosis (with Controls unchanged)
density_micro <- ggplot(gene_long_micro_filtered, aes(x = Expression_z, fill = fill_group)) +
  geom_density(alpha = 0.7, color = "black") +
  facet_wrap(~ Gene, scales = "free_y") +
  scale_fill_manual(values = color_mapping_micro) +
  labs(title = "Density Plot of ISG Expression Z-Scores by micro_diagnosis",
       x = "Expression Z-Score",
       y = "Density") +
  theme_pubr()

print(density_micro)

#########################################################################

# AUROC for each ISG in differentiating between Viral vs Non-Viral cases

#########################################################################
# Prepare data: binary Viral label
roc_data <- gene_long_matched %>%
  mutate(
    Viral_binary = if_else(micro_diagnosis == "Viral", 1, 0)
  ) %>%
  filter(Gene %in% c("IFI27", "MX1", "IFI44L"))

# Create ROC objects and extract data for plotting
roc_list <- list()
roc_plot_data <- data.frame()

for (gene in c("IFI27", "MX1", "IFI44L")) {
  gene_df <- roc_data %>% filter(Gene == gene)
  roc_obj <- roc(gene_df$Viral_binary, gene_df$Expression_z, direction = "<", quiet = TRUE)
  
  # Store the AUC
  auc_val <- round(auc(roc_obj), 3)
  
  # Extract plot data from roc object
  coords_df <- data.frame(
    FPR = 1 - roc_obj$specificities,
    TPR = roc_obj$sensitivities,
    Gene = paste0(gene, " (AUC = ", auc_val, ")")
  )
  
  roc_plot_data <- rbind(roc_plot_data, coords_df)
}

# Plot all ROC curves together
roc_plot <- ggplot(roc_plot_data, aes(x = FPR, y = TPR, color = Gene)) +
  geom_line(size = 1.2) +
  geom_abline(linetype = "dashed", color = "gray50") +
  theme_pubr() +
  labs(
    title = "ROC Curves for IFI27, MX1, and IFI44L",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)",
    color = "Gene"
  ) +
  theme(legend.position = "bottom")

print(roc_plot)
