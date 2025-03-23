##############################
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
gene_long <- gene_long %>%
  group_by(Gene) %>%
  mutate(Expression_z = (Expression - mean(Expression[Group == "Controls"], na.rm = TRUE)) /
           sd(Expression[Group == "Controls"], na.rm = TRUE)) %>%
  ungroup()

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
