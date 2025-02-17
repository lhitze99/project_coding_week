#### set global options and load libraries ####

DATA_DIR = "/cellfile/cellnet/erkTorQtl/data"
setwd("/cellfile/cellnet/erkTorQtl/data/BXD")

library(dplyr)
library(tidyr)
library(ggplot2)
library(openxlsx)


#### Load data ####

# Load RNA-seq data
#rna_data_filtered <- readRDS("/cellfile/cellnet/erkTorQtl/data/processed/bxd_rna_data_filtered.rds")

# Load the data that only contains the genes that are not missing in more than 90 % of the samples
rna_data_filtered <- readRDS(paste0(DATA_DIR, "/BXD/processed_data/bxd_rna_data_filtered_more_90_per_samples.rds"))

# Load the data that only contains the genes that are not missing in more than 90 % of the samples and only female samples
#rna_data_filtered <- readRDS(paste0(DATA_DIR,"/BXD/processed_data/bxd_rna_data_filtered_female.rds"))

# Load the rna_data_filtered_female_corrected_age.rds data
#rna_data_filtered <- readRDS(paste0(DATA_DIR,"/BXD/corrected_data/bxd_rna_data_filtered_female_corrected_age.rds"))

# Load the rna_data_filtered_female_corrected_age_diet.rds data
#rna_data_filtered <- readRDS(paste0(DATA_DIR,"/BXD/corrected_data/bxd_rna_data_filtered_female_corrected_age_diet.rds"))

# Load the rna_data_filtered_female_corrected_age_diet_plus_interaction.rds data
#rna_data_filtered <- readRDS(paste0(DATA_DIR,"/BXD/corrected_data/bxd_rna_data_filtered_female_corrected_age_diet_plus_interaction.rds"))

# Load metadata
metadata <- readRDS(paste0(DATA_DIR,"/BXD/processed_data/bxd_metadata.rds"))

# Load marker lists 
mtor_positive_translation <- readRDS(paste0(DATA_DIR,"/marker_lists_mouse/mtor_markers_positive_mouse_without_Jade2.rds"))
mtor_negative_translation <- readRDS(paste0(DATA_DIR,"/marker_lists_mouse/mtor_markers_negative_mouse.rds"))
erk_positive_translation <- readRDS(paste0(DATA_DIR,"/marker_lists_mouse/erk_markers_positive_mouse_without_Txnrd1.rds"))
erk_negative_translation <- readRDS(paste0(DATA_DIR,"/marker_lists_mouse/erk_markers_negative_mouse_without_EIF4EBP1.rds"))


#### Score calculation ####

## for mTOR ##

# rank genes for each sample with higher expression getting higher rank
rank_matrix <- apply(rna_data_filtered, 2, function(x){
  # rank genes, na.last = keep to keep NAs in original position and not include them in ranking process
  ranks <- rank(x, na.last = "keep")
  # normalize ranks (NAs ignored)
  ranks <- ranks / sum(!is.na(ranks))
  return(ranks)
})

# Get positive and negative markers for each sample
mtor_pos_rank_matrix <- rank_matrix[unlist(mtor_positive_translation),]
mtor_neg_rank_matrix <- rank_matrix[unlist(mtor_negative_translation),]

# Calculate mean rank for positive and negative markers for each sample (NAs ignored)
mtor_pos_mean_rank <- colMeans(mtor_pos_rank_matrix, na.rm = TRUE)
mtor_neg_mean_rank <- colMeans(mtor_neg_rank_matrix, na.rm = TRUE)

# Calculate score: mean rank of pos markers - mean rank of negative markers
mtor_score_matrix <- mtor_pos_mean_rank - mtor_neg_mean_rank
# we have na values here because not every sample has an entry for the genes


## for ERK ##

# rank genes for each sample with higher expression getting higher rank
rank_matrix <- apply(rna_data_filtered, 2, function(x){
  # rank genes, na.last = keep to keep NAs in original position and not include them in ranking process
  ranks <- rank(x, na.last = "keep")
  # normalize ranks (NAs ignored)
  ranks <- ranks / sum(!is.na(ranks))
  return(ranks)
})

# Get positive and negative markers for each sample
erk_pos_rank_matrix <- rank_matrix[unlist(erk_positive_translation),]
erk_neg_rank_matrix <- rank_matrix[unlist(erk_negative_translation),]

# Calculate mean rank for positive and negative markers for each sample (NAs ignored)
erk_pos_mean_rank <- colMeans(erk_pos_rank_matrix, na.rm = TRUE)
erk_neg_mean_rank <- colMeans(erk_neg_rank_matrix, na.rm = TRUE)

# Calculate score: mean rank of pos markers - mean rank of negative markers
erk_score_matrix <- erk_pos_mean_rank - erk_neg_mean_rank
# we have na values here because not every sample has an entry for the genes




### Score plots across age ####
library(ggplot2)

### for mTor ###

# matrix to dataframe
mtor_score_df <- as.data.frame(mtor_score_matrix)

# Merge scores and metadata
mtor_score_df$OmicsEarTag <- rownames(mtor_score_df)
mtor_score_df_with_metadata <- merge(mtor_score_df, metadata, by.x = "OmicsEarTag", by.y = "OmicsEarTag", all.x = TRUE)

# Filter out rows with NA scores
mtor_score_df_with_metadata <- mtor_score_df_with_metadata[!is.na(mtor_score_df_with_metadata$mtor_score_matrix), ]
#saveRDS(mtor_score_df_with_metadata, paste0("score_data/bxd_mtor_score_df_with_metadata_corrected_age.rds"))
#saveRDS(mtor_score_df_with_metadata, paste0("score_data/bxd_mtor_score_df_with_metadata_corrected_age_diet.rds"))
#saveRDS(mtor_score_df_with_metadata, paste0("score_data/bxd_mtor_score_df_with_metadata_corrected_age_diet_plus_interaction.rds"))

# Fit a linear model: score as a function of Age
lm_model <- lm(mtor_score_matrix ~ Age, data = mtor_score_df_with_metadata)

# Get the summary of the model to extract the p-value
summary_lm <- summary(lm_model)

# Extract the p-value for the Age variable (second coefficient)
p_value <- summary_lm$coefficients["Age", "Pr(>|t|)"]

# Plot with p-value annotation
ggplot(mtor_score_df_with_metadata, aes(x = Age, y = mtor_score_matrix)) + 
  geom_point(size = 3, alpha = 0.7) +  # scatter plot of scores by age
  geom_smooth(method = "lm", se = FALSE) +  # add linear regression line
  theme_minimal() + 
  labs(title = "mTOR Scores by Age", x = "Age (Days)", y = "mTOR Score") +
  scale_x_continuous(
    breaks = seq(from = 150, 
                 to = max(mtor_score_df_with_metadata$Age, na.rm = TRUE), 
                 by = 50),
    limits = c(min(mtor_score_df_with_metadata$Age, na.rm = TRUE), 
               max(mtor_score_df_with_metadata$Age, na.rm = TRUE))
  ) +
  annotate("text", x = Inf, y = Inf, label = paste("P-value:", format(p_value, digits = 3)), 
           hjust = 1.05, vjust = 1, size = 5)
#ggsave("score_data/Plots/mtor_score_by_age_uncorrected_data.png", width = 8, height = 6, dpi = 100)
#ggsave("score_data/Plots/mtor_score_by_age_corrected_data.png", width = 8, height = 6, dpi = 100)


# Visualize each diet separately
ggplot(mtor_score_df_with_metadata, aes(x = Age, y = mtor_score_matrix)) + 
  geom_point(size = 2, alpha = 0.7) + 
  geom_smooth(method = "lm", se = FALSE) + 
  theme_minimal() + 
  labs(title = "mTOR Scores across Age separated by Diet", x = "Age (Days)", y = "mTOR Score") +
  facet_wrap(~ Diet) +  # Facet by Strain
  theme(
    plot.title = element_text(size = 20),      # Title font size and bold
    axis.title = element_text(size = 18),                    # Axis labels font size
    axis.text = element_text(size = 18),                     # Tick labels font size
    strip.text = element_text(size = 16)      # Facet labels font size and bold
  )
#ggsave("score_data/Plots/mtor_score_by_age_separated_by_diet_uncorrected_data.png", width = 8, height = 6, dpi = 100)
#ggsave("score_data/Plots/mtor_score_by_age_separated_by_diet_age_corrected_data.png", width = 8, height = 6, dpi = 100)

# Check p values

# Subset the data for the "chow" (CD) diet group
lm_model_cd <- lm(mtor_score_matrix ~ Age, data = mtor_score_df_with_metadata[mtor_score_df_with_metadata$Diet == "CD", ])
summary_cd <- summary(lm_model_cd)

# Extract the p-value for the "CD" group
p_value_cd <- summary_cd$coefficients["Age", "Pr(>|t|)"]
cat("P-value for CD (chow) diet group:", p_value_cd, "\n")

# Subset the data for the "high-fat" (HFD) diet group
lm_model_hfd <- lm(mtor_score_matrix ~ Age, data = mtor_score_df_with_metadata[mtor_score_df_with_metadata$Diet == "HF", ])
summary_hfd <- summary(lm_model_hfd)

# Extract the p-value for the "HFD" group
p_value_hfd <- summary_hfd$coefficients["Age", "Pr(>|t|)"]
cat("P-value for HFD (high-fat) diet group:", p_value_hfd, "\n")


### Score plots across diet ####

# distribution of scores across different diets
ggplot(mtor_score_df_with_metadata, aes(x = Diet, y = mtor_score_matrix, fill = Diet)) + 
  geom_boxplot(alpha = 0.6) + 
  geom_jitter(width = 0.2, size = 1.5, color = "darkgreen", alpha = 0.5) +
  theme_minimal() + 
  labs(title = "mTOR Scores by Diet", x = "Diet", y = "mTOR Score")


### Score plots across strains ####

# Filter strains that have more than 3 samples
mtor_score_df_with_metadata_filtered <- mtor_score_df_with_metadata %>%
  group_by(Strain) %>%
  filter(n() > 3) %>%
  ungroup()  # Ungroup to avoid issues with ggplot

# Distribution of scores across different strains ordered based on median Score
ggplot(mtor_score_df_with_metadata_filtered, aes(x = reorder(Strain, -mtor_score_matrix, FUN = median), y = mtor_score_matrix, fill = Strain)) + 
  geom_boxplot(alpha = 0.6) + 
  geom_jitter(width = 0.2, size = 1.5, color = "darkgreen", alpha = 0.5) +
  theme_minimal() + 
  labs(title = "Distribution of mTOR Scores by Strain (with more than 3 samples)", x = "Strain", y = "mTOR Score") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none")
#ggsave("score_data/Plots/mtor_score_by_strain_uncorrected_data.png", width = 8, height = 6, dpi = 100)


### Score plots by sex ####

# distribution of scores by sex
ggplot(mtor_score_df_with_metadata, aes(x = Sex, y = mtor_score_matrix, fill = Sex)) + 
  geom_boxplot(alpha = 0.6) + 
  geom_jitter(width = 0.2, size = 1.5, color = "darkgreen", alpha = 0.5) +
  theme_minimal() + 
  labs(title = "Distribution of mTOR Scores by Sex", x = "Sex", y = "mTOR Score")



## Differences between young and old within the strains with more than 3 samples weighted by number of quantified genes ##

# Calculate the number of non-zero genes for all samples in rna_data_filtered
num_non_zero_genes <- colSums(rna_data_filtered != 0, na.rm = TRUE)

# Create a vector of the OmicsEarTag to match with the column names of rna_data_filtered
sample_names <- mtor_score_df_with_metadata_filtered$OmicsEarTag

# Extract the non-zero gene counts only for the samples present in mtor_score_df_with_metadata_filtered
num_non_zero_genes_filtered <- num_non_zero_genes[sample_names]

# Add the non-zero gene counts as a new column in the filtered DataFrame
mtor_score_df_with_metadata_filtered$NumNonZeroGenes <- num_non_zero_genes_filtered

# Calculate the weighted mean score for each Age_Group, Strain, and Diet
mtor_diff_data <- mtor_score_df_with_metadata_filtered %>%
  group_by(Strain, Age_Group, Diet) %>%
  summarise(mean_mtor_score = weighted.mean(mtor_score_matrix, NumNonZeroGenes, na.rm = TRUE), .groups = "drop")

# Calculate the difference in score between Age Groups within each Strain and Diet
mtor_diff_data_wide <- mtor_diff_data %>%
  pivot_wider(names_from = Age_Group, values_from = mean_mtor_score) %>%
  mutate(Difference = `>= 450 days` - `< 450 days`)  # Calculate difference between old and young

# Reshape to get differences for both diets (CD and HF)
diet_diff <- mtor_diff_data_wide %>%
  filter(Diet %in% c("CD", "HF")) %>%
  dplyr::select(Strain, Diet, Difference) %>%
  pivot_wider(names_from = Diet, values_from = Difference)

# Plot differences between Age Groups for each Strain
ggplot(diet_diff, aes(x = CD, y = HF, label = Strain)) +
  geom_point(size = 3, color = "darkblue") +  # Points for each strain
  geom_text(vjust = 1.5, hjust = 1.5) +  # Add strain labels next to the points
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +  # Diagonal line (y = x)
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Horizontal line at y = 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # Vertical line at x = 0
  theme_minimal() +
  labs(title = "Difference in mTOR Scores Between Younger and Older Mice by Strain (with Weighted Mean Score)",
       subtitle = "X-axis: Difference in Low-Fat Diet (CD), Y-axis: Difference in High-Fat Diet (HF)",
       x = "Difference in mTOR Score (CD)",
       y = "Difference in mTOR Score (HF)")
#ggsave("score_data/Plots/mtor_score_difference_by_strain_weighted_uncorrected_data.png", width = 8, height = 6, dpi = 100)





### for ERK ###

# matrix to dataframe
erk_score_df <- as.data.frame(erk_score_matrix)

# Merge scores and metadata
erk_score_df$OmicsEarTag <- rownames(erk_score_df)
erk_score_df_with_metadata <- merge(erk_score_df, metadata, by.x = "OmicsEarTag", by.y = "OmicsEarTag", all.x = TRUE)

# Filter out rows with NA scores
erk_score_df_with_metadata <- erk_score_df_with_metadata[!is.na(erk_score_df_with_metadata$erk_score_matrix), ]
#saveRDS(erk_score_df_with_metadata, paste0("score_data/bxd_erk_score_df_with_metadata_corrected_age.rds"))
#saveRDS(erk_score_df_with_metadata, paste0("score_data/bxd_erk_score_df_with_metadata_corrected_age_diet.rds"))
#saveRDS(erk_score_df_with_metadata, paste0("score_data/bxd_erk_score_df_with_metadata_corrected_age_diet_plus_interaction.rds"))

# Fit a linear model: score as a function of Age
lm_model <- lm(erk_score_matrix ~ Age, data = erk_score_df_with_metadata)

# Get the summary of the model to extract the p-value
summary_lm <- summary(lm_model)

# Extract the p-value for the Age variable (second coefficient)
p_value <- summary_lm$coefficients["Age", "Pr(>|t|)"]

# Plot with p-value annotation
ggplot(erk_score_df_with_metadata, aes(x = Age, y = erk_score_matrix)) + 
  geom_point(size = 3, alpha = 0.7) +  # scatter plot of scores by age
  geom_smooth(method = "lm", se = FALSE) +  # add linear regression line
  theme_minimal() + 
  labs(title = "ERK Scores by Age", x = "Age (Days)", y = "ERK Score") +
  scale_x_continuous(
    breaks = seq(from = 150, 
                 to = max(erk_score_df_with_metadata$Age, na.rm = TRUE), 
                 by = 50),
    limits = c(min(erk_score_df_with_metadata$Age, na.rm = TRUE), 
               max(erk_score_df_with_metadata$Age, na.rm = TRUE))
  ) +
  annotate("text", x = Inf, y = Inf, label = paste("P-value:", format(p_value, digits = 3)), 
           hjust = 1.05, vjust = 2, size = 5)
#ggsave("score_data/Plots/erk_score_by_age_uncorrected_data.png", width = 8, height = 6, dpi = 100)
#ggsave("score_data/Plots/erk_score_by_age_corrected_data.png", width = 8, height = 6, dpi = 100)


# Visualize each diet separately
ggplot(erk_score_df_with_metadata, aes(x = Age, y = erk_score_matrix)) + 
  geom_point(size = 2, alpha = 0.7) + 
  geom_smooth(method = "lm", se = FALSE) + 
  theme_minimal() + 
  labs(title = "ERK Scores across Age separated by Diet", x = "Age (Days)", y = "ERK Score") +
  facet_wrap(~ Diet) +  # Facet by Strain
  theme(
    plot.title = element_text(size = 20),      # Title font size and bold
    axis.title = element_text(size = 18),                    # Axis labels font size
    axis.text = element_text(size = 18),                     # Tick labels font size
    strip.text = element_text(size = 16)      # Facet labels font size and bold
  )
#ggsave("score_data/Plots/erk_score_by_age_separated_by_diet_uncorrected_data.png", width = 8, height = 6, dpi = 100)
#ggsave("score_data/Plots/erk_score_by_age_separated_by_diet_age_corrected_data.png", width = 8, height = 6, dpi = 100)

# Check p values

# Subset the data for the "chow" (CD) diet group
lm_model_cd <- lm(erk_score_matrix ~ Age, data = erk_score_df_with_metadata[erk_score_df_with_metadata$Diet == "CD", ])
summary_cd <- summary(lm_model_cd)

# Extract the p-value for the "CD" group
p_value_cd <- summary_cd$coefficients["Age", "Pr(>|t|)"]
cat("P-value for CD (chow) diet group:", p_value_cd, "\n")

# Subset the data for the "high-fat" (HFD) diet group
lm_model_hfd <- lm(erk_score_matrix ~ Age, data = erk_score_df_with_metadata[erk_score_df_with_metadata$Diet == "HF", ])
summary_hfd <- summary(lm_model_hfd)

# Extract the p-value for the "HFD" group
p_value_hfd <- summary_hfd$coefficients["Age", "Pr(>|t|)"]
cat("P-value for HFD (high-fat) diet group:", p_value_hfd, "\n")


### Score plots across diet ####

# distribution of scores across different diets
ggplot(erk_score_df_with_metadata, aes(x = Diet, y = erk_score_matrix, fill = Diet)) + 
  geom_boxplot(alpha = 0.6) + 
  geom_jitter(width = 0.2, size = 1.5, color = "darkgreen", alpha = 0.5) +
  theme_minimal() + 
  labs(title = "ERK Scores by Diet", x = "Diet", y = "ERK Score")


### Score plots across strains ####

# Filter strains that have more than 3 samples
erk_score_df_with_metadata_filtered <- erk_score_df_with_metadata %>%
  group_by(Strain) %>%
  filter(n() > 3) %>%
  ungroup()  # Ungroup to avoid issues with ggplot

# distribution of scores across different strains ordered based on median Score
ggplot(erk_score_df_with_metadata_filtered, aes(x = reorder(Strain, -erk_score_matrix, FUN = median), y = erk_score_matrix, fill = Strain)) + 
  geom_boxplot(alpha = 0.6) + 
  geom_jitter(width = 0.2, size = 1.5, color = "darkgreen", alpha = 0.5) +
  theme_minimal() + 
  labs(title = "Distribution of ERK Scores by Strain (with more than 3 samples)", x = "Strain", y = "ERK Score") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none")
#ggsave("score_data/Plots/erk_score_by_strain_uncorrected_data.png", width = 8, height = 6, dpi = 100)


### Score plots by sex ####

# distribution of scores by sex
ggplot(erk_score_df_with_metadata, aes(x = Sex, y = erk_score_matrix, fill = Sex)) + 
  geom_boxplot(alpha = 0.6) + 
  geom_jitter(width = 0.2, size = 1.5, color = "darkgreen", alpha = 0.5) +
  theme_minimal() + 
  labs(title = "Distribution of ERK Scores by Sex", x = "Sex", y = "ERK Score")



## Differences between young and old within the strains with more than 3 samples weighted by number of quantified genes ##

# Calculate the number of non-zero genes for all samples in rna_data_filtered
num_non_zero_genes <- colSums(rna_data_filtered != 0, na.rm = TRUE)

# Create a vector of the OmicsEarTag to match with the column names of rna_data_filtered
sample_names <- erk_score_df_with_metadata_filtered$OmicsEarTag

# Extract the non-zero gene counts only for the samples present in erk_score_df_with_metadata_filtered
num_non_zero_genes_filtered <- num_non_zero_genes[sample_names]

# Add the non-zero gene counts as a new column in the filtered DataFrame
erk_score_df_with_metadata_filtered$NumNonZeroGenes <- num_non_zero_genes_filtered

# Calculate the weighted mean score for each Age_Group, Strain, and Diet
erk_diff_data <- erk_score_df_with_metadata_filtered %>%
  group_by(Strain, Age_Group, Diet) %>%
  summarise(mean_erk_score = weighted.mean(erk_score_matrix, NumNonZeroGenes, na.rm = TRUE), .groups = "drop")

# Calculate the difference in score between Age Groups within each Strain and Diet
erk_diff_data_wide <- erk_diff_data %>%
  pivot_wider(names_from = Age_Group, values_from = mean_erk_score) %>%
  mutate(Difference = `>= 450 days` - `< 450 days`)  # Calculate difference between old and young

# Reshape to get differences for both diets (CD and HF)
diet_diff <- erk_diff_data_wide %>%
  filter(Diet %in% c("CD", "HF")) %>%
  dplyr::select(Strain, Diet, Difference) %>%
  pivot_wider(names_from = Diet, values_from = Difference)

# Plot differences between Age Groups for each Strain
ggplot(diet_diff, aes(x = CD, y = HF, label = Strain)) +
  geom_point(size = 3, color = "darkblue") +  # Points for each strain
  geom_text(vjust = 1.5, hjust = 1.5) +  # Add strain labels next to the points
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +  # Diagonal line (y = x)
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Horizontal line at y = 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # Vertical line at x = 0
  theme_minimal() +
  labs(title = "Difference in ERK Scores Between Younger and Older Mice by Strain (with Weighted Mean Score)",
       subtitle = "X-axis: Difference in Low-Fat Diet (CD), Y-axis: Difference in High-Fat Diet (HF)",
       x = "Difference in ERK Score (CD)",
       y = "Difference in ERK Score (HF)")
#ggsave("score_data/Plots/erk_score_difference_by_strain_weighted_uncorrected_data.png", width = 8, height = 6, dpi = 100)
