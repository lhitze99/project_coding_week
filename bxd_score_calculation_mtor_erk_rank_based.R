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

# Load Lifespan Information Data (treat empty cells as NA values)
survive_raw <- read.xlsx(paste0(DATA_DIR,"/BXD/raw_data/aTablesS1_AgingDB.xlsx", na.strings=""))
survive <- subset(survive_raw, survive_raw$CauseOfDeath=="Natural" | survive_raw$CauseOfDeath=="Euthanized") # & survive_raw$Sex=="F"
# Subset data to include only strains with at least 2 samples
survive <- survive %>%
  group_by(StrainNameCurrent) %>%
  filter(n() >= 2) %>%
  ungroup()  # Remove the grouping structure


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


### Score plots across strains for each diet separately - only necessary for age corrected data (for QTL) ####

# Filter strains that have more than 3 samples and that are on CD diet
mtor_score_df_with_metadata_filtered_cd <- mtor_score_df_with_metadata %>%
  group_by(Strain) %>%
  filter(n() > 3 & Diet == "CD") %>%
  ungroup()  # Ungroup to avoid issues with ggplot

# Convert Strain to character to avoid issues with colors
mtor_score_df_with_metadata_filtered_cd$Strain <- as.character(mtor_score_df_with_metadata_filtered_cd$Strain)

# Color all strains grey and the parents in different colors
mtor_score_df_with_metadata_filtered_cd <- mtor_score_df_with_metadata_filtered_cd %>%
  mutate(Color = ifelse(Strain %in% c("C57BL6J", "DBA2J"), Strain, "Other"))

# distribution of scores across different strains ordered based on median Score for CD diet
ggplot(mtor_score_df_with_metadata_filtered_cd, aes(x = reorder(Strain, -mtor_score_matrix, FUN = median), y = mtor_score_matrix, fill = Color)) + 
  geom_boxplot(alpha = 0.6) + 
  geom_jitter(width = 0.2, size = 1.5, color = "darkgreen", alpha = 0.5) +
  theme_minimal() + 
  ylim(-0.35, 0.35) +
  labs(title = "CD", x = "Strain", y = "mTOR Score") + 
  scale_fill_manual(values = c("C57BL6J" = "red", "DBA2J" = "blue", "Other" = "grey")) + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )
#ggsave("score_data/Plots/mtor_score_by_strain_cd_age_corrected_data.png", width = 8, height = 6, dpi = 100)

# Filter strains that have more than 3 samples and that are on HF diet
mtor_score_df_with_metadata_filtered_hf <- mtor_score_df_with_metadata %>%
  group_by(Strain) %>%
  filter(n() > 3 & Diet == "HF") %>%
  ungroup()  # Ungroup to avoid issues with ggplot

# Convert Strain to character to avoid issues with colors
mtor_score_df_with_metadata_filtered_hf$Strain <- as.character(mtor_score_df_with_metadata_filtered_hf$Strain)

# Color all strains grey and the parents in different colors
mtor_score_df_with_metadata_filtered_hf <- mtor_score_df_with_metadata_filtered_hf %>%
  mutate(Color = ifelse(Strain %in% c("C57BL6J", "DBA2J"), Strain, "Other"))

# distribution of scores across different strains ordered based on median Score for HF diet
ggplot(mtor_score_df_with_metadata_filtered_hf, aes(x = reorder(Strain, -mtor_score_matrix, FUN = median), y = mtor_score_matrix, fill = Color)) + 
  geom_boxplot(alpha = 0.6) + 
  geom_jitter(width = 0.2, size = 1.5, color = "darkgreen", alpha = 0.5) +
  theme_minimal() + 
  ylim(-0.35, 0.35) +
  labs(title = "HF", x = "Strain", y = "mTOR Score") + 
  scale_fill_manual(values = c("C57BL6J" = "red", "DBA2J" = "blue", "Other" = "grey")) + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

#ggsave("score_data/Plots/mtor_score_by_strain_hf_age_corrected_data.png", width = 8, height = 6, dpi = 100)


### Score plots by sex ####

# distribution of scores by sex
ggplot(mtor_score_df_with_metadata, aes(x = Sex, y = mtor_score_matrix, fill = Sex)) + 
  geom_boxplot(alpha = 0.6) + 
  geom_jitter(width = 0.2, size = 1.5, color = "darkgreen", alpha = 0.5) +
  theme_minimal() + 
  labs(title = "Distribution of mTOR Scores by Sex", x = "Sex", y = "mTOR Score")


#### Difference between young and old ####

## Scores by diet for young and old ##

# Define the age threshold
age_threshold <- 450

# Add new column based on the age threshold
mtor_score_df_with_metadata_filtered <- mtor_score_df_with_metadata_filtered %>%
  mutate(Age_Group = ifelse(Age < age_threshold, "< 450 days", ">= 450 days"))

# Ensure Age_Group is treated as a factor
mtor_score_df_with_metadata_filtered$Age_Group <- factor(mtor_score_df_with_metadata_filtered$Age_Group, levels = c("< 450 days", ">= 450 days"))

# Plot distribution of mTOR scores by Diet and Age_Group
ggplot(mtor_score_df_with_metadata_filtered, aes(x = Age_Group, y = mtor_score_matrix, fill = Diet)) + 
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9), outliers = FALSE) +  # Boxplot for score distribution
  geom_jitter(size = 1.3, color = "black", alpha = 0.5, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9)) +  # Jittered points for individual scores
  theme_minimal() + 
  labs(title = "mTOR Scores by Age Group and Diet", x = "Age Group", y = "mTOR Score") +
  scale_fill_manual(values = c("blue", "darkgreen"), name = "Diet")  # Adjust colors for the fill



## Differences between young and old within the strains with more than 3 samples ##

# Calculate the mean score for each Age_Group, Strain, and Diet
mtor_diff_data <- mtor_score_df_with_metadata_filtered %>%
  group_by(Strain, Age_Group, Diet) %>%
  summarise(mean_mtor_score = mean(mtor_score_matrix, na.rm = TRUE), .groups = "drop")

# Calculate difference in score between Age Groups within each Strain and Diet
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
  labs(title = "Difference in mTOR Scores Between Age Groups by Strain",
       subtitle = "X-axis: Difference in Low-Fat Diet (CD), Y-axis: Difference in High-Fat Diet (HFD)",
       x = "Difference in mTOR Score (CD)",
       y = "Difference in mTOR Score (HF)")
#ggsave("score_data/Plots/mtor_score_difference_by_strain_uncorrected_data.png", width = 8, height = 6, dpi = 100)



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


#### Difference between the 4 age groups ####

# Define the age thresholds
age_threshold1 <- 300
age_threshold2 <- 450
age_threshold3 <- 650

# Add new column based on the age thresholds
mtor_score_df_with_metadata_filtered <- mtor_score_df_with_metadata_filtered %>%
  mutate(Age_Group = case_when(
    Age < age_threshold1 ~ "< 300 days",
    Age >= age_threshold1 & Age < age_threshold2 ~ "300-450 days",
    Age >= age_threshold2 & Age < age_threshold3 ~ "450-650 days",
    Age >= age_threshold3 ~ ">= 650 days"
  ))

# Ensure Age_Group is treated as a factor
mtor_score_df_with_metadata_filtered$Age_Group <- factor(mtor_score_df_with_metadata_filtered$Age_Group, levels = c("< 300 days", "300-450 days", "450-650 days", ">= 650 days"))

# Plot distribution of scores by Diet and Age_Group and add a line for the trends
ggplot(mtor_score_df_with_metadata_filtered, aes(x = Age_Group, y = mtor_score_matrix, fill = Diet)) + 
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9), outliers = FALSE) +  # Boxplot for score distribution
  geom_jitter(size = 1.3, color = "black", alpha = 0.5, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9)) +  # Jittered points for individual scores
  theme_minimal() + 
  labs(title = "mTOR Scores by 4 Age Groups and Diet", x = "Age Groups", y = "mTOR Score") +
  scale_fill_manual(values = c("blue", "darkgreen"), name = "Diet")  # Adjust colors for the fill
#ggsave("score_data/Plots/mtor_score_by_age_group_and_diet_uncorrected_data.png", width = 8, height = 6, dpi = 100)



#### Pairwise y-axis differences and record if same strain ####

# Generate all pairs of rows
pairwise_scores <- expand.grid(
  Row1 = 1:nrow(mtor_score_df_with_metadata_filtered),
  Row2 = 1:nrow(mtor_score_df_with_metadata_filtered)
)

# Remove pairs where Row1 and Row2 are the same
pairwise_scores <- pairwise_scores %>%
  filter(Row1 != Row2)

# Calculate pairwise distances and check if same strain
pairwise_scores <- pairwise_scores %>%
  mutate(
    Score1 = mtor_score_df_with_metadata_filtered$mtor_score_matrix[Row1],
    Score2 = mtor_score_df_with_metadata_filtered$mtor_score_matrix[Row2],
    Strain1 = mtor_score_df_with_metadata_filtered$Strain[Row1],
    Strain2 = mtor_score_df_with_metadata_filtered$Strain[Row2],
    # Calculate absolute difference in scores (y-axis distance)
    distance = abs(Score1 - Score2),
    # Check whether they are from same strain
    same_strain = Strain1 == Strain2
  )

# Boxplot of distances between same and different strains
ggplot(pairwise_scores, aes(x = same_strain, y = distance, fill = same_strain)) +
  geom_boxplot(alpha = 0.6) +
  theme_minimal() +
  labs(title = "Pairwise mTOR Score Differences (Same vs Different Strain)",
       x = "Same Strain",
       y = "Pairwise Distance (mTOR Scores)") +
  scale_x_discrete(labels = c("Different Strain", "Same Strain")) +
  scale_fill_manual(values = c("skyblue", "orange")) +
  theme(legend.position = "none")

# Distribution of distances for same vs different strains
ggplot(pairwise_scores, aes(x = distance, fill = same_strain)) +
  geom_density(alpha = 0.6) +
  theme_minimal() +
  labs(title = "Density of Pairwise mTOR Score Differences",
       x = "Pairwise Distance (mTOR Scores)",
       fill = "Same Strain") +
  scale_fill_manual(values = c("skyblue", "orange")) +
  theme(legend.position = "top")






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


### Score plots across strains for each diet separately - only necessary for age corrected data (for QTL) ####

# Filter strains that have more than 3 samples and that are on CD diet
erk_score_df_with_metadata_filtered_cd <- erk_score_df_with_metadata %>%
  group_by(Strain) %>%
  filter(n() > 3 & Diet == "CD") %>%
  ungroup()  # Ungroup to avoid issues with ggplot

# Convert Strain to character to avoid issues with colors
erk_score_df_with_metadata_filtered_cd$Strain <- as.character(erk_score_df_with_metadata_filtered_cd$Strain)

# Color all strains grey and the parents in different colors
erk_score_df_with_metadata_filtered_cd <- erk_score_df_with_metadata_filtered_cd %>%
  mutate(Color = ifelse(Strain %in% c("C57BL6J", "DBA2J"), Strain, "Other"))

# distribution of scores across different strains ordered based on median Score for CD diet
ggplot(erk_score_df_with_metadata_filtered_cd, aes(x = reorder(Strain, -erk_score_matrix, FUN = median), y = erk_score_matrix, fill = Color)) + 
  geom_boxplot(alpha = 0.6) + 
  geom_jitter(width = 0.2, size = 1.5, color = "darkgreen", alpha = 0.5) +
  theme_minimal() + 
  ylim(-0.35,0.35) +
  labs(title = "CD", x = "Strain", y = "ERK Score") + 
  scale_fill_manual(values = c("C57BL6J" = "red", "DBA2J" = "blue", "Other" = "grey")) + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )
#ggsave("score_data/Plots/erk_score_by_strain_cd_age_corrected_data.png", width = 8, height = 6, dpi = 100)

# Filter strains that have more than 3 samples and that are on HF diet
erk_score_df_with_metadata_filtered_hf <- erk_score_df_with_metadata %>%
  group_by(Strain) %>%
  filter(n() > 3 & Diet == "HF") %>%
  ungroup()  # Ungroup to avoid issues with ggplot

# Convert Strain to character to avoid issues with colors
erk_score_df_with_metadata_filtered_hf$Strain <- as.character(erk_score_df_with_metadata_filtered_hf$Strain)

# Color all strains grey and the parents in different colors
erk_score_df_with_metadata_filtered_hf <- erk_score_df_with_metadata_filtered_hf %>%
  mutate(Color = ifelse(Strain %in% c("C57BL6J", "DBA2J"), Strain, "Other"))

# distribution of scores across different strains ordered based on median Score for HF diet and remove legend
ggplot(erk_score_df_with_metadata_filtered_hf, aes(x = reorder(Strain, -erk_score_matrix, FUN = median), y = erk_score_matrix, fill = Color)) + 
  geom_boxplot(alpha = 0.6) + 
  geom_jitter(width = 0.2, size = 1.5, color = "darkgreen", alpha = 0.5) +
  theme_minimal() + 
  ylim(-0.35,0.35) +
  labs(title = "HF", x = "Strain", y = "ERK Score") + 
  scale_fill_manual(values = c("C57BL6J" = "red", "DBA2J" = "blue", "Other" = "grey")) + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )
#ggsave("score_data/Plots/erk_score_by_strain_hf_age_corrected_data.png", width = 8, height = 6, dpi = 100)


### Score plots by sex ####

# distribution of scores by sex
ggplot(erk_score_df_with_metadata, aes(x = Sex, y = erk_score_matrix, fill = Sex)) + 
  geom_boxplot(alpha = 0.6) + 
  geom_jitter(width = 0.2, size = 1.5, color = "darkgreen", alpha = 0.5) +
  theme_minimal() + 
  labs(title = "Distribution of ERK Scores by Sex", x = "Sex", y = "ERK Score")



#### Difference between young and old ####

## Scores by diet for young and old ##

# Define the age threshold
age_threshold <- 450

# Add new column based on the age threshold
erk_score_df_with_metadata_filtered <- erk_score_df_with_metadata_filtered %>%
  mutate(Age_Group = ifelse(Age < age_threshold, "< 450 days", ">= 450 days"))

# Ensure Age_Group is treated as a factor
erk_score_df_with_metadata_filtered$Age_Group <- factor(erk_score_df_with_metadata_filtered$Age_Group, levels = c("< 450 days", ">= 450 days"))

# Plot distribution of scores by Diet and Age_Group
ggplot(erk_score_df_with_metadata_filtered, aes(x = Age_Group, y = erk_score_matrix, fill = Diet)) + 
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9), outliers = FALSE) +  # Boxplot for score distribution
  geom_jitter(size = 1.3, color = "black", alpha = 0.5, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9)) +  # Jittered points for individual scores
  theme_minimal() + 
  labs(title = "ERK Scores by Age Group and Diet", x = "Age Group", y = "ERK Score") +
  scale_fill_manual(values = c("blue", "darkgreen"), name = "Diet")  # Adjust colors for the fill



## Differences between young and old within the strains with more than 3 samples ##

# Calculate the mean score for each Age_Group, Strain, and Diet
erk_diff_data <- erk_score_df_with_metadata_filtered %>%
  group_by(Strain, Age_Group, Diet) %>%
  summarise(mean_erk_score = mean(erk_score_matrix, na.rm = TRUE), .groups = "drop")

# Calculate difference in score between Age Groups within each Strain and Diet
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
  labs(title = "Difference in ERK Scores Between Younger and Older Mice by Strain",
       subtitle = "X-axis: Difference in Low-Fat Diet (CD), Y-axis: Difference in High-Fat Diet (HFD)",
       x = "Difference in ERK Score (CD)",
       y = "Difference in ERK Score (HF)")
#ggsave("score_data/Plots/erk_score_difference_by_strain_uncorrected_data.png", width = 8, height = 6, dpi = 100)



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


#### Difference between the 4 age groups #####

# Define the age thresholds
age_threshold1 <- 300
age_threshold2 <- 450
age_threshold3 <- 650

# Add new column based on the age thresholds
erk_score_df_with_metadata_filtered <- erk_score_df_with_metadata_filtered %>%
  mutate(Age_Group = case_when(
    Age < age_threshold1 ~ "< 300 days",
    Age >= age_threshold1 & Age < age_threshold2 ~ "300-450 days",
    Age >= age_threshold2 & Age < age_threshold3 ~ "450-650 days",
    Age >= age_threshold3 ~ ">= 650 days"
  ))

# Ensure Age_Group is treated as a factor
erk_score_df_with_metadata_filtered$Age_Group <- factor(erk_score_df_with_metadata_filtered$Age_Group, levels = c("< 300 days", "300-450 days", "450-650 days", ">= 650 days"))

# Plot distribution of scores by Diet and Age_Group and add a line for the trends
ggplot(erk_score_df_with_metadata_filtered, aes(x = Age_Group, y = erk_score_matrix, fill = Diet)) + 
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9), outliers = FALSE) +  # Boxplot for score distribution
  geom_jitter(size = 1.3, color = "black", alpha = 0.5, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9)) +  # Jittered points for individual scores
  theme_minimal() + 
  labs(title = "ERK Scores by 4 Age Groups and Diets", x = "Age Groups", y = "ERK Score") +
  scale_fill_manual(values = c("blue", "darkgreen"), name = "Diet")  # Adjust colors for the fill
#ggsave("score_data/Plots/erk_score_by_age_group_and_diet_uncorrected_data.png", width = 8, height = 6, dpi = 100)




#### Pairwise y-axis differences and record if same strain ####

# Generate all pairs of rows
pairwise_scores <- expand.grid(
  Row1 = 1:nrow(erk_score_df_with_metadata_filtered),
  Row2 = 1:nrow(erk_score_df_with_metadata_filtered)
)

# Remove pairs where Row1 and Row2 are the same
pairwise_scores <- pairwise_scores %>%
  filter(Row1 != Row2)

# Calculate pairwise distances and check if same strain
pairwise_scores <- pairwise_scores %>%
  mutate(
    Score1 = erk_score_df_with_metadata_filtered$erk_score_matrix[Row1],
    Score2 = erk_score_df_with_metadata_filtered$erk_score_matrix[Row2],
    Strain1 = erk_score_df_with_metadata_filtered$Strain[Row1],
    Strain2 = erk_score_df_with_metadata_filtered$Strain[Row2],
    # Calculate absolute difference in scores (y-axis distance)
    distance = abs(Score1 - Score2),
    # Check whether they are from same strain
    same_strain = Strain1 == Strain2
  )

# Boxplot of distances between same and different strains
ggplot(pairwise_scores, aes(x = same_strain, y = distance, fill = same_strain)) +
  geom_boxplot(alpha = 0.6) +
  theme_minimal() +
  labs(title = "Pairwise ERK Score Differences (Same vs Different Strain)",
       x = "Same Strain",
       y = "Pairwise Distance (ERK Scores)") +
  scale_x_discrete(labels = c("Different Strain", "Same Strain")) +
  scale_fill_manual(values = c("skyblue", "orange")) +
  theme(legend.position = "none")

# Distribution of distances for same vs different strains
ggplot(pairwise_scores, aes(x = distance, fill = same_strain)) +
  geom_density(alpha = 0.6) +
  theme_minimal() +
  labs(title = "Density of Pairwise ERK Score Differences",
       x = "Pairwise Distance (ERK Scores)",
       fill = "Same Strain") +
  scale_fill_manual(values = c("skyblue", "orange")) +
  theme(legend.position = "top")




#### Scores vs Lifespan ####

# Calculate median score for each strain with at least 2 samples
mtor_median_scores <- mtor_score_df_with_metadata %>%
  group_by(Strain) %>%
  filter(n() >= 3) %>%
  summarise(median_score = median(mtor_score_matrix, na.rm = TRUE), .groups = "drop")
erk_median_scores <- erk_score_df_with_metadata %>%
  group_by(Strain) %>%
  filter(n() >= 3) %>%
  summarise(median_score = median(erk_score_matrix, na.rm = TRUE), .groups = "drop")

# Calculate mean lifespan for each strain in survive
mean_lifespan <- survive %>%
  group_by(StrainNameCurrent) %>%
  summarise(mean_lifespan = mean(AgeAtDeath, na.rm = TRUE), .groups = "drop")

# Merge the median scores and median lifespan data
mtor_score_lifespan_data <- merge(mtor_median_scores, mean_lifespan, by.x = "Strain", by.y = "StrainNameCurrent")
erk_score_lifespan_data <- merge(erk_median_scores, mean_lifespan, by.x = "Strain", by.y = "StrainNameCurrent")

# Plot mean lifespan vs median score with strain labels
ggplot(mtor_score_lifespan_data, aes(x = mean_lifespan, y = median_score, label = Strain)) +
  geom_point(size = 3, color = "darkblue") +  # Points for each strain
  #geom_text(vjust = 1.7, hjust = 0.5, size.unit = "mm") +  # Add strain labels next to the points and make fint size small
  geom_smooth(method = "lm", se = FALSE) +  # Add linear regression line
  theme_minimal() +
  labs(title = "Mean Lifespan vs Median mTOR Score by Strain",
       x = "Mean Lifespan (Days)",
       y = "Median mTOR Score")
ggplot(erk_score_lifespan_data, aes(x = mean_lifespan, y = median_score, label = Strain)) +
  geom_point(size = 3, color = "darkblue") +  # Points for each strain
  #geom_text(vjust = 1.7, hjust = 0.5, size.unit = "mm") +  # Add strain labels next to the points and make fint size small
  geom_smooth(method = "lm", se = FALSE) +  # Add linear regression line
  theme_minimal() +
  labs(title = "Mean Lifespan vs Median ERK Score by Strain",
       x = "Mean Lifespan (Days)",
       y = "Median ERK Score")

# Calculate Correlation between Lifespan and Scores
correlation_mtor <- cor.test(mtor_score_lifespan_data$mean_lifespan, mtor_score_lifespan_data$median_score)
correlation_erk <- cor.test(erk_score_lifespan_data$mean_lifespan, erk_score_lifespan_data$median_score)





#### Scores vs Lifespan per Diet ####

# Calculate median score for each strain with at least 2 samples and diet
mtor_median_scores_diet <- mtor_score_df_with_metadata %>%
  group_by(Strain, Diet) %>%
  filter(n() >= 3) %>%
  summarise(median_score = median(mtor_score_matrix, na.rm = TRUE), .groups = "drop")
erk_median_scores_diet <- erk_score_df_with_metadata %>%
  group_by(Strain, Diet) %>%
  filter(n() >= 3) %>%
  summarise(median_score = median(erk_score_matrix, na.rm = TRUE), .groups = "drop")

# Merge the median scores and mean lifespan data
mtor_score_lifespan_data_diet <- merge(mtor_median_scores_diet, mean_lifespan, by.x = "Strain", by.y = "StrainNameCurrent")
erk_score_lifespan_data_diet <- merge(erk_median_scores_diet, mean_lifespan, by.x = "Strain", by.y = "StrainNameCurrent")

# Plot mean lifespan vs median score with strain labels
ggplot(mtor_score_lifespan_data_diet, aes(x = mean_lifespan, y = median_score, label = Strain)) +
  geom_point(size = 3, color = "darkblue") +  # Points for each strain
  #geom_text(vjust = 1.7, hjust = 0.5, size.unit = "mm") +  # Add strain labels next to the points and make fint size small
  geom_smooth(method = "lm", se = FALSE) +  # Add linear regression line
  theme_minimal() +
  labs(title = "Mean Lifespan vs Median mTOR Score by Strain and Diet",
       x = "Mean Lifespan (Days)",
       y = "Median mTOR Score") +
  facet_wrap(~ Diet)  # Separate plots by Diet
ggplot(erk_score_lifespan_data_diet, aes(x = mean_lifespan, y = median_score, label = Strain)) +
  geom_point(size = 3, color = "darkblue") +  # Points for each strain
  #geom_text(vjust = 1.7, hjust = 0.5, size.unit = "mm") +  # Add strain labels next to the points and make fint size small
  geom_smooth(method = "lm", se = FALSE) +  # Add linear regression line
  theme_minimal() +
  labs(title = "Mean Lifespan vs Median ERK Score by Strain and Diet",
       x = "Mean Lifespan (Days)",
       y = "Median ERK Score") +
  facet_wrap(~ Diet)  # Separate plots by Diet

# Calculate Correlation between Lifespan and Scores per Diet separately
correlation_mtor_CD_Diet <- cor.test(mtor_score_lifespan_data_diet$mean_lifespan[mtor_score_lifespan_data_diet$Diet == "CD"], mtor_score_lifespan_data_diet$median_score[mtor_score_lifespan_data_diet$Diet == "CD"])
correlation_mtor_HF_Diet <- cor.test(mtor_score_lifespan_data_diet$mean_lifespan[mtor_score_lifespan_data_diet$Diet == "HF"], mtor_score_lifespan_data_diet$median_score[mtor_score_lifespan_data_diet$Diet == "HF"])
correlation_erk_CD_Diet <- cor.test(erk_score_lifespan_data_diet$mean_lifespan[erk_score_lifespan_data_diet$Diet == "CD"], erk_score_lifespan_data_diet$median_score[erk_score_lifespan_data_diet$Diet == "CD"])
correlation_erk_HF_Diet <- cor.test(erk_score_lifespan_data_diet$mean_lifespan[erk_score_lifespan_data_diet$Diet == "HF"], erk_score_lifespan_data_diet$median_score[erk_score_lifespan_data_diet$Diet == "HF"])

# Print p-values of the correlation results for the different Diets
cat("Correlation between mean Lifespan and median mTOR Score for CD Diet:", correlation_mtor_CD_Diet$p.value, "\n")
cat("Correlation between mean Lifespan and median mTOR Score for HF Diet:", correlation_mtor_HF_Diet$p.value, "\n")
cat("Correlation between mean Lifespan and median ERK Score for CD Diet:", correlation_erk_CD_Diet$p.value, "\n")
cat("Correlation between mean Lifespan and median ERK Score for HF Diet:", correlation_erk_HF_Diet$p.value, "\n")




#### Heritability of the Scores ####

## Using mixed-effects model like in  https://doi.org/10.1093/sleep/zsz278 ##
### and Jan ###

library(lme4)

# Fit a mixed-effects model for scores with Strain as a random effect
mtor_mixed_model <- lmer(mtor_score_matrix ~ (1 | Strain), data = mtor_score_df_with_metadata)
erk_mixed_model <- lmer(erk_score_matrix ~ (1 | Strain), data = erk_score_df_with_metadata)

# View the variance components
summary(mtor_mixed_model)
summary(erk_mixed_model)

# Extract the variance components (genetic and residual variance)
mtor_var_components <- as.data.frame(VarCorr(mtor_mixed_model))
erk_var_components <- as.data.frame(VarCorr(erk_mixed_model))

# Calculate broad-sense heritability
mtor_genetic_variance <- mtor_var_components$vcov[1]  # Genetic variance
erk_genetic_variance <- erk_var_components$vcov[1]  # Genetic variance
mtor_residual_variance <- mtor_var_components$vcov[2]  # Residual variance
erk_residual_variance <- erk_var_components$vcov[2]  # Residual variance
mtor_H2 <- mtor_genetic_variance / (mtor_genetic_variance + mtor_residual_variance)
erk_H2 <- erk_genetic_variance / (erk_genetic_variance + erk_residual_variance)

print(paste("Broad-sense heritability mTOR Score (H²):", mtor_H2))
print(paste("Broad-sense heritability ERK Score (H²):", erk_H2))


#### Heritability of the Scores per Diet ####

# Load Score Dataframes of data corrected for age
mtor_pheno_bxd <- readRDS(paste0(DATA_DIR,"/BXD/score_data/bxd_mtor_score_df_with_metadata_corrected_age.rds"))
erk_pheno_bxd <- readRDS(paste0(DATA_DIR,"/BXD/score_data/bxd_erk_score_df_with_metadata_corrected_age.rds"))

# Separate the score data into two dataframes (one for each diet as we will perform QTL for each diet separately)
mtor_pheno_bxd_diet_CD <- mtor_pheno_bxd[mtor_pheno_bxd$Diet == "CD", ]
mtor_pheno_bxd_diet_HF <- mtor_pheno_bxd[mtor_pheno_bxd$Diet == "HF", ]
erk_pheno_bxd_diet_CD <- erk_pheno_bxd[erk_pheno_bxd$Diet == "CD", ]
erk_pheno_bxd_diet_HF <- erk_pheno_bxd[erk_pheno_bxd$Diet == "HF", ]

# Fit a mixed-effects model for scores with Strain as a random effect for each diet
mtor_mixed_model_CD <- lmer(mtor_score_matrix ~ (1 | Strain), data = mtor_pheno_bxd_diet_CD)
mtor_mixed_model_HF <- lmer(mtor_score_matrix ~ (1 | Strain), data = mtor_pheno_bxd_diet_HF)
erk_mixed_model_CD <- lmer(erk_score_matrix ~ (1 | Strain), data = erk_pheno_bxd_diet_CD)
erk_mixed_model_HF <- lmer(erk_score_matrix ~ (1 | Strain), data = erk_pheno_bxd_diet_HF)

# View the variance components
summary(mtor_mixed_model_CD)
summary(mtor_mixed_model_HF)
summary(erk_mixed_model_CD)
summary(erk_mixed_model_HF)

# Extract the variance components (genetic and residual variance) for each diet
mtor_var_components_CD <- as.data.frame(VarCorr(mtor_mixed_model_CD))
mtor_var_components_HF <- as.data.frame(VarCorr(mtor_mixed_model_HF))
erk_var_components_CD <- as.data.frame(VarCorr(erk_mixed_model_CD))
erk_var_components_HF <- as.data.frame(VarCorr(erk_mixed_model_HF))

# Calculate broad-sense heritability for each diet
mtor_genetic_variance_CD <- mtor_var_components_CD$vcov[1]  # Genetic variance
mtor_residual_variance_CD <- mtor_var_components_CD$vcov[2]  # Residual variance
mtor_H2_CD <- mtor_genetic_variance_CD / (mtor_genetic_variance_CD + mtor_residual_variance_CD)
mtor_genetic_variance_HF <- mtor_var_components_HF$vcov[1]  # Genetic variance
mtor_residual_variance_HF <- mtor_var_components_HF$vcov[2]  # Residual variance
mtor_H2_HF <- mtor_genetic_variance_HF / (mtor_genetic_variance_HF + mtor_residual_variance_HF)
erk_genetic_variance_CD <- erk_var_components_CD$vcov[1]  # Genetic variance
erk_residual_variance_CD <- erk_var_components_CD$vcov[2]  # Residual variance
erk_H2_CD <- erk_genetic_variance_CD / (erk_genetic_variance_CD + erk_residual_variance_CD)
erk_genetic_variance_HF <- erk_var_components_HF$vcov[1]  # Genetic variance
erk_residual_variance_HF <- erk_var_components_HF$vcov[2]  # Residual variance
erk_H2_HF <- erk_genetic_variance_HF / (erk_genetic_variance_HF + erk_residual_variance_HF)

print(paste("Broad-sense heritability mTOR Score for CD Diet (H²):", mtor_H2_CD))
print(paste("Broad-sense heritability mTOR Score for HF Diet (H²):", mtor_H2_HF))
print(paste("Broad-sense heritability ERK Score for CD Diet (H²):", erk_H2_CD))
print(paste("Broad-sense heritability ERK Score for HF Diet (H²):", erk_H2_HF))

# The difference in heritability between the two diets can be due to the different strains that are present
