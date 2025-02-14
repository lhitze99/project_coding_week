#### set global options and load libraries ####

DATA_DIR = "/cellfile/cellnet/erkTorQtl/data"
setwd("/cellfile/cellnet/erkTorQtl/data/BXD/processed_data")

library(openxlsx)
library(dplyr)
library(tidyr)


#### Load data ####

# Data reprocessed by the authors:
# Read counts were normalized to RPKM values using gene lengths from ENSEMBL82 v2015-10-02. 
# All RNA-seq data were then scaled by adding 1 to the normalized counts and then taking the log2.

# Load metadata about samples (eg.,  age, strain, diet, and sex) (8 x 627)
bxd_header <- read.table(paste0(DATA_DIR,"/BXD/raw_data/aData_S1_AllOmicsandPhenotypeData.csv"), nrows = 8, row.names = 1, header = T, sep = ",")
bxd_header <- bxd_header[,3:ncol(bxd_header)]

# Load data (rows = different features (e.g., mRNA counts), columns = samples) (76566 x 627)
# processed and normalized set of all omics data
bxd_data <- read.table(paste0(DATA_DIR,"/BXD/raw_data/aData_S1_AllOmicsandPhenotypeData.csv"), row.names = 1, header = T, na.strings = "", sep = ",", skip = 8, stringsAsFactors = FALSE)
bxd_data <- bxd_data[1:nrow(bxd_data),3:ncol(bxd_data)]

# Load marker lists 
mtor_positive_translation <- readRDS(paste0(DATA_DIR,"/marker_lists_mouse/mtor_markers_positive_mouse_without_Jade2.rds"))
mtor_negative_translation <- readRDS(paste0(DATA_DIR,"/marker_lists_mouse/mtor_markers_negative_mouse.rds"))
erk_positive_translation <- readRDS(paste0(DATA_DIR,"/marker_lists_mouse/erk_markers_positive_mouse_without_Txnrd1.rds"))
erk_negative_translation <- readRDS(paste0(DATA_DIR,"/marker_lists_mouse/erk_markers_negative_mouse_without_EIF4EBP1.rds"))



#### Prepare Data ####

## Extract mRNA data ##

# Identify which rows correspond to transcriptomic data (rownames starting with "mRNA_")
rna_data <- bxd_data[grep("^mRNA_", rownames(bxd_data)), ]

# Remove first rows because those are no genes
rna_data <- rna_data[-(1:70),]

# There are genes appearing multiple times but with another number

# Extract gene names from row names (the gene name follows the last underscore)
gene_names <- sub(".*_", "", rownames(rna_data))

# Add gene names as a new column in the data
rna_data_with_genes <- cbind(GeneName = gene_names, rna_data)

# Find the row with the highest expression (highest row sum) for each gene
rna_data_with_genes <- rna_data_with_genes[order(rna_data_with_genes$GeneName, -rowSums(rna_data_with_genes[, -1], na.rm = TRUE)),]

# Remove duplicates by keeping only the first occurrence of each gene name
rna_data_filtered <- rna_data_with_genes[!duplicated(rna_data_with_genes$GeneName),]

# Set row names to gene names (in the first column) 
rownames(rna_data_filtered) <- rna_data_filtered$GeneName

# Remove the GeneName column
rna_data_filtered <- rna_data_filtered[ , -1] 

# Function to check presence of marker genes
check_markers <- function(marker_list, gene_data) {
  present <- marker_list %in% gene_data
  if (all(present)) {
    print("All marker genes are present in the dataset.")
  } else {
    missing_genes <- marker_list[!present]
    print("The following marker genes are missing:")
    print(missing_genes)
  }
}

# All marker genes are present
check_markers(mtor_positive_translation, rownames(rna_data_filtered))
check_markers(mtor_negative_translation, rownames(rna_data_filtered))
check_markers(erk_positive_translation, rownames(rna_data_filtered))
check_markers(erk_negative_translation, rownames(rna_data_filtered))

# Remove samples (columns) that have only NA values
rna_data_filtered <- rna_data_filtered[, colSums(!is.na(rna_data_filtered)) > 0]

# Check if there are genes with only NA values
rownames(rna_data_filtered)[rowSums(is.na(rna_data_filtered)) == ncol(rna_data_filtered)]

# Save the filtered RNA data
# saveRDS(rna_data_filtered, paste0("bxd_rna_data_filtered.rds"))


## Extract relevant metadata ##

# Extract relevant metadata rows and transpose the header data
metadata <- t(bxd_header[c("OmicsEarTag", "Age", "Sex", "Strain", "Diet", "Order", "SacDate"), ])

# Convert to a data frame for easier manipulation
metadata <- as.data.frame(metadata)

# Convert data types
metadata$Age <- as.numeric(as.character(metadata$Age))
metadata$Sex <- as.factor(metadata$Sex)
metadata$Strain <- as.factor(metadata$Strain)
metadata$Diet <- as.factor(metadata$Diet)
metadata$SacDate <- as.factor(metadata$SacDate)
#metadata$OmicsEarTag <- as.factor(metadata$OmicsEarTag)

# Save the metadata
# saveRDS(metadata, paste0("bxd_metadata.rds"))




#### inspect the data ####

## missing samples per gene ##

# Number of non-zero samples for each gene
num_non_zero_samples <- rowSums(rna_data_filtered != 0)

# Number of zero samples for each gene
num_zero_samples <- rowSums(rna_data_filtered == 0)

# Percentage of non-zero samples for each gene
percentage_non_zero_samples <- (num_non_zero_samples / ncol(rna_data_filtered))

# Percentage of missing samples for each gene
percentage_zero_samples <- (num_zero_samples / ncol(rna_data_filtered))

# Commonly missing genes (genes that are missing in more than 10 % of the samples)
print(names(percentage_zero_samples[percentage_zero_samples > 0.1]))
print(length(percentage_zero_samples[percentage_zero_samples > 0.1]))

# Save genes that are not missing in more than 90 % of the samples 
gene_filter <- names(percentage_non_zero_samples[percentage_non_zero_samples>0.9])

# Make dataframe for easier plotting
df_non_zero_samples <- data.frame(Gene = rownames(rna_data_filtered), NumNonZeroSamples = num_non_zero_samples, 
                                  PercentageNonZeroSamples = percentage_non_zero_samples)

# Histogram of number of non-zero samples per gene
hist(df_non_zero_samples$NumNonZeroSamples, breaks = 200, main = "Distribution of Non-Zero Samples per Gene", 
     xlab = "Number of Non-Zero Samples", ylab = "Frequency")

# Histogram of number of zero samples per gene
hist(num_zero_samples, breaks = 200, main = "Distribution of Zero Samples per Gene", 
     xlab = "Number of Zero Samples", ylab = "Frequency")

# Keep only the genes that are not missing in more than 90 % of the samples
rna_data_filtered_sub <- rna_data_filtered[rownames(rna_data_filtered) %in% gene_filter,]

dim(rna_data_filtered)
dim(rna_data_filtered_sub)

# Check if all marker genes are still present
check_markers(mtor_positive_translation, rownames(rna_data_filtered_sub)) # 3 genes missing
check_markers(mtor_negative_translation, rownames(rna_data_filtered_sub))
check_markers(erk_positive_translation, rownames(rna_data_filtered_sub)) # 5 genes missing
check_markers(erk_negative_translation, rownames(rna_data_filtered_sub)) # 3 genes missing

# Save the names of the missing marker genes
missing_mtor_positive <- mtor_positive_translation[!mtor_positive_translation %in% rownames(rna_data_filtered_sub)]
# missing_mtor_negative <- mtor_negative_translation[!mtor_negative_translation %in% rownames(rna_data_filtered_sub)] # no gene missing
missing_erk_positive <- erk_positive_translation[!erk_positive_translation %in% rownames(rna_data_filtered_sub)]
missing_erk_negative <- erk_negative_translation[!erk_negative_translation %in% rownames(rna_data_filtered_sub)]

# Add the missing marker genes to rna_data_filtered_sub

# Add missing mTOR positive genes with their entries of rna_data_filtered
missing_mtor_positive_data <- rna_data_filtered[missing_mtor_positive,]
rna_data_filtered_sub <- rbind(rna_data_filtered_sub, missing_mtor_positive_data)

# Add missing ERK positive genes with their entries of rna_data_filtered
missing_erk_positive_data <- rna_data_filtered[missing_erk_positive,]
rna_data_filtered_sub <- rbind(rna_data_filtered_sub, missing_erk_positive_data)

# Add missing ERK negative genes with their entries of rna_data_filtered
missing_erk_negative_data <- rna_data_filtered[missing_erk_negative,]
rna_data_filtered_sub <- rbind(rna_data_filtered_sub, missing_erk_negative_data)

dim(rna_data_filtered_sub)

# Save the data that only contains the genes that are not missing in more than 90 % of the samples
# saveRDS(rna_data_filtered_sub, paste0("bxd_rna_data_filtered_more_90_per_samples.rds"))




## for the marker genes only ##
# change name for mtor_negative_translation, erk_positive_translation and erk_negative_translation
bxd_data_sub <- rna_data_filtered[rownames(rna_data_filtered) %in% unlist(mtor_positive_translation),]
num_non_zero_samples_mtor_pos <- rowSums(bxd_data_sub != 0)
num_zero_samples_mtor_pos <- rowSums(bxd_data_sub == 0)
percentage_non_zero_samples_mtor_pos <- (num_non_zero_samples_mtor_pos / ncol(bxd_data_sub))
df_non_zero_samples_mtor_pos <- data.frame(Gene = rownames(bxd_data_sub), NumNonZeroSamplesmTORpos = num_non_zero_samples_mtor_pos, 
                                           PercentageNonZeroSamplesmTORpos = percentage_non_zero_samples_mtor_pos)
barplot(df_non_zero_samples_mtor_pos$NumNonZeroSamplesmTORpos, names.arg = df_non_zero_samples_mtor_pos$Gene,
        las = 2, main = "Number of Non-Zero Samples per mTOR Positive Marker Genes", xlab = "Genes", 
        ylab = "Number of Non-Zero Samples", cex.names = 0.7)
hist(df_non_zero_samples_mtor_pos$NumNonZeroSamplesmTORpos, main = "Distribution of Non-Zero Samples per mTOR Positive Marker Gene", 
     xlab = "Number of Non-Zero Samples", ylab = "Frequency")





## missing genes per sample ##

# Number of non-zero genes for each sample
num_non_zero_genes <- colSums(rna_data_filtered != 0)

# Number of zero genes for each sample
num_zero_genes <- colSums(rna_data_filtered == 0)

# Percentage of non-zero genes for each sample
percentage_non_zero_genes <- (num_non_zero_genes / nrow(rna_data_filtered))

# Percentage of missing genes for each sample
percentage_zero_genes <- (num_zero_genes / nrow(rna_data_filtered))

# Commonly missing samples (samples that are missing in more than 10 % of the genes)
print(names(percentage_zero_genes[percentage_zero_genes > 0.1]))
print(length(percentage_zero_genes[percentage_zero_genes > 0.1]))

# Make dataframe for easier plotting
df_non_zero_genes <- data.frame(Samples = colnames(rna_data_filtered), NumNonZeroGenes = num_non_zero_genes, 
                                PercentageNonZeroGenes = percentage_non_zero_genes)

# Histogram of number of non-zero samples per gene
hist(df_non_zero_genes$NumNonZeroGenes, breaks = 200, main = "Distribution of Quantified (Non-Zero) Genes per Sample", 
     xlab = "Number of Non-Zero Genes", ylab = "Frequency")

# Histogram of number of zero genes per sample
hist(num_zero_genes, breaks = 200, main = "Distribution of Zero Genes per Sample", 
     xlab = "Number of Zero Genes", ylab = "Frequency")




## Explore potential Batch Effects ##

library(ggplot2)

# Merge metadata with dataframe of non-zero genes
df_non_zero_genes_with_metadata <- merge(df_non_zero_genes, metadata, by.x = "Samples", by.y = "OmicsEarTag")

# Plot number of non-zero genes per sample, colored by diet
ggplot(df_non_zero_genes_with_metadata, aes(x = NumNonZeroGenes, fill = Diet)) +
  geom_histogram(binwidth = 100, position = "identity", alpha = 0.5) +
  labs(title = "Distribution of Quantified (Non-Zero) Genes by Diet", x = "Number of Quantified Genes", y = "Frequency") +
  theme_minimal()

# Scatter plot of number of age as a function of non-zero genes
ggplot(df_non_zero_genes_with_metadata, aes(x = NumNonZeroGenes, y = Age)) +
  geom_point(alpha = 0.7) +  # Scatter plot of non-zero genes by age
  scale_y_continuous(breaks = seq(min(df_non_zero_genes_with_metadata$Age), 
                                  max(df_non_zero_genes_with_metadata$Age), 
                                  by = 20)) +  # Show every 10th day
  labs(title = "Number of Quantified (Non-Zero) Genes by Age", x = "Number of Quantified Genes", y = "Age (Days)") +
  theme_minimal() +
  theme(axis.text.y = element_text(vjust = 0.5, hjust=1))

# Boxplot of number of non-zero genes by strain
ggplot(df_non_zero_genes_with_metadata, aes(x = Strain, y = NumNonZeroGenes)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5, color = "blue") +  # Add jitter for individual points
  labs(title = "Distribution of Quantified Genes by Strain", x = "Strain", y = "Number of Quantified Genes") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels for readability

# Convert SacDate to a proper date format if it's not already
df_non_zero_genes_with_metadata$SacDate <- as.Date(df_non_zero_genes_with_metadata$SacDate, format = "%d.%m.%Y")
summary(df_non_zero_genes_with_metadata$SacDate)

# Scatter plot of non-zero genes by sacrifice date (ordered chronologically)
ggplot(df_non_zero_genes_with_metadata, aes(x = NumNonZeroGenes, y = SacDate)) +
  geom_point(alpha = 0.7) +  # Scatter plot of non-zero genes by sacrifice date
  labs(title = "Number of Quantified (Non-Zero) Genes by Sacrifice Date", x = "Number of Quantified Genes", y = "Sacrifice Date") +
  scale_y_date(date_labels = "%b-%Y", date_breaks = "1 month") +  # Full date format (Day-Month-Year)
  theme_minimal() +
  theme(axis.text.y = element_text(hjust = 1))  # Rotate dates for readability
