#### set global options and load libraries ####

DATA_DIR = "/cellfile/cellnet/erkTorQtl/data/BXD/processed_data"
setwd("/cellfile/cellnet/erkTorQtl/data/BXD")

library(RcppEigen)
library(parallel)

#### Load data ####

# Load RNA-seq data
#rna_data_filtered <- readRDS("/cellfile/cellnet/erkTorQtl/data/processed/bxd_rna_data_filtered.rds") # Genes x Samples

# Load the RNA-seq data that only contains the genes that are not missing in more than 90 % of the samples
rna_data_filtered <- readRDS(paste0(DATA_DIR, "/bxd_rna_data_filtered_more_90_per_samples.rds")) # Genes x Samples

# Load metadata
metadata <- readRDS(paste0(DATA_DIR, "/bxd_metadata.rds")) # Samples x Variables
# Make OmicsEarTag as rownames and remove the OmicsEarTag column
rownames(metadata) <- metadata$OmicsEarTag
metadata <- metadata[,-1]
# Only keep the samples that are present in the RNA-seq data
metadata <- metadata[colnames(rna_data_filtered),]



#### Inspect metadata ####

# Check the metadata
head(metadata)

# Convert variables to appropriate types
metadata$Sex <- as.factor(metadata$Sex)            
metadata$Strain <- as.factor(metadata$Strain)     
metadata$Diet <- as.factor(metadata$Diet)          
metadata$Order <- as.numeric(metadata$Order)       
metadata$Age <- as.numeric(metadata$Age)           
metadata$SacDate <- as.Date(metadata$SacDate, format = "%d.%m.%y")

# Plot the distribution of Age
hist(metadata$Age, breaks = 30, main = "Distribution of age", xlab = "Age")

# Plot the distribution of Diet
table(metadata$Diet)
barplot(table(metadata$Diet), main = "Distribution of diet", xlab = "Diet")

# Plot the distribution of Sex
table(metadata$Sex)
barplot(table(metadata$Sex), main ="Distribution of Sex", xlab ="Sex")

# Plot the distributioin of Strain
table(metadata$Strain)
barplot(table(metadata$Strain), main = "Distribution of Strain", xlab = "Strain", las = 2)

# Plot the distribution of SacDate
table(metadata$SacDate)
barplot(table(metadata$SacDate), main = "Distribution of SacDate", xlab = "SacDate", las = 2)

# Plot Diet by Age
boxplot(metadata$Age ~ metadata$Diet, main = "Age by Diet", xlab = "Diet", ylab = "Age")
#barplot(table(metadata$Diet, metadata$Age), main = "Age by Diet", xlab = "Diet", ylab = "Age", las = 2)

# Plot Sex by Age
boxplot(metadata$Age ~ metadata$Sex, main = "Age by Sex", xlab = "Sex", ylab = "Age (Days)")

# Plot Strain by Age
boxplot(metadata$Age ~ metadata$Strain, main = "Age by Strain", xlab = "Strain", ylab = "Age", las = 2)

# Plot SacDate by Age
boxplot(metadata$Age ~ metadata$SacDate, main = "Age by SacDate", xlab = "SacDate", ylab = "Age", las = 2)

# Plot Diet by Sex
table(metadata$Diet, metadata$Sex)
barplot(table(metadata$Diet, metadata$Sex), main = "Diet by Sex", xlab ="Sex", beside = TRUE, legend.text = TRUE, las = 1)

# Plot Diet by Strain
table(metadata$Diet, metadata$Strain)
barplot(table(metadata$Diet, metadata$Strain), main = "Diet by Strain", xlab = "Diet", beside = TRUE, las = 2)

# Plot Diet by SacDate
table(metadata$Diet, metadata$SacDate)
barplot(table(metadata$Diet, metadata$SacDate), main = "Diet by SacDate", xlab = "Diet", beside = TRUE, las = 2)

# Plot Sex by Strain
table(metadata$Sex, metadata$Strain)
barplot(table(metadata$Sex, metadata$Strain), main = "Sex by Strain", xlab = "Sex", beside = TRUE, las = 2)

# Get number of male samples
nrow(metadata[metadata$Sex == "M",])


#### Remove the male samples ####

# Remove the male samples from rna_data_filtered and metadata
rna_data_filtered <- rna_data_filtered[,metadata$Sex == "F"]
metadata <- metadata[metadata$Sex == "F",]

# Save the data with only the female samples
saveRDS(rna_data_filtered, paste0("processed_data/bxd_rna_data_filtered_female.rds"))



#### Data Correction ####

# # Function to fit linear model and return residuals for a single gene
# get_residuals <- function(data, mm){
#   out <- fastLm(mm, data, method=3)$residuals
#   gc()
#   return(out) 
# }

source("/cellfile/cellnet/erkTorQtl/Scripts_Lara/get_residuals.R")

## Correct for Age and Diet ##

# Check for na values 
any(is.na(metadata$Sex))
any(is.na(metadata$Diet))
any(is.na(metadata$Age))
any(is.na(metadata$Strain))
any(is.na(metadata$SacDate))

# Factorize Diet and make Age as numeric
metadata$Diet <- as.factor(metadata$Diet)
metadata$Age <- as.numeric(metadata$Age)


# Create Design matrix for Age
mm <- model.matrix(~Age, data = metadata)
# Create Design matrix for Age and Diet
#mm <- model.matrix(~Age + Diet, data = metadata)
## Create Design matrix for Age and Diet and the interaction between them
#mm <- model.matrix(~Age + Diet + Age:Diet, data = metadata)

# Calculate number of cores
no_cores <- detectCores() - 1

# Initiate cluster, export "mm" to worker nodes and load library on each node
cl <- makeCluster(no_cores)
parallel::clusterExport(cl, varlist = c("mm", "get_residuals"))
invisible(clusterEvalQ(cl, {library(RcppEigen)}))


# Correct for Age
rna_data_filtered_corrected_age <- t(parApply(cl = cl, rna_data_filtered, 1, FUN=function(r){
  get_residuals(r, mm)
}))
# Correct for Age and Diet
#rna_data_filtered_corrected_age_diet <- t(parApply(cl = cl, rna_data_filtered, 1, FUN=function(r){
#  get_residuals(r, mm)
#}))
## Correct for Age and Diet and the interaction between them
#rna_data_filtered_corrected_age_diet_plus_interaction <- t(parApply(cl = cl, rna_data_filtered, 1, FUN=function(r){
#  get_residuals(r, mm)
#}))

# Stop cluster
stopCluster(cl)

# Set column names as in rna_data_filtered
colnames(rna_data_filtered_corrected_age) <- colnames(rna_data_filtered)
#colnames(rna_data_filtered_corrected_age_diet) <- colnames(rna_data_filtered)
#colnames(rna_data_filtered_corrected_age_diet_plus_interaction) <- colnames(rna_data_filtered)

# Save the corrected data
saveRDS(rna_data_filtered_corrected_age, paste0("corrected_data/bxd_rna_data_filtered_female_corrected_age.rds"))
#saveRDS(rna_data_filtered_corrected_age_diet, paste0("corrected_data/bxd_rna_data_filtered_female_corrected_age_diet.rds"))
#saveRDS(rna_data_filtered_corrected_age_diet_plus_interaction, paste0("corrected_data/bxd_rna_data_filtered_female_corrected_age_diet_plus_interaction.rds"))

