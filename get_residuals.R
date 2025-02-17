#### Function to perform Data Correction ####

# Load the necessary package
# if (!require("RcppEigen")) install.packages("RcppEigen", dependencies=TRUE)
# library(RcppEigen)

# Function to fit a linear model and return residuals for a single gene
# Arguments:
#   - data: A numeric vector containing the data for a single gene.
#   - mm: The model matrix containing covariates.
# Returns:
#   - A numeric vector of residuals from the linear model.

get_residuals <- function(data, mm){
  # Fit a linear model using the fastLm function and extract residuals
  out <- fastLm(mm, data, method=3)$residuals
  # Free up memory
  gc()
  return(out) 
}



# Example of how to use the function

# source("/.../get_residuals.R")
# 
# mm <- model.matrix(~ A + B, data = metadata)
# 
# # Calculate number of cores
# no_cores <- detectCores() - 1
# 
# # Initiate cluster, export "mm" to worker nodes and load library on each node
# cl <- makeCluster(no_cores)
# parallel::clusterExport(cl, varlist = c("mm", "get_residuals"))
# invisible(clusterEvalQ(cl, {library(RcppEigen)}))
# 
# # Correct for Age and Diet
# rna_data_filtered_corrected_a_b <- t(parApply(cl = cl, rna_data_filtered, 1, FUN=function(r){
#   get_residuals(r, mm)
# }))
# 
# # Stop cluster
# stopCluster(cl)
# 
# # Set column names as in rna_data_filtered
# colnames(rna_data_filtered_corrected_age_diet) <- colnames(rna_data_filtered)
# 
# # Save corrected data