### R script to perform Spatial adjustment
# 1) Identify outliers and replace them with NA
# 2) Spatial analysis using SpATS package. 
#  Outcomes: - PDF file with plots for raw values, fitted values, residuals, BLUPs, and spatial pattern in one page for each trait.
#            - A csv file with the predicted values of the traits per genotype (BLUPs) after spatial correction.
#            - A csv file with the observed and fitted values for each plot (position: row x column) in the field

library(SpATS)
library(RColorBrewer)
library(plyr)
library(ggplot2)
library(dplyr)


SpatialCorrectionSpATS <- function(dfr, pheno='HD', gt='Genotype', rn='row', cn='col', pmain='') {
  if (pmain == '') {
    pmain <- pheno
  }
  colnames(dfr)[colnames(dfr) == rn] <- 'rn'
  colnames(dfr)[colnames(dfr) == cn] <- 'cn'
  dfr$cnf <- as.factor(dfr$cn)
  dfr$rnf <- as.factor(dfr$rn)
  spatD <- SpATS(response = pheno, 
                 spatial = ~ SAP(cn, rn, nseg = c(5,10)),
                 genotype = gt, 
                 random = ~ rnf + cnf,
                 data = dfr, 
                 control = list(tolerance = 1e-04), 
                 genotype.as.random = TRUE)
  plot(spatD, main = pmain)
  print(summary(spatD))
  fitted_values <- spatD$fitted
  # Combine fitted values with original data (gt, rn, cn, and phenotype)
  dfrF <- data.frame(
    Genotype = dfr[[gt]],  # genotype
    row = dfr$rn,           # row
    col = dfr$cn,           # column
    Observed = dfr[[pheno]],# observed phenotype
    Fitted = fitted_values  # fitted phenotype
  )
  colnames(dfrF) <- c("Genotype", "row", "col", paste0(pheno, "_Observed"), paste0(pheno, "_Fitted"))
  # Predicted values
  dfrP <- predict(spatD, which = gt)
  varhd <- variogram(spatD)
  plot(varhd)
  return(list(dfrF = dfrF, dfrP = dfrP))
}

# Identify and replace outliers using IQR method, and count them
replace_outliers_with_na <- function(x) {
  Q1 <- quantile(x, 0.25, na.rm = TRUE)
  Q3 <- quantile(x, 0.75, na.rm = TRUE)
  IQR_value <- Q3 - Q1
  # Define lower and upper bounds
  lower_bound <- Q1 - 1.5 * IQR_value
  upper_bound <- Q3 + 1.5 * IQR_value
  # Find outliers
  outliers <- x < lower_bound | x > upper_bound
  outlier_values <- x[outliers]
  num_outliers <- sum(outliers, na.rm = TRUE)
  # Print number of outliers replaced
  cat("Number of outliers replaced with NA:", num_outliers, "\n")
  # Print list of outliers
  if (num_outliers > 0) {
    cat("Outlier values:", outlier_values, "\n")
  } else {
    cat("No outliers found.\n")
  }
  # Replace outliers with NA
  x[outliers] <- NA
  return(x)
}

## define file names
fn <- "Phenotype_data_for_SpATS.csv"
fnpdf <- "SpATS_plots.pdf"
report <- "Outliers_SpATS_report.txt"
fnFitted <- "SpATS_fitted_values.csv" 
fnPredicted <- "SpATS_BLUPs.csv"  


dfr <- read.csv(fn)
dim(dfr)  
head(dfr)
colnames(dfr)

# add a column with the genotype and its position (gt+row+col)
dfr$genotype_pos <- paste(dfr$Genotype, dfr$Row, dfr$Column, sep = "_")

trFc <- 5 # number of the column with the first trait in dfr
trLc <- (length(colnames(dfr))-1)
traits <- colnames(dfr)[trFc:trLc]
print(traits)

pdf(fnpdf)
newdfr <- c()
combined_fitted_data <- c()  # To collect fitted data for all traits

################################################################################
### 1) Remove outliers #########################################################
################################################################################

# Open a connection to the logfile report
report <- file(report, open = "a")

# Start sinking both output and messages to the log file connection
sink(report, type = "output")
sink(report, type = "message")

for (tr in traits) {
  mytr <- sub('_.*', '', tr)
  print(mytr)
  
  pheno <- tr
  gt <- 'Genotype'
  rn <- 'Row'
  cn <- 'Column'
  
  # Subset the necessary columns from dfr
  mydfr <- dfr[, c(pheno, gt, rn, cn)]
  head(mydfr)
  
  cat("Identifying outliers for: ", pheno, "\n")
  # Identify and replace outliers using IQR method
  mydfr[, pheno] <- replace_outliers_with_na(mydfr[, pheno])
  
  # Convert phenotype column to numeric
  mydfr[, pheno] <- as.numeric(mydfr[, pheno])
  head(mydfr)
  
  # Handle exceptions for missing or invalid data
  exceptions <- dfr[, pheno][is.na(mydfr[, pheno])]
  exceptions <- exceptions[!is.na(exceptions)]
  exceptions <- exceptions[exceptions != '']
  
  print(exceptions)
  
  pmain <- paste(study, year, mytr)
  
  
  ##########################################################
  ### 2) Spatial correction BLUPs ##########################
  ##########################################################  
  cat("Start spatial correction for: ", pheno, "\n")
  
  # Perform Spatial Correction and BLUPs calculation
  res <- try(SpatialCorrectionSpATS_BLUPs(mydfr, pheno = tr, gt = gt, rn = rn, cn = cn, pmain = pmain))
  
  # If no error occurred in the SpatialCorrectionSpATS_BLUPs function
  if (!any(class(res) %in% "try-error")) {
    
    # Check if 'predicted' or 'errors' columns exist in res$dfrP
    pred_cols <- grep('predicted', colnames(res$dfrP))
    err_cols <- grep('errors', colnames(res$dfrP))
    
    if (length(pred_cols) > 0 && length(err_cols) > 0) {
      colnames(res$dfrP)[pred_cols] <- tr
      colnames(res$dfrP)[err_cols] <- paste(tr, 'StdErr', sep = '-')
    } else {
      cat("No 'predicted' or 'errors' columns found for trait:", tr, "\n")
      next  # Skip this iteration if columns are missing
    }
    
    # Combine fitted values data for each trait
    fitted_data_trait <- res$dfrF
    
    # Ensure consistent column names for fitted values
    observed_col <- paste0(tr, "_Observed")
    fitted_col <- paste0(tr, "_Fitted")
    
    # If it's the first iteration, initialize combined_fitted_data
    if (is.null(nrow(combined_fitted_data))) {
      combined_fitted_data <- fitted_data_trait
    } else {
      # Combine fitted data for all traits using cbind
      combined_fitted_data <- cbind(
        combined_fitted_data, 
        fitted_data_trait[, c(observed_col, fitted_col)]
      )
    }
    
    # Collect predicted values (BLUEs)
    pred_cols <- grep(tr, colnames(res$dfrP))
    if (length(pred_cols) > 0) {
      if (is.null(nrow(newdfr))) {
        newdfr <- res$dfrP[, c(1, pred_cols)]  # First time initialize newdfr
      } else {
        # Merge with previous data if there are valid predicted columns
        newdfr <- merge(newdfr, res$dfrP[, c(1, pred_cols)], by = 'Genotype', all = TRUE)
      }
    } else {
      cat("No predicted values found for trait:", tr, "\n")
    }
  } else {
    cat("Error occurred for trait:", tr, "\n")
  }
}

# Stop sinking
sink(type = "message")
sink(type = "output")

# Close the log file report
close(report)

dev.off()

# Add the column genotype_pos to the fitted data
head(combined_fitted_data)
combined_fitted_data$genotype_pos <- paste(combined_fitted_data$Genotype, combined_fitted_data$row, combined_fitted_data$col, sep = "_")

subdfr <- dfr[, c("genotype_pos", "Replicate")]
head(subdfr)

merged_fitted_data <- merge(subdfr, combined_fitted_data, by.x = "genotype_pos", by.y = "genotype_pos", all.x =T, all.y =T) #all column and rows from both tables will be includede.
head(merged_fitted_data)

# Remove the column genotype_pos
merged_fitted_data <- merged_fitted_data[, !names(merged_fitted_data) %in% "genotype_pos"]
head(merged_fitted_data)

# Move "Genotype" to be in the first column
merged_fitted_data <- merged_fitted_data[, c("Genotype", setdiff(names(merged_fitted_data), "Genotype"))]
head(merged_fitted_data)

# Save the fitted data for all traits to a CSV file
write.csv(merged_fitted_data, fnFitted, row.names = FALSE, na = '')

# Process and save the final predicted (BLUPs) data
finaldfrP <- data.frame(apply(newdfr[, -1], 2, signif, 3))
finaldfrP <- cbind(newdfr[, 1], finaldfrP)
colnames(finaldfrP)[1] <- 'Genotype'

fnout <- paste0(fnPredicted)
write.csv(finaldfrP, fnout, row.names = FALSE, na = '')

