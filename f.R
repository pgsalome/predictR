library(ROSE)
library(smotefamily)
library(caret)

handle_imbalance <- function(csv_path, clinical_csv) {
  # Read the original data
  data <- read.csv(csv_path)

  # Read the clinical data
  clinical_data <- read.csv(clinical_csv)

  # Assuming there's a common identifier like 'patient_id' in both datasets
  common_id <- "sub"  # Replace with the actual identifier

  # Merge the data based on the common identifier
  merged_data <- merge(data, clinical_data, by = common_id, all.x = TRUE)

  # Check if the "test" column exists in the dataframe
  if ("test" %in% colnames(merged_data)) {
    merged_data <- merged_data[, !colnames(merged_data) %in% "test"]
  }

  # Identify time_, status_ columns and categorical variables
  time_cols <- grep("^time_", names(merged_data), value = TRUE)
  status_cols <- grep("^status_", names(merged_data), value = TRUE)

  outcomes <- list()

  # Apply SMOTE from smotefamily
  for (status_col in status_cols) {
    target <- as.factor(merged_data[, status_col])
    predictors <- merged_data[, !names(merged_data) %in% c(time_cols, status_cols, common_id)]

    # Correctly using SMOTE
    tryCatch({
      smote_data <- SMOTE(target, predictors, K = 5)$data
      # Save only the balanced dataset
      outcomes[[paste0("SMOTE_", status_col)]] <- smote_data
    }, error = function(e) {
      message("Error in SMOTE for ", status_col, ": ", e$message)
    })
  }

  return(outcomes)
}

# Usage
csv_path <- '/path/to/your/original

# Usage
csv_path <- '/path/to/your/original/data.csv'
clinical_csv <- '/path/to/your/clinical/data.csv'
result_datasets <- handle_imbalance(csv_path, clinical_csv)



# Usage
csv_path<-'/home/pgsalome/R/toolbx/organized_folders/ricci_outcomes.csv'

clinical_csv <- '/home/pgsalome/R/toolbx/organized_folders/NTCP_TCP/clidmeraddos_features/ricci-images/features_cli_ohe.csv'
result_datasets <- handle_imbalance(csv_path,clinical_csv )

