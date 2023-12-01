#' Feature Selection Signature Analysis
#'
#' This function performs feature selection signature analysis on a set of data files. It loads data, performs feature filtering based on frequency, and combines all features to create a comprehensive dataframe. The function then saves the processed dataframes into the specified output directories.
#'

#' @param outputDir Directory where the output dataframes will be saved.

#'
#' @return A list of dataframes processed based on the provided feature set and criteria.
#' @export

generate_models <- function(outputDir) {
  significant_features_files <- Sys.glob(paste(outputDir,"signficant_features_results/*/*.csv", sep = "/"))
  significant_features_files <- filter_empty_csvs(significant_features_files)
  # Function to process a list of CSV files

  result <- group_dataframes_for_modelling(significant_features_files, N=10)
  missing_data_files <- result$missing_untouched

  imputed <- impute_data(missing_data_files)

}
