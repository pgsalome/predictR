#' Impute Missing Data
#'
#' This function takes a list of dataframes and performs multiple imputation methods
#' on each dataframe. It returns a list of dataframes where each dataframe is imputed
#' using a different method.
#'
#' @param outputDir directory where the significant feature results is stored from which a list of dataframes with missing values to be imputed.
#' @examples
#' # Assuming df1 and df2 are dataframes with missing values
#' my_dfs <- list(df1, df2)
#' imputed_dfs <- impute_data(my_dfs)
#' names(imputed_dfs) # Lists the names of imputed dataframes
#' @importFrom mice mice complete
#' @importFrom Amelia amelia
#' @importFrom VIM kNN
#' @importFrom Hmisc impute aregImpute
#' @importFrom caret findCorrelation
#' @importFrom mi mi
#' @export
fsImpute <- function(outputDir) {

  significant_features_files <- Sys.glob(paste(outputDir,"*/*/*.csv", sep = "/"))
  significant_features_files <- filter_empty_csvs(significant_features_files)
  # Function to process a list of CSV files

  results <- group_dataframes_for_modelling(significant_features_files, N=10)
  for (outcome in names(results)) {
    message(outcome)
    # Access both groups and combine them into one list
    missing_data_files_low_ratio <- results[[outcome]]$Original_with_Missing_Low_Ratio
    missing_data_files_high_ratio <- results[[outcome]]$Original_with_Missing_High_Ratio
    missing_data_files <- c(missing_data_files_low_ratio, missing_data_files_high_ratio)

    imputed_dfs <- list()

    for (i in 1:length(missing_data_files)) {
      outputDir_md <- paste(outputDir,"signficant_features_results",outcome, sep = "/")
      name <- names(missing_data_files[i])
      message(missing_data_files[i])
      message(name)
      df <- missing_data_files[[i]]


      ###### remove correalted features using a median imputation ####
      correlated_indices <- caret::findCorrelation(cor(impute(df, median), use = "ev", method = "pearson"),
                                            cutoff = 0.85, verbose = FALSE, names = FALSE, exact = TRUE)
      df_cor <- df[, -correlated_indices, drop = FALSE]

      # Mean/Median/Mode Imputation
      message("Mean/Median/Mode Imputation")
      imputed_dfs[[paste0(name,'_mean')]] <- impute(df, mean)
      imputed_dfs[[paste0(name,'_median')]] <- impute(df, median)
      imputed_dfs[[paste0(name,'_mode')]] <- impute(df, mode)

      # KNN Imputation
      message("KNN Imputation")

      original_col_names <- names(df)
      imputed_dfs[[paste0(name,'_knn')]] <- kNN(df)[original_col_names]
      original_col_names <- names(df_cor)
      imputed_dfs[[paste0(name,'_knncor')]] <- kNN(df_cor)[original_col_names]

      # Bayes
      #message("Byes")
      # bayes_model <- impute_with_feature_elimination(df,method="bayes")
      # for (j in 1:length(bayes_model)) {
      #   imputed_dfs[[paste0(name,'_bayes', j)]] <- bayes_model[[j]]
      # }
      # bayescor_model <- impute_with_feature_elimination(df_cor,method="bayes")
      # for (j in 1:length(bayescor_model)) {
      #   imputed_dfs[[paste0('bayescor', j, '_', i)]] <- bayes_model[[j]]
      # }
      #

      # tis failed we can check it out later
      #mlmcor_model <- impute_with_feature_elimination(df_cor,method="mlm",n_datasets=2)
      #imputed_dfs[[paste0('mllcor_', i)]] <- impute_with_feature_elimination(df_cor,method="mlm",n_datasets=2)

      # Multiple Imputation
      message("MICE")
      # imputed_dfs[[paste0(name,'_mice')]] <- impute_with_feature_elimination(df,method="mice")
      # imputed_dfs[[paste0(name,'_micecor')]] <- impute_with_feature_elimination(df_cor,method="mice")

    }
    save_dataframes(imputed_dfs, outputDir_md)
  }


}
#' Impute Missing Data with Feature Elimination
#'
#' Attempts to impute missing data in a dataframe using the Amelia package. If Amelia
#' fails due to too many variables, the function iteratively removes the least important
#' feature (based on an importance number appended to the feature name) and retries
#' the imputation. The process continues until successful imputation or all features are exhausted.
#'
#' @param df A dataframe with missing values to be imputed.
#' @param max_attempts The maximum number of attempts to impute the data by
#'        eliminating features. Defaults to the number of columns in `df`.
#' @return A dataframe with imputed values. If imputation is successful, it returns
#'         the imputed dataframe. If all features are eliminated without successful
#'         imputation, the function stops with an error.
#' @examples
#' # Assuming df is a dataframe with missing values and appended importance numbers
#' imputed_df <- impute_with_feature_elimination(df)
#' @importFrom Amelia amelia
#' @export
impute_with_feature_elimination <- function(df, method, n_datasets = 1) {
  attempt <- 1
  max_attempt <- ncol(df)

  while (attempt <= max_attempt) {
    result <- NULL

    # Choose the imputation method based on the method argument
    if (method == "mlm") {
      # Try to run Amelia
      result_s <- tryCatch({
        amelia(df, m = n_datasets)$imputations$imp1
      }, error = function(e) {
        NULL
      })
    } else if (method == "bayes") {
      # Try to run Bayesian imputation
      result_s <- tryCatch({
        mi_result <- mi(df)
        complete(mi_result)
      }, error = function(e) {
        NULL
      })
    } else if (method == "mice") {
      # Try to run Bayesian imputation
      result_s <- tryCatch({
        mice_result <- mice(df)
        complete(mice_result)
      }, error = function(e) {
        NULL
      })
    }

    # Check if imputation ran successfully
    if (!is.null(result_s)) {
      return(result_s)
    }

    # If imputation failed, remove the least important feature
    feature_importance <- sapply(names(df), function(name) {
      # Extracting the importance number from the feature name
      as.numeric(sub(".*_", "", name))
    })

    # Identify the feature with the lowest importance
    least_important_feature <- names(df)[which.min(feature_importance)]

    # Remove the least important feature
    message("imputation failed removing lowest ranked feature and trying again")
    df <- df[, !names(df) %in% least_important_feature, drop = FALSE]
    # Increment the attempt counter
    attempt <- attempt + 1
  }

  stop("All features have been eliminated, but imputation was not successful.")
}
