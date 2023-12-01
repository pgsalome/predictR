# helper_functions.R

#' Save Feature Files
#'
#' Saves feature data to CSV files with names derived from the provided names.
#'
#' @param feature_data A named list of data frames or lists of data frames.
#' @param outputDir Directory to save the output files.
#' @export
save_feature_files <- function(feature_data, outputDir) {
  if (!dir.exists(outputDir)) {
    dir.create(outputDir, recursive = TRUE)
  }

  for (var_name in names(feature_data)) {
    # Extract the part of the name after 'features_'
    modified_name <- if (grepl("features_", var_name)) {
      sub("features_", "", var_name)
    } else {
      var_name
    }

    if (is.data.frame(feature_data[[var_name]])) {
      # Handle a single data frame
      file_path <- file.path(outputDir, paste0(modified_name, ".csv"))
      write.csv(feature_data[[var_name]], file_path, row.names = FALSE)
    } else if (is.list(feature_data[[var_name]])) {
      # Handle a list of data frames
      list_of_dfs <- feature_data[[var_name]]
      for (sub_df_name in names(list_of_dfs)) {
        # Create a unique file name using the parent name and child name
        unique_file_name <- paste0(sub_df_name, ".csv")
        file_path <- file.path(outputDir, unique_file_name)
        write.csv(list_of_dfs[[sub_df_name]], file_path, row.names = FALSE)
      }
    }
  }
}

#' Save Dataframes to CSV Files
#'
#' This function iterates over a list of dataframes, each associated with a method, and saves them as CSV files in a specified output directory. The filenames are generated using a template 'sigfeats_method_name.csv', where 'sigfeats' is an acronym for 'significant features', and 'method' is replaced with the actual method name.
#'
#' @param dfs A list of dataframes, where each dataframe is associated with a method.
#' @param output_dir A string representing the directory path where the CSV files will be saved. The function will create the directory if it does not exist.
#' @param name An optional string used as part of the filename for each saved dataframe. If provided, it is concatenated with the method name to generate the full filename. If not provided, the filename will consist only of the method name.
#'
#' @return No return value; this function is called for its side effect of saving files.
#' @export
save_dataframes <- function(dfs, output_dir, name=NULL) {
  # Ensure the output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Iterate over each method and save its dataframe
  for (method in names(dfs)) {
    file_name <- if (!is.null(name) && nzchar(name)) {
      paste(output_dir, paste0( method, "_", name, ".csv"), sep = "/")
    } else {
      paste(output_dir, paste0(method, ".csv"), sep = "/")
    }
    write.csv(dfs[[method]], file_name, row.names = FALSE)
  }
}




# filter_correlated_features.R
#' Filter Highly Correlated Features
#'
#' Removes features that are highly correlated based on the provided cutoff.
#'
#' @param DT Data frame to perform correlation filtering on.
#' @param cutoff Correlation cutoff value to determine feature removal.
#' @return Data frame with highly correlated features removed.
#' @importFrom caret findCorrelation
#' @export
filter_correlated_features <- function(df, cutoff) {

  corMatrix =  cor(df, y = NULL, use = "ev", method = c("pearson"))
  correlated_indices = findCorrelation(
    corMatrix,
    cutoff = cutoff, # correlation coefficient
    verbose = FALSE,
    names = FALSE,
    exact = TRUE
  )
  df <- df[, -correlated_indices]
  return(df)
}




# create_correlation_heatmap.R
#' Create Correlation Heatmap
#'
#' Generates a heatmap for visualizing correlations in the data.
#'
#' @param DT Data frame to visualize correlations.
#' @param plotName Name of the output plot file.
#' @param method Correlation method to use.
#' @export
create_correlation_heatmap <- function(DT, plotName, method = "pearson") {
  png(plotName, res = 250, width = 4000, height = 6000)
  heatmap(cor(DT, use = "ev", method = method), col = colorRampPalette(c("green", "white", "red"))(20), symm = TRUE)
  dev.off()
}


#' Extract name information from a CSV file path
#'
#' This function extracts name information from a CSV file path by splitting
#' the filename based on underscores and returning a concatenated string.
#'
#' @param datacsv A character string representing the CSV file path.
#' @return A character string containing the extracted name information.
#' @export
get_mdroiname <- function(datacsv) {
  # Extract the base name of the file (i.e., without path)
  filename <- basename(datacsv)
  # Remove the file extension from the filename
  filename_no_ext <- sub("\\.[^.]*$", "", filename)
  # Split the filename based on underscores
  split_path <- strsplit(filename_no_ext, "_")[[1]]
  # Concatenate the desired parts of the split filename
  return(paste(split_path[2], split_path[3], sep='_'))
}


#' Check for NA or Empty String
#'
#' This function checks each element of a vector to determine if it is an NA value or an empty string. It is useful for identifying missing values in a dataset where missing data might be represented in multiple ways.
#' @param x A vector where each element is to be checked for being NA or an empty string.
#' @return A logical vector, where TRUE indicates that an element is either NA or an empty string, and FALSE otherwise.
#' @export
#'
#'
is_missing <- function(x) {
  is.na(x) | x == ""
}
#' Check for Missing Values in a CSV File
#'
#' This function reads a CSV file from the specified file path, checks if it contains any missing values (NAs), and returns TRUE if any missing values are found.
#' @param file_path The file path to the CSV file to be checked.
#' @return TRUE if the CSV file contains missing values (NAs), FALSE otherwise.

#' @export
# Create a function to check for missing values in a CSV file
check_missing_values <- function(file_path) {
  df <- read.csv(file_path)
  has_missing_values <- any(is.na(df))
  return(has_missing_values)
}



#' Check Features to Observations Ratio
#'
#' This function checks if the number of features in a dataframe is less than 1/N times the number of observations, where N is a user-defined parameter.
#'
#' @param df A dataframe to be checked.
#' @param N A positive numeric value representing the threshold ratio (1/N).
#'
#' @return TRUE if the number of features is less than 1/N times the number of observations, FALSE otherwise.
#'

#' @export
check_features_to_observations_ratio <- function(df, N=10) {
  num_features <- ncol(df)
  num_observations <- nrow(df)

  # Calculate the threshold for the feature-to-observations ratio
  threshold <- num_observations / N

  # Check if the number of features is less than the threshold
  if (num_features < threshold) {
    return(TRUE)  # Number of features is less than 1/N times the number of observations
  } else {
    return(FALSE)  # Number of features is equal to or greater than 1/N times the number of observations
  }
}



#' Filter Out Empty CSV Files
#'
#' Iterates over a list of CSV file paths, checks if each file is empty or only contains a header,
#' and excludes empty files. Returns the list of files that are not empty.
#'
#' @param file_paths A character vector of file paths to check.
#' @return A list of file paths that are not empty.
#' @examples
#' list_of_files = c("path/to/your/first_file.csv", "path/to/your/second_file.csv")
#' non_empty_files = filter_empty_csvs(list_of_files)
#' @export
filter_empty_csvs <- function(file_paths) {
  multi_column_files <- c()  # To keep track of files with more than one column

  for (file_path in file_paths) {
    if (file.exists(file_path)) {
      con <- file(file_path, open = "r")
      first_line <- readLines(con, n = 1)
      close(con)

      num_columns <- length(unlist(strsplit(first_line, ",")))
      if (num_columns > 1) {
        multi_column_files <- c(multi_column_files, file_path)
      }
    }
  }

  return(multi_column_files)  # Return the list of files with more than one column
}

#' Filter and Rename DataFrame Columns Based on Maximum Count
#'
#' This function processes a dataframe by identifying columns with similar base names and different numeric suffixes. It retains only the column with the highest numeric suffix for each base name and renames it to include this maximum count.
#'
#' @param df A dataframe with column names that include a base name and a numeric suffix.
#' @return A dataframe with filtered and renamed columns.
#' @examples
#' df <- data.frame(
#'   "feature1_123" = 1:10,
#'   "feature1_456" = 1:10,
#'   "feature2_789" = 21:30
#' )
#' filtered_df <- filter_columns(df)
#' @export
remove_duplicate_signfeat <- function(df) {
  # Split the column names into base name and number suffix
  split_names <- strsplit(names(df), "_")
  base_names <- sapply(split_names, function(x) if(x[1] != "sub") paste(head(x, -1), collapse = "_") else x[1])
  suffixes <- sapply(split_names, function(x) if(x[1] != "sub") tail(x, 1) else "")

  # Initialize an empty list for new column names
  new_col_names <- list()

  for (base_name in unique(base_names)) {
    indices <- which(base_names == base_name)
    suffix_vals <- suffixes[indices]

    if (base_name == "sub") {
      # Directly use 'sub' without modification
      new_col_names[["sub"]] <- df[["sub"]]
    } else {
      # Check if suffixes are numeric and convert
      if (all(grepl("^[0-9]+$", suffix_vals))) {
        numeric_suffixes <- as.numeric(suffix_vals)
      } else {
        # Handle non-numeric suffixes by assigning them a low value
        numeric_suffixes <- ifelse(grepl("^[0-9]+$", suffix_vals), as.numeric(suffix_vals), -Inf)
      }

      max_index <- indices[which.max(numeric_suffixes)]
      max_suffix <- suffixes[max_index]
      new_col_name <- paste(base_name, max_suffix, sep = "_")
      new_col_names[[new_col_name]] <- df[[names(df)[max_index]]]
    }
  }

  new_df <- do.call(cbind, new_col_names)
  return(new_df)
}
#' extract_modality_from_filename
#'
#' This function is designed to extract the modality from a filename. It assumes that the filename
#' format includes the modality as the second last part, followed by a file extension (e.g., .csv).
#'
#' @param filename A string representing the full path of the file, including its name and extension.
#'
#' @return A string representing the modality extracted from the filename. For example, if the input
#' filename is "path/to/feature_modality.csv", the function will return "modality".
#'
#' @examples
#' filename <- "/path/to/your/file/feature_modality.csv"
#' modality <- group_methods(filename)
#' print(modality)  # Should print 'modality'
#'
#' @export
# Function to extract modality from filename
extract_modality_from_filename <- function(filename) {
  parts <- unlist(strsplit(basename(filename), "_"))
  modality_with_extension <- parts[length(parts)-1]
  return(sub("\\.csv$", "", modality_with_extension))
}

#' Process and Save Files
#'
#' This function filters significant feature files based on a regex pattern,
#' groups them by methods, merges the grouped files, and saves the merged
#' dataframes to the specified output directory.
#'
#' @param regex_pattern A regular expression pattern to filter significant feature files.
#' @param outputDir_md  The directory where merged dataframes will be saved.
#' @param significant_features_files A character vector containing file paths of significant features.
#'
#' @return None (dataframes are saved to the specified directory)
#'
#' @examples
#' process_and_save_files(regex_pattern_mr, "output_dir", significant_features_files)
#'
#' @export
group_and_save_dfs <- function(regex_pattern,pattern_name, outputDir_md, significant_features_files) {
  significant_features_files_filtered <- grep(regex_pattern, significant_features_files, value = TRUE)
  method_groups <- split(significant_features_files_filtered,
                         sapply(significant_features_files_filtered, extract_modality_from_filename))
  merged_dfs <- lapply(method_groups, merge_grouped_files)
  save_dataframes(merged_dfs, outputDir_md,pattern_name)
}

#' Merge Grouped Files
#'
#' This function takes a vector of file paths, reads the CSV files, and merges them by a common column ('sub').
#' @param file_group A character vector containing file paths to be merged.
#' @return A merged dataframe containing data from all the input files.

#' @export
merge_grouped_files <- function(file_group) {
  merged_dataframe <- NULL

  # Iterate over each file path in the vector of file_group
  for (file_path in file_group) {
    if (file.exists(file_path)) {
      df <- read.csv(file_path)

      if (is.null(merged_dataframe)) {
        merged_dataframe <- df
      } else {
        merged_dataframe <- merge(merged_dataframe, df, by = 'sub', all = TRUE)
      }
    } else {
      warning(paste("File not found and will be skipped:", file_path))
    }
  }

  return(merged_dataframe)
}

#' Update Feature Counts in a List of Data Frames
#'
#' This function iterates over a list of feature data frames and a corresponding list of new features.
#' It updates the feature count in each data frame based on the new features provided.
#'
#' @param pca_feature_counts A list of data frames, each representing feature counts for a specific PCA method.
#' @param pca_features A list of new features corresponding to each PCA method.
#' @return A list of data frames with updated feature counts.
#' @examples
#' pca_features <- perform_pca_feature_selection(DT[train_ind, ])
#' pca_feature_counts <- initialize_feature_counts(pca_features)
#' updated_counts <- update_features_df_list(pca_feature_counts, pca_features)
# Helper function to update the features data frame
update_features_df_list <- function(pca_feature_counts, pca_features) {
  # Check if lengths are the same
  if (length(pca_feature_counts) != length(pca_features)) {
    stop("Length of pca_feature_counts and pca_features must be the same")
  }

  for (i in seq_along(pca_features)) {
    feature_df <- pca_feature_counts[[i]]
    features_list <- pca_features[[i]]

    # Update feature count
    for (f in features_list) {
      feature_df <- update_features_df(feature_df,f)
    }

    pca_feature_counts[[i]] <- feature_df
  }

  return(pca_feature_counts)
}
#' Update Feature Counts in a Data Frame
#'
#' Updates the count of each feature in a data frame. If a feature is new, it is added to the data frame
#' with a count of 1. If a feature already exists, its count is incremented.
#'
#' @param feature_df A data frame with columns 'f' for features and 'count' for their counts.
#' @param new_features A vector of new features to be added or updated in the feature_df.
#' @return The updated data frame with counts of features.
#' @examples
#' feature_df <- data.frame(f = character(), count = integer())
#' new_features <- c("feature1", "feature2")
#' updated_df <- update_features_df(feature_df, new_features)
update_features_df  <- function(feature_df, new_features) {

  for (f in new_features) {
    if (f %in% feature_df$f) {
      # Increase count for existing features
      feature_df$count[feature_df$f == f] <- feature_df$count[feature_df$f == f] + 1
    } else {
      # Add new feature with a count of 1
      row <- data.frame(f=f, c=1)
      names(row) <- c("f", "count")
      feature_df <- rbind(feature_df, row)
    }
  }

  return(feature_df)
}

#' Group Data Frames for Modelling
#'
#' This function processes a list of CSV files and categorizes them into four lists based on specific criteria:
#'   1. No missing values and less than 1/N features-to-observations ratio.
#'   2. Missing values removed and less than 1/N features-to-observations ratio.
#'   3. Missing values untouched.
#'   4. Missing values removed and more than 1/N features-to-observations ratio.
#'
#' @param file_paths A character vector of file paths to CSV files to be processed.
#' @param N A positive numeric value representing the threshold ratio (1/N) for the feature-to-observations ratio.
#'
#' @return A list containing four sublists with the following categories:
#'   - no_missing_less_than_N: Data frames with no missing values and less than 1/N features-to-observations ratio.
#'   - missing_removed_less_than_N: Data frames with missing values removed and less than 1/N features-to-observations ratio.
#'   - missing_untouched: Data frames with missing values untouched.
#'   - missing_removed_more_than_N: Data frames with missing values removed and more than 1/N features-to-observations ratio.
#'
#' @examples
#' # Specify a list of CSV files and the threshold value N
#' csv_files <- c("file1.csv", "file2.csv", "file3.csv")
#' result <- group_dataframes_for_modelling(csv_files, 5)
#' print(result)
#'
#' @export
group_dataframes_for_modelling <- function(file_paths, N) {
  outcomes <- list()

  for (file_path in file_paths) {
    # Check if the file exists
    if (!file.exists(file_path)) {
      warning(paste("File not found:", file_path))
      next
    }

    # Read the CSV file
    df <- read.csv(file_path)
    file_name <- gsub("\\.csv$", "", basename(file_path))

    # Extract the outcome from the directory name
    outcome <- basename(dirname(file_path))
    if (!is.list(outcomes[[outcome]])) {
      outcomes[[outcome]] <- list(
        "Original_with_Missing_Low_Ratio" = list(),
        "Original_with_Missing_High_Ratio" = list(),
        "Cleaned_Low_Ratio" = list(),
        "Cleaned_High_Ratio" = list(),
        "No_Missing_Low_Ratio" = list(),
        "No_Missing_High_Ratio" = list()
      )
    }

    # Check for missing values and calculate the threshold
    has_missing_values <- anyNA(df)
    num_features <- ncol(df)
    num_observations <- nrow(df)
    threshold <- num_observations / N

    # Classify the dataframe into one of the six groups
    if (has_missing_values) {
      df_cleaned <- na.omit(df)

      if (num_features < threshold) {
        outcomes[[outcome]]$Original_with_Missing_Low_Ratio[[file_name]] <- df
        outcomes[[outcome]]$Cleaned_Low_Ratio[[file_name]] <- df_cleaned
      } else {
        outcomes[[outcome]]$Original_with_Missing_High_Ratio[[file_name]] <- df
        outcomes[[outcome]]$Cleaned_High_Ratio[[file_name]] <- df_cleaned
      }
    } else {
      if (num_features < threshold) {
        outcomes[[outcome]]$No_Missing_Low_Ratio[[file_name]] <- df
      } else {
        outcomes[[outcome]]$No_Missing_High_Ratio[[file_name]] <- df
      }
    }
  }

  return(outcomes)
}

#' Retrieve a Subset of CSV Lists Based on Included Names
#'
#' This function filters a list of CSV file paths, returning only those
#' whose associated names are included in a specified list.
#'
#' @param datacsvlists A list of CSV file paths.
#' @param to_include A vector of names to be included in the returned list.
#' @return A filtered list of CSV file paths.
#' @export

get_list <- function(datacsvlists, to_include) {
  Filter(function(csv_list) {
    name <- get_mdroiname(csv_list)
    name %in% to_include
  }, datacsvlists)
}

