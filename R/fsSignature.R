#' Feature Selection Signature Analysis
#'
#' This function performs feature selection signature analysis on a set of data files. It loads data, performs feature filtering based on frequency, and combines all features to create a comprehensive dataframe. The function then saves the processed dataframes into the specified output directories.
#'
#' @param features_folder String path to the folder containing feature data files.
#' @param outputDir Directory where the output dataframes will be saved.
#' @param freq Frequency threshold for feature filtering, default is 850.
#'
#' @return A list of dataframes processed based on the provided feature set and criteria.
#' @export

fsSignature <- function(features_folder, outputDir, freq = 500) {

  datacsvlists <- Sys.glob(paste(features_folder, "/*", sep = ""))
  files <- Sys.glob(paste(outputDir,"feature_selection_results/*/*/*.rds", sep = "/"))

  # Extract feature names from files and filter datacsvlists
  feature_names_from_files <- unique(sapply(files, function(file_path) basename(dirname(file_path))))
  datacsvlists <- Filter(function(csv_file) basename(gsub("features_|\\.csv$", "", csv_file))
                         %in% feature_names_from_files, datacsvlists)

  get_mdroidf <- function(csv_list, result_folder, outcometype, freq) {
    print(csv_list)
    name <- get_mdroiname(csv_list)

    # Load data
    df <- load_data(csv_list, corrm = FALSE, scaledf = TRUE, appendname = TRUE)

    # Get a list of .rds files in the result folder for the specific outcometype

    files <- Sys.glob(paste(outputDir,"feature_selection_results/*", name, "*.rds", sep = "/"))

    # Check if no files are present and display a warning message
    if (length(files) == 0) {
      warning(paste("No feature selection was performed for the feature set:", name))
    }
    extract_outcome <- function(file_path) {
      path_parts <- strsplit(file_path, "/")[[1]]
      path_parts[length(path_parts) - 2]  # Assuming 'outcome' is always at this fixed position in the path
    }

    # Split files based on outcome
    files_by_outcome <- split(files, sapply(files, extract_outcome))

    for (outcome in names(files_by_outcome)) {
      files <- files_by_outcome[[outcome]]
      # Initialize an empty list to store feature data and dataframes
      feature_data_list <- list()
      dfs <- list()

      # Read feature data for each method and create dataframes
      for (file_path in files) {
        # Get names for saving
        method <- sub("\\.rds$", "", basename(file_path))
        features <- readRDS(file_path)

        # Filter features based on frequency and remove "eccentricity" from features
        valid_features <- features$f[features$count > freq & features$f != "eccentricity"]
        valid_counts <- features$count[features$count > freq & features$f != "eccentricity"]

        # Check if there are any features left after filtering
        if (length(valid_features) > 0) {
          # Create a temporary dataframe with the filtered features
          temp_df <- df[, c(paste(name,valid_features,sep="_"), 'sub')]

          # Rename the columns by appending the counts
          colnames(temp_df)[colnames(temp_df) != 'sub'] <- paste(name,valid_features, valid_counts, sep = "_")

          # Store the new dataframe in dfs
          dfs[[method]] <- temp_df
        } else {
          dfs[[method]] <- data.frame(sub = df$sub)  # Assuming 'sub' is a required column
        }
      }
      combined_df <- Reduce(function(x, y) merge(x, y, by = "sub", all = TRUE), dfs)
      # Create the 'all' dataframe
      dfs$all <- remove_duplicate_signfeat(combined_df)
      outputDir_md <- paste(outputDir, "/signficant_features_results/", outcome, sep = '')
      save_dataframes(dfs, outputDir_md, name)
    }
    }

  datalist <- lapply(datacsvlists, get_mdroidf, result_folder = outputDir, freq = freq)

  # ######################################
  #  do another round
  outcomes <- basename(Sys.glob(paste(outputDir, "/signficant_features_results/*", sep = '')))
  for (outcome in outcomes) {
    outputDir_md <- paste(outputDir, "/signficant_features_results/", outcome, sep = '')
    significant_features_files <- Sys.glob(paste(outputDir, "/signficant_features_results/", outcome, "/*", sep = ''))
    significant_features_files <- filter_empty_csvs(significant_features_files)

    # Splitting the basenames at underscores and extracting the first two parts
    split_names <- strsplit(basename(significant_features_files), "_")
    first_two_names <- sapply(split_names, function(x) paste(x[1:min(length(x), 2)], collapse = "_"))

    # Grouping the files based on the first two names
    grouped_files <- split(significant_features_files, first_two_names)
    # Apply the function to each group
    merged_dfs <- lapply(grouped_files, merge_grouped_files)
    save_dataframes(merged_dfs, outputDir_md)

    ###### now for the last merge #######
    # Construct a regular expression to match files that end with one of the modalities and .csv

    significant_features_files <- Sys.glob(paste(outputDir, "/signficant_features_results/", outcome, "/*", sep = ''))

    # mr_modalities <- c("tt1","ct1","swi","tt2","flr","adc")
    # modality_pattern_mr <- paste(mr_modalities, collapse = "|")
    # regex_pattern_mr <- paste0("(", modality_pattern_mr, ")\\.csv$")

    # Define your regex patterns
    regex_patterns <- list(
      mr = "(tt1|ct1|swi|tt2|flr|adc)\\.csv$",
      mrcli = "(tt1|ct1|swi|tt2|flr|adc|cli)\\.csv$",
      ddem = "(lem|dds)\\.csv$",
      ddemcli = "(lem|dds|cli)\\.csv$",
      rt = "(lem|dds|ct)\\.csv$",
      rtcli = "(lem|dds|ct|cli)\\.csv$",
      rad = "(tt1|ct1|swi|tt2|flr|adc|ctt)\\.csv$",
      radcli = "(tt1|ct1|swi|tt2|flr|adc|ctt|cli)\\.csv$",
      all = ".*/[^_/]+_[^_/]+\\.csv$"
    )
    # Iterate over regex patterns and apply the processing function
    for (pattern_name in names(regex_patterns)) {
      group_and_save_dfs(regex_patterns[[pattern_name]],pattern_name, outputDir_md, significant_features_files)
    }
  }

}
