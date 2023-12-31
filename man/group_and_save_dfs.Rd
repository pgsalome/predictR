% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_functions.R
\name{group_and_save_dfs}
\alias{group_and_save_dfs}
\title{Process and Save Files}
\usage{
group_and_save_dfs(
  regex_pattern,
  pattern_name,
  outputDir_md,
  significant_features_files
)
}
\arguments{
\item{regex_pattern}{A regular expression pattern to filter significant feature files.}

\item{outputDir_md}{The directory where merged dataframes will be saved.}

\item{significant_features_files}{A character vector containing file paths of significant features.}
}
\value{
None (dataframes are saved to the specified directory)
}
\description{
This function filters significant feature files based on a regex pattern,
groups them by methods, merges the grouped files, and saves the merged
dataframes to the specified output directory.
}
\examples{
process_and_save_files(regex_pattern_mr, "output_dir", significant_features_files)

}
