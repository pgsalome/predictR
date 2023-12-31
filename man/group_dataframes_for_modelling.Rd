% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_functions.R
\name{group_dataframes_for_modelling}
\alias{group_dataframes_for_modelling}
\title{Group Data Frames for Modelling}
\usage{
group_dataframes_for_modelling(file_paths, N)
}
\arguments{
\item{file_paths}{A character vector of file paths to CSV files to be processed.}

\item{N}{A positive numeric value representing the threshold ratio (1/N) for the feature-to-observations ratio.}
}
\value{
A list containing four sublists with the following categories:
  - no_missing_less_than_N: Data frames with no missing values and less than 1/N features-to-observations ratio.
  - missing_removed_less_than_N: Data frames with missing values removed and less than 1/N features-to-observations ratio.
  - missing_untouched: Data frames with missing values untouched.
  - missing_removed_more_than_N: Data frames with missing values removed and more than 1/N features-to-observations ratio.
}
\description{
This function processes a list of CSV files and categorizes them into four lists based on specific criteria:
  1. No missing values and less than 1/N features-to-observations ratio.
  2. Missing values removed and less than 1/N features-to-observations ratio.
  3. Missing values untouched.
  4. Missing values removed and more than 1/N features-to-observations ratio.
}
\examples{
# Specify a list of CSV files and the threshold value N
csv_files <- c("file1.csv", "file2.csv", "file3.csv")
result <- group_dataframes_for_modelling(csv_files, 5)
print(result)

}
