% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fsSignature.R
\name{fsSignature}
\alias{fsSignature}
\title{Feature Selection Signature Analysis}
\usage{
fsSignature(features_folder, outputDir, freq = 500)
}
\arguments{
\item{features_folder}{String path to the folder containing feature data files.}

\item{outputDir}{Directory where the output dataframes will be saved.}

\item{freq}{Frequency threshold for feature filtering, default is 850.}
}
\value{
A list of dataframes processed based on the provided feature set and criteria.
}
\description{
This function performs feature selection signature analysis on a set of data files. It loads data, performs feature filtering based on frequency, and combines all features to create a comprehensive dataframe. The function then saves the processed dataframes into the specified output directories.
}
