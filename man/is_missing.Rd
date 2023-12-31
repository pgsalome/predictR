% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_functions.R
\name{is_missing}
\alias{is_missing}
\title{Check for NA or Empty String}
\usage{
is_missing(x)
}
\arguments{
\item{x}{A vector where each element is to be checked for being NA or an empty string.}
}
\value{
A logical vector, where TRUE indicates that an element is either NA or an empty string, and FALSE otherwise.
}
\description{
This function checks each element of a vector to determine if it is an NA value or an empty string. It is useful for identifying missing values in a dataset where missing data might be represented in multiple ways.
}
