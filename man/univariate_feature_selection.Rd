% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FeatureSelection.R
\name{univariate_feature_selection}
\alias{univariate_feature_selection}
\title{Unified Univariate Feature Selection}
\usage{
univariate_feature_selection(DT, S, CS)
}
\arguments{
\item{DT}{Data frame containing the features and outcome.}

\item{S}{Survival time or binary indicator for survival.}

\item{CS}{Censoring status.}
}
\value{
A data frame of selected features and their counts.
}
\description{
Chooses the correct univariate feature selection method based on the outcome type.
}
