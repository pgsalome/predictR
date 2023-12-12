# fsSurv.R

#' Main Survival Analysis Function
#'
#' This function performs survival analysis based on the provided dataset and outcome type.
#' It integrates various feature selection methods like univariate analysis, Boruta, and Lasso.
#'
#' @param datacsv Path to the CSV file containing the data.
#' @param outcomecsv Path to the CSV file containing outcome data.
#' @param outcometype The type of outcome to analyze ('os', 'pfs', 'tox', or 'group').
#' @param outputDir Directory to save the output files.

#' @export
fsSurv <- function (datacsv, outcomecsv, outcometype, outputDir, nrep,run_parallel ) {

  Data <- load_data(datacsv, scaledf = TRUE, corrm = TRUE, appendname = FALSE)
  outcome <- read.csv(outcomecsv, header = TRUE, sep = ",", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
  outcome <- switch(outcometype,
                    "os" = data.frame(sub = outcome$sub, CS = outcome$status_os, S = outcome$time_os),
                    "pfs" = data.frame(sub = outcome$sub, CS = outcome$status_pfs, S = outcome$time_pfs),
                    "tox" = data.frame(sub = outcome$sub, CS = outcome$status_tox, S = outcome$time_tox),
                    "group" = data.frame(sub = outcome$sub, CS = outcome$group, S = 'binary'))
  DT <- merge(Data, outcome, by = 'sub')

  feature_selection_loop(DT, outputDir_md, nrep,run_parallel )
}

