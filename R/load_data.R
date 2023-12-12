
# load_data.R

#' Load and Preprocess Data
#'
#' Loads and preprocesses the data from a given CSV file.
#'
#' @param datacsv Path to the CSV file containing the data.
#' @param scaledf Boolean indicating whether to scale the data.
#' @param corrm Boolean indicating whether to apply a correlation matrix.
#' @param rmvna Boolean indicating whether to remove columns with NA or infinite values.
#' @param appendname Boolean indicating whether to append names to the data.
#' @importFrom caret nearZeroVar
#' @return A preprocessed data frame.
#' @export

load_data <- function(datacsv, scaledf=FALSE, corrm=TRUE, appendname=FALSE) {

  Data <- read.csv(datacsv, header=TRUE, sep=",", quote="\"", dec=".", fill=TRUE, comment.char="")

  # # Removing rows with most NAs
  # Calculate the percentage of missing (NA or empty) data for each patient (row)
  missing_percentage <- rowSums(apply(Data, 2, is_missing)) / ncol(Data)
  # Filter out patients (rows) with more than 50% missing data
  Data <- Data[missing_percentage < 0.9, ]


  # Store and remove patient IDs
  subs <- Data$sub
  Data$sub <- NULL

  name_to_append <- get_mdroiname(datacsv)
  name_split <- strsplit(name_to_append, "_")[[1]]
  name <- name_split[1]

  raddos <- c('adc', 'dds', 'ctt', 'ct1', 'tt2',
              'swi', 'flr', 'tt1', 'smt', 'sht',
              'letd', 'letdfirst', 'lemd','lemdfirst',
              "mkm2d","mkm5d","mkm6d","mkm7d","mkm8d","mkm9d","mkm10d")

  # Conditional data processing based on file name
  if (name %in% c('mkm2m', 'mkm5m','mkm6m','mkm7m','mkm8m',
                  'mkm9m','mkm10m','letm','lemm')) {
    numeric_columns <- sapply(Data, is.numeric)
    #Data <- subset(Data, select=c(mean:max, std, d2:d98, v20, v15))
    Data <- Data[numeric_columns]
  } else if (name %in% raddos) {
    Data <- subset(Data, select=-c(1:which(colnames(Data) == "original_shape_Elongation") - 1))
  }
  #else {
  #Data <- subset(Data, select=-c(1:which(colnames(Data) == "cliohe_age") - 1))
  #  }


  # Removing columns with NA or infinite values
  Data <- Data[, sapply(Data, function(x) !any(is.na(x) | is.infinite(x)))]

  # Generating modified column names
  if (appendname) {
    modreg <- paste(name_to_append, "_", sep='')
    colnames(Data) <- paste0(modreg, colnames(Data))
  }

  # # Remove non-variant features and highly correlated features
  badCols <- nearZeroVar(Data)
  if (length(badCols > 0))
    Data <- Data[,-badCols]
  if (corrm) {
    Data <- filter_correlated_features(Data, 0.8)
  }

  # Scaling data
  if (scaledf) {
    Data <- as.data.frame(scale(Data))
  }

  Data$sub <- subs
  return(Data)
}

