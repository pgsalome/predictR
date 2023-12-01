#' @title Transform longitudinal data into time-to-event format
#'
#' @import dplyr
#' @import tidyr
#' @import stringr
#' @param file_name The name of the file containing the longitudinal data.
#' @param start_date The name of the variable containing the start date of the observation period.
#' @param fu_date The prefix used for the variables containing the follow-up dates.
#' @param var The prefix used for the variables containing the variable of interest.
#' @param id_loc The column index or name containing the unique identifier for each observation (default is 1).
#' @param units The time units used for calculating the time to event (default is "days").
#' @param date_format The date format used in the input data (default is "%m/%d/%Y").
#'
#' @return A data frame in time-to-event format with one row per timepoint.
#' @export

transform_fuData <- function(file_name, start_date, fu_date, var, id_loc = 1, units = "days", date_format = "%m/%d/%Y") {
  df <- read.csv(file_name)

  # Convert date columns to Date type and calculate the difference
  df <- df %>%
    mutate(across(starts_with(fu_date), ~ as.Date(.x, format = date_format)) - as.Date(.[[start_date]], format = date_format), .names = "diff_{.col}")

  # Convert to long format
  df_long <- df %>%
    pivot_longer(cols = starts_with(var), names_to = "timepoint", values_to = "value") %>%
    pivot_longer(cols = starts_with("diff_"), names_to = "diff_timepoint", values_to = "diff_value")

  # Extract and match timepoint numbers
  df_long <- df_long %>%
    mutate(timepoint_num = str_extract(timepoint, "\\d+"),
           diff_timepoint_num = str_extract(diff_timepoint, "\\d+")) %>%
    filter(timepoint_num == diff_timepoint_num)

  # Select and rename columns
  df_final <- df_long %>%
    select({{ id_loc }}, value, diff_value) %>%
    rename(time_to_event = diff_value, variable_value = value)

  return(df_final)
}
