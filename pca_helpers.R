find_flattening_point <- function(variance, threshold = 0.01) {
  # Calculate incremental increases in explained variance
  increments <- diff(variance)

  # Find the first point where the increment is below the threshold
  flattening_point <- which(increments < threshold)[1]

  # In case no flattening point is found, return the length of the variance vector
  if (is.na(flattening_point)) {
    return(length(variance))
  } else {
    return(flattening_point)
  }
}
find_elbow_point <- function(values) {
  # Creating a line from the first to the last point
  n_points <- length(values)
  line <- seq(values[1], values[n_points], length.out = n_points)

  # Calculating distances from the actual values to the line
  distances <- abs(values - line)

  # The index of the maximum distance is the "elbow"
  return(which.max(distances))
}
perform_pca_analysis <- function(DT) {
  # Standardizing the data
  DT_scaled <- scale(DT)

  # Performing PCA
  pca_result <- prcomp(DT_scaled, scale. = TRUE)
  explained_variance <- pca_result$sdev^2 / sum(pca_result$sdev^2)
  cum_var <- cumsum(explained_variance)

  # Initialize an empty list to store results
  results <- list()

  # Looping through variance thresholds from 0.5 to 0.9
  for (variance_threshold in seq(0.5, 0.9, by = 0.05)) {
    # Finding the number of components that explain at least the current threshold of variance
    num_components_variance <- which(cum_var >= variance_threshold)[1]
    results[[paste("VarianceThreshold", variance_threshold)]] <- num_components_variance
  }

  # Method 5: Kaiser’s Rule
  eigenvalues <- pca_result$sdev^2
  num_components_kaiser <- sum(eigenvalues > 1)
  results[["KaisersRule"]] <- num_components_kaiser

  # Method 4: Automated Elbow Detection in Scree Plot
  num_components_elbow <- find_elbow_point(pca_result$sdev)
  results[["ElbowMethod"]] <- num_components_elbow

  # Method 3: Automated Flattening Point Detection
  num_components_flattening <- find_flattening_point(cum_var)
  results[["FlatteningPointMethod"]] <- num_components_flattening

  # Returning the results
  return(results)
}
