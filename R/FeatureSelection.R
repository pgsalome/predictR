# FeatureSelection.R

#' Boruta Feature Selection
#'
#' Applies the Boruta feature selection algorithm to the dataset.
#'
#' @param DT Data frame containing the features and outcome.
#' @param CS Censoring status.
#' @param S Survival time or binary indicator for survival.
#' @importFrom survival Surv
#' @importFrom Boruta Boruta getSelectedAttributes getImpFerns
#' @return A data frame of selected features and their counts.
#' @export
boruta_feature_selection <- function(DT, CS, S) {
  features_boruta <- data.frame(f=character(), count=integer(), stringsAsFactors=FALSE)


  DT$CS <- CS
  DT$S <- S
  if (S[1]=='binary'){
    boruta <- Boruta(CS ~ ., data = DT, pValue=0.5, mcAdj = TRUE, doTrace = 0,  holdHistory = FALSE,maxRuns = 500)
  } else {
    boruta <- Boruta(survival::Surv(S, CS) ~ ., data = DT, pValue=0.5, mcAdj = TRUE, doTrace = 0,
                     holdHistory = FALSE,maxRuns = 500)
  }

  features1 <- getSelectedAttributes(boruta, withTentative = FALSE)
  features1 <- features1[!is.na(features1)]


  return(features1)
}



#' Univariate Feature Selection for Generalized Linear Models
#'
#' Performs univariate feature selection using generalized linear models.
#' @param DT Data frame containing the features and group.
#' @param group The group variable for analysis.
#' @param p The p-value for significance testing.
#' @param NF Number of features to select.
#' @return A list containing the selected features.
#' @export
#'
get_featuresUnivGlm <- function(DT, group, p, NF) {
  covariates <- colnames(DT)
  features_all <- data.frame(f=character(), p=numeric())
  for (i in 1:ncol(DT)) {
    x <- covariates[i]
    e <- if (is.factor(DT[[x]])) {
      as.formula(paste('group ~ factor(', x, ')'))
    } else {
      as.formula(paste('group ~', x))
    }
    fit <- glm(e, data=DT)
    if (coef(summary(fit))[, 4][2] < p) {
      row <- data.frame(f=covariates[i], p=coef(summary(fit))[, 4][2])
      features_all <- rbind(features_all, row)
    }
  }
  x <- features_all[order(-features_all$p), ][1:min(NF, nrow(features_all)), ]
  features <- x$f
  res <- list(features = features)
  return(res)
}

#' Univariate Feature Selection for Cox Proportional Hazards Model
#'
#' Performs univariate feature selection using Cox proportional hazards model.
#' @param DT Data frame containing the features and surv object.
#' @param surv Survival object (created using Surv function).
#' @param p The p-value for significance testing.
#' @param NF Number of features to select.
#' @return A list containing the selected features.
#' @export
get_featuresUnivCox <- function(DT, surv, p, NF) {
  covariates <- colnames(DT)
  univ_formulas <- sapply(covariates, function(x) {
    if (is.factor(DT[[x]])) {
      as.formula(paste('surv ~ factor(', x, ')'))
    } else {
      as.formula(paste('surv ~', x))
    }
  })
  univ_models <- lapply(univ_formulas, function(x) { coxph(x, data = DT) })
  univ_results <- lapply(univ_models, function(x) {
    x <- summary(x)
    p.value <- signif(x$waldtest["pvalue"], digits = 2)
    beta <- signif(x$coef[1], digits = 2)
    HR <- signif(x$coef[2], digits = 5)
    HR.confint.lower <- signif(x$conf.int[, "lower .95"], 2)
    HR.confint.upper <- signif(x$conf.int[, "upper .95"], 2)
    wald.test <- paste0(signif(x$wald["test"], digits = 2),
                        " (", HR.confint.lower, "-", HR.confint.upper, ")")
    res <- c(beta, HR, wald.test, p.value)
    names(res) <- c("beta", "HR", "wald.test", "p.value")
    return(res)
  })
  res <- t(as.data.frame(univ_results, check.names = FALSE))
  univariate_results <- as.data.frame(res)
  filtered_uv_results <- univariate_results[as.numeric(as.character(univariate_results$p.value)) < p, ]

  # Print HR to verify
  HR <- as.numeric(as.character(filtered_uv_results$HR))

  # Calculating Inverse or same HR
  IHR <- sapply(HR, function(x) if (x < 1) 1 / x else x)
  sorted_uv_results <- sort(IHR, index.return = TRUE, decreasing = TRUE)
  features <- names(filtered_uv_results$beta[sorted_uv_results$ix[1:NF]])
  res <- list(fs = filtered_uv_results, features = features)
  return(res)
}

#' Unified Univariate Feature Selection
#'
#' Chooses the correct univariate feature selection method based on the outcome type.
#' @param DT Data frame containing the features and outcome.
#' @param S Survival time or binary indicator for survival.
#' @param CS Censoring status.
#' @importFrom survival Surv
#' @return A data frame of selected features and their counts.
#' @export
univariate_feature_selection <- function(DT, S, CS) {
  features_uv <- data.frame(f=character(), count=integer(), stringsAsFactors=FALSE)
  if (S[1] == 'binary') {
    filtered_uv_results <- try(get_featuresUnivGlm(DT, CS, p=0.05, NF=nrow(DT)/10), TRUE)
  } else {
    surv <- Surv(S, CS)
    filtered_uv_results <- get_featuresUnivCox(DT, surv, p=0.05, NF=nrow(DT)/10)
  }
  if (!inherits(filtered_uv_results, "try-error")) {
    features0 <- as.character(filtered_uv_results$features)
    features0 <- features0[!is.na(features0)]

  }
  return(features0)
}




#' Minimum Redundancy Maximum Relevance Feature Selection
#'
#' Performs feature selection using the Minimum Redundancy Maximum Relevance (mRMR) method.
#' mRMR selects features that are highly relevant to the target response while also being minimally redundant with each other. This method is particularly useful in high-dimensional data scenarios, such as genomics and radiomics.
#'
#' @param DT Data frame containing the features and the response variable.
#' @param response_col Name of the column in DT that contains the response variable.
#' @return A vector of selected feature names based on the mRMR criteria.
#' @importFrom mRMRe mRMR.data mRMR.ensemble solutions scores set.thread.count
#' @importFrom survival Surv
#' @examples
#' # Assuming DT is your data frame and 'response' is the response variable
#' selected_features <- mrmr_feature_selection(DT, "response")
#' @export
mrmr_feature_selection <- function(DT, CS, S, N=10, solution_count=100, run_parallel=TRUE) {

  if (S[1] == 'binary') {
    DT <-  data.frame(target=CS, DT)
    data <- mRMR.data(data = DT)

  } else {
    DT <- data.frame(target=Surv(S, CS), DT)
    data <- mRMR.data(data = DT )

  }
  total_features <- ncol(DT) - 1

  max_feature_count <- 1  # Start with the lowest possible value

  # Iteratively find the maximum feasible feature_count
  for (fc in 1:total_features) {
    if (choose(total_features, fc) >= solution_count) {
      max_feature_count <- fc
    } else {
      next
    }
  }
  if (run_parallel) {
    # Set the number of threads for parallel computation
    numThreads <- parallel::detectCores()
    set.thread.count(numThreads)
  }
  filter <- mRMR.ensemble(data = data, target_indices = 1,
                feature_count = max_feature_count , solution_count = solution_count)

  if (run_parallel) {
    # Reset thread count to default (optional)
    set.thread.count(1)
  }
  # Assuming 'filter' is your mRMR filter object
  solution_matrices <- solutions(filter)
  solution_scores <- scores(filter)
  # Initialize a vector to store the cumulative scores for each feature
  cumulative_scores <- numeric(ncol(DT))
  # Initialize a vector to count the occurrences of each feature
  feature_counts <- numeric(ncol(DT))

  # Iterate through each solution matrix and its corresponding scores matrix
  for (i in 1:length(solution_matrices)) {
    sol_matrix <- solution_matrices[[i]]
    score_matrix <- solution_scores[[i]]

    # Iterate through each feature in the solution
    for (j in 1:ncol(sol_matrix)) {
      feature_indices <- sol_matrix[, j]
      feature_scores <- score_matrix[, j]

      # Update cumulative scores and counts
      for (k in 1:length(feature_indices)) {
        feature_idx <- feature_indices[k]
        cumulative_scores[feature_idx] <- cumulative_scores[feature_idx] + feature_scores[k]
        feature_counts[feature_idx] <- feature_counts[feature_idx] + 1
      }
    }
  }

  # Calculate the average score for each feature
  avg_scores <- ifelse(feature_counts > 0, cumulative_scores / feature_counts, 0)
  # Rank the features based on average scores
  ranked_features <- rank(-avg_scores)  # Negative sign for descending order
  # Select the top 'nrow(DT)/N' features
  top_feature_indices <- which(ranked_features <= nrow(DT)/N)
  # Get the feature names
  selected_features <- names(DT)[top_feature_indices]

  return(selected_features)
}


#' Random Forest Feature Selection
#'
#' This function implements the Random Forest Importance method for feature selection. It utilizes the Random Forest algorithm to assess the importance of each feature in the dataset. The function is designed to handle both binary and survival outcome data.
#'
#' @param DT A data frame containing the features along with survival time (`S`) and censoring status (`CS`).
#' @param CS A vector indicating the censoring status. If your data is survival data, this parameter should be provided. For binary classification, it can be ignored.
#' @param S A vector indicating the survival time or the binary outcome. Survival time is expected if `CS` is provided, otherwise a binary outcome is assumed.
#'
#' @return A vector of selected feature names ranked by their importance.
#'
#' @importFrom randomForestSRC var.select
#' @importFrom doParallel registerDoParallel
#'
#' @examples
#' # For binary classification
#' selected_features_binary <- random_forest_feature_selection(DT = your_data, S = your_binary_outcome)
#'
#' # For survival data
#' selected_features_survival <- random_forest_feature_selection(DT = your_data, CS = your_censoring_status, S = your_survival_time)
#'
#' @export

rf_feature_selection <- function(DT, CS, S,run_paralell = TRUE) {

  # Initialize list to store different sets of features
  features_rf <- list()
  if (run_paralell) {
    # Register parallel backend
    options(rf.cores = parallel::detectCores(), mc.cores = parallel::detectCores())
  }
  # Check if the outcome is binary or survival
  if (S[1] == 'binary') {
    # Binary outcome
    # Ensuring S is part of the data frame
    DT$Outcome <- as.factor(S)  # Ensure the outcome is treated as a factor for classification
    var_select_result_md <- var.select(Outcome ~ ., data = DT, method = "md")
    #var_select_result_vh <- var.select(Outcome ~ ., data = DT, method = "vh")
    var_select_result_vhvimp <- var.select(Outcome ~ ., data = DT, method = "vh.vimp")

    #rf_model <- rfsrc(Outcome ~ ., data = DT, ntree = 100, na.action = "na.omit",parallel = parallel)
  } else {
    # Survival outcome
    # Ensuring S and CS are part of the data frame
    DT$Time <- S
    DT$Status <- CS
    #rf_model <- rfsrc(Surv(Time, Status) ~ ., data = DT, ntree = ntree,
    #                  na.action = "na.omit", parallel = parallel)

    var_select_result_md <- var.select(Surv(Time, Status) ~ ., data = DT, method = "md")
    #var_select_result_vh <- var.select(Surv(Time, Status) ~ ., data = DT, method = "vh")
    var_select_result_vhvimp <- var.select(Surv(Time, Status) ~ ., data = DT, method = "vh.vimp")
  }

  features_rf[["rfmd"]] <- var_select_result_md$topvars
  #features_rf[["rfvh"]] <- var_select_result_vh$topvars
  features_rf[["rfvhvimp"]] <- var_select_result_vhvimp$topvars

  if (run_paralell) {
    stopImplicitCluster()  # Stop parallel backend
  }
  return(features_rf)
}


#' XGBoost Feature Selection
#'
#' Performs feature selection using XGBoost and selects the top 1/N features.
#'
#' @param DT Data frame containing the features and outcome.
#' @param CS Censoring status.
#' @param S Survival time or binary indicator for survival.
#' @param nrounds Number of boosting rounds.
#' @param N The divisor to determine the top percentage of features to select.
#' @return A data frame of top selected features based on their importance scores.
#' @importFrom xgboost xgb.DMatrix xgboost xgb.importance
#' @export
xgboost_feature_selection <- function(DT, CS, S, nrounds = 100, N = 10) {
  DT <- DT[, !(colnames(DT) %in% c("S", "CS")), drop = FALSE]

  if (S[1] == 'binary') {
    dtrain <- xgb.DMatrix(data = as.matrix(DT), label = S)
    objective <- "binary:logistic"
  } else {
    survival_labels <- ifelse(CS == 1, -S, S)
    dtrain <- xgb.DMatrix(data = as.matrix(DT), label = survival_labels)
    objective <- "survival:cox"
  }

  params <- list(
    objective = objective,
    booster = "gbtree",
    eval_metric = "logloss"
  )

  bst <- xgboost(params = params, data = dtrain, nrounds = nrounds, verbose = 0)
  importance_matrix <- xgb.importance(model = bst)
  importance_matrix <- importance_matrix[order(-importance_matrix$Gain),]

  # Initialize list to store different sets of features
  features_xgboost <- list()

  # Select top 1/N features
  num_features_fixed <- max(1, floor(nrow(DT) / N))
  features_xgboost[["xgbTopN"]] <- head(importance_matrix, num_features_fixed)$Feature

  # Percentage-based feature selection
  percentages <- seq(1, 10, by = 1)
  for (pct in percentages) {
    num_features_pct <- max(1, floor(nrow(DT) * pct / 100))
    feature_set_name <- paste0("xgb", pct, "Pct")
    features_xgboost[[feature_set_name]] <- head(importance_matrix, num_features_pct)$Feature
  }

  return(features_xgboost)
}


#' Lasso Feature Selection
#'
#' Performs Lasso feature selection on the given dataset.
#'
#' @param DT Data frame containing the features and outcome.
#' @param CS Censoring status.
#' @param S Survival time or binary indicator for survival.
#' @return A data frame of selected features and their counts.
#' @export
lasso_feature_selection <- function(DT, CS, S, cv_iter=100,alpha=1,run_parallel=TRUE) {

  if (S[1]=='binary'){
    measure = "auc"
    family  = "binomial"
    y = S[1]
  } else {
    measure = "deviance"
    family  = "cox"
    y <- Surv(S, CS)
  }
  glmnetSTRUCT = list()
  x <- model.matrix(y~.,data=DT)
  glmnetSTRUCT$glmnet = lasso_cv(x=x, y=y, cv_iter=cv_iter,pen=NULL,alpha=alpha, measure=measure, family=family,run_parallel=run_parallel)
  if (!typeof(glmnetSTRUCT$glmnet) == "character"){
    # x = variables, y = outcome
    ## Store selected features for reference: lambda.min is the value of lambda that gives minimum mean-cross-validated error

    glmnetSTRUCT$selection = coef(glmnetSTRUCT$glmnet$glmnet.fit,s=glmnetSTRUCT$glmnet$lambda.min) # selected features
    # Use sparse matrix properties to find non-zero elements
    non_zero_indices = glmnetSTRUCT$selection@i + 1  # Adding 1 because R is 1-indexed
    non_zero_values = glmnetSTRUCT$selection@x

    # Extracting feature names
    feature_names = rownames(glmnetSTRUCT$selection)
    non_zero_feature_names = feature_names[non_zero_indices]

    # Saving the non-zero feature names
    features2 = non_zero_feature_names
    } else {
    message( "No Lasso features were returned")
    features2 = NULL
    }

  return(features2)
}
#' Extract Significant Features from PCA
#'
#' This function extracts significant features from the PCA result based on the specified number of components and threshold. It can optionally ignore a specific column during the extraction process.
#'
#' @param pca_result The result of PCA computation using prcomp or similar methods.
#' @param num_components The number of principal components to consider for feature extraction.
#' @param threshold The threshold for determining significance of features based on their loadings.
#' @param ignore_column Optional; the name of the column to be ignored during feature extraction, defaults to "sub".
#' @return A vector of significant feature names.
#' @examples
#' pca_res <- prcomp(DT, scale. = TRUE)
#' significant_features <- extract_significant_features(pca_res, num_components = 5, threshold = 0.3)
#' @export
extract_significant_features <- function(pca_result, num_components, threshold, ignore_column = "sub") {
  loadings <- abs(pca_result$rotation[, 1:num_components])
  if (ignore_column %in% rownames(loadings)) {
    loadings <- loadings[-which(rownames(loadings) == ignore_column), ]
  }
  significant <- apply(loadings, 1, function(x) any(x > threshold))
  return(rownames(loadings)[significant])
}

#' Find the flattening point in a vector of variance values
#'
#' This function calculates the flattening point in a vector of variance values
#' based on the specified threshold.
#'
#' @param variance A numeric vector of variance values.
#' @param threshold A numeric threshold defining the minimum incremental increase
#'                 in variance to be considered as a non-flattening point (default: 0.01).
#'
#' @return An integer indicating the flattening point index.
#'
#' @examples
#' variance_values <- c(0.1, 0.2, 0.3, 0.31, 0.32, 0.33, 0.35, 0.4)
#' find_flattening_point(variance_values, threshold = 0.05)
#'
#' @export
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
#' Find the elbow point in a vector of values
#'
#' This function calculates the elbow point in a vector of numeric values.
#'
#' @param values A numeric vector of values.
#'
#' @return An integer indicating the index of the elbow point.
#'
#' @examples
#' data_values <- c(10, 15, 20, 25, 30, 35, 40, 45, 50)
#' find_elbow_point(data_values)
#'
#' @export
find_elbow_point <- function(values) {
  # Creating a line from the first to the last point
  n_points <- length(values)
  line <- seq(values[1], values[n_points], length.out = n_points)

  # Calculating distances from the actual values to the line
  distances <- abs(values - line)

  # The index of the maximum distance is the "elbow"
  return(which.max(distances))
}


#' Perform PCA Feature Selection
#'
#' This function performs PCA on the given dataset and extracts significant features using various methods including Kaiser's rule, the elbow method, and flattening point method, at different significance thresholds.
#'
#' @param DT Data frame containing the features for PCA.
#' @return A list of significant features for each method and threshold.
#' @examples
#' pca_features <- perform_pca_feature_selection(DT)
#' @export
#'
perform_pca_feature_selection <- function(DT,apply_varimax=FALSE) {
  pca_result <- prcomp(DT, scale. = TRUE)
  varimax_append = ""
  # Apply Varimax rotation if requested
  if (apply_varimax) {
    rotation <- varimax(pca_result$rotation)
    pca_result$rotation <- rotation$loadings
    varimax_append = "varimax"
  }

  explained_variance <- pca_result$sdev^2 / sum(pca_result$sdev^2)
  cum_var <- cumsum(explained_variance)

  num_components_kaiser <- sum(pca_result$sdev^2 > 1)
  num_components_elbow <- find_elbow_point(pca_result$sdev)
  num_components_flattening <- find_flattening_point(cum_var)

  features_pca <- list()

  significance_thresholds <- seq(0.15, 0.3, by = 0.05)
  variance_thresholds <- c(0.85, 0.9, 0.95)

  for (threshold in significance_thresholds) {
    threshold_label <- as.integer(threshold * 100)

    # Kaiser Rule
    features_pca[[paste0("pcaKaiser",varimax_append, sprintf("%02d", threshold_label))]] <- extract_significant_features(pca_result, num_components_kaiser, threshold)

    # Elbow Method
    features_pca[[paste0("pcaElbow",varimax_append, sprintf("%02d", threshold_label))]] <- extract_significant_features(pca_result, num_components_elbow, threshold)

    # Flattening Point Method
    features_pca[[paste0("pcaFlattening", varimax_append, sprintf("%02d", threshold_label))]] <- extract_significant_features(pca_result, num_components_flattening, threshold)

    # Variance Thresholds
    for (variance_threshold in variance_thresholds) {
      num_components_variance <- which(cum_var >= variance_threshold)[1]
      features_pca[[paste0("pcaVariance",varimax_append, variance_threshold * 100, sprintf("%02d", threshold_label))]] <- extract_significant_features(pca_result, num_components_variance, threshold)
    }
  }

  return(features_pca)
}


#' Feature Selection Loop
#'
#' This function executes a feature selection loop for a given dataset. It applies univariate feature selection, Boruta feature selection, and Lasso feature selection iteratively for a specified number of iterations.
#'
#' @param DT Data frame containing the features and outcome variables.
#' @param outputDir Directory to save the output files.
#' @return The function does not return a value but saves feature selection results to files in the specified output directory.
#' @export
#' @importFrom caret createDataPartition
#
feature_selection_loop <- function(DT, outputDir, nrep=1000,run_parallel = TRUE) {
  S <- DT$S
  CS <- DT$CS
  DT$sub <- NULL
  DT$S <- NULL
  DT$CS <- NULL

  # Check if the iterations file exists and read the last iteration count
  iter_file <- file.path(outputDir, "iterations.rds")
  if(file.exists(iter_file)) {
    start_iter <- readRDS(iter_file) + 1
    message(paste("Resuming from iteration", start_iter))
  } else {
    start_iter <- 1
  }
  if (start_iter>nrep)
    return(NULL)
  # Initialize separate data frames for each method and a combined data frame
  features_uv <- data.frame(f = character(), count = integer(), stringsAsFactors = FALSE)
  features_boruta <- data.frame(f = character(), count = integer(), stringsAsFactors = FALSE)
  features_lasso <- data.frame(f = character(), count = integer(), stringsAsFactors = FALSE)
  features_enet <- data.frame(f = character(), count = integer(), stringsAsFactors = FALSE)
  features_mrmr <- data.frame(f = character(), count = integer(), stringsAsFactors = FALSE)

  for(j in start_iter:nrep) {

    message(paste('Fold #', j))
    train_ind <- createDataPartition(CS, times = 1, p = 0.9, list = FALSE)

    message(paste('Running the univariate feature selection #', j))
    selected_uv <- univariate_feature_selection(DT[train_ind, ], S[train_ind], CS[train_ind])
    features_uv <- update_features_df(features_uv, selected_uv)

    message(paste('Running the boruta feature selection #', j))
    selected_boruta <- boruta_feature_selection(DT[train_ind, ], CS[train_ind], S[train_ind])
    features_boruta <- update_features_df(features_boruta, selected_boruta)

    message(paste('Running the lasso feature selection #', j))
    selected_lasso <- lasso_feature_selection(DT[train_ind, ], CS[train_ind], S[train_ind],alpha=1,run_parallel =run_parallel)
    features_lasso <- update_features_df(features_lasso, selected_lasso)

    message(paste('Running the elastic net feature selection #', j))
    selected_enet <- lasso_feature_selection(DT[train_ind, ], CS[train_ind], S[train_ind],alpha=0.5,run_parallel =run_parallel)
    features_enet <- update_features_df(features_enet, selected_enet)

    message(paste('Running the Minimum Redundancy Maximum Relevance feature selection #', j))
    selected_mrmr <- mrmr_feature_selection(DT[train_ind, ], CS[train_ind], S[train_ind],run_parallel =run_parallel)
    features_mrmr <- update_features_df(features_mrmr, selected_mrmr)

    message(paste('Running xgboost feature selection #', j))
    selected_xgboost <- xgboost_feature_selection(DT[train_ind, ], CS[train_ind], S[train_ind])
    if (j == start_iter)
      features_xgboost <- setNames(
        replicate(length(selected_xgboost), data.frame(f = character(), count = integer(), stringsAsFactors = FALSE), simplify = FALSE),
        names(selected_xgboost)
      )
    features_xgboost <- update_features_df_list(features_xgboost, selected_xgboost)

    message(paste('Running rf feature selection #', j))

    selected_rf <- rf_feature_selection(DT[train_ind, ], CS[train_ind], S[train_ind])
    if (j == start_iter )
    features_rf <- setNames(
      replicate(length(selected_rf ), data.frame(f = character(), count = integer(), stringsAsFactors = FALSE), simplify = FALSE),
      names(selected_rf)
    )
    features_rf <- update_features_df_list(features_rf, selected_rf)

    message(paste('Running the PCA feature selection #', j))
    selected_pca <- perform_pca_feature_selection(DT[train_ind, ])
    if (j == start_iter)
      features_pca <- setNames(
        replicate(length(selected_pca), data.frame(f = character(), count = integer(), stringsAsFactors = FALSE), simplify = FALSE),
        names(selected_pca)
      )
    features_pca <- update_features_df_list(features_pca, selected_pca)

    message(paste('Running the PCA feature selection with varimax #', j))
    selected_pcavarimax <- perform_pca_feature_selection(DT[train_ind, ],apply_varimax=TRUE)
    if (j == start_iter)
      features_pcavarimax<- setNames(
        replicate(length(selected_pcavarimax), data.frame(f = character(), count = integer(), stringsAsFactors = FALSE), simplify = FALSE),
        names(selected_pcavarimax)
      )

    features_pcavarimax<- update_features_df_list(features_pcavarimax, selected_pcavarimax)


    feature_data <- list(
      features_uv = features_uv,
      features_boruta = features_boruta,
      features_lasso = features_lasso,
      features_enet = features_enet,
      features_rf = features_rf,
      features_xgboost = features_xgboost,
      features_pca = features_pca,
      features_pcavarimax = features_pcavarimax,
      features_mrmr = features_mrmr)

    save_feature_files(feature_data, outputDir)
    saveRDS(j, iter_file)

  }
  message("Done")
}

#' Stepwise Variable Selection in Cox Proportional Hazards Model
#'
#' Performs stepwise variable selection (both forward and backward) in a Cox proportional hazards model.
#' It handles right-censored data (using Time and Status) or interval-censored data (using T1, T2, and Status).
#' The function returns the names of the variables selected by the stepwise procedure.
#'
#' @param Time The time variable name for right-censored data (time to event).
#' @param T1 The start time variable name for interval-censored data.
#' @param T2 The end time variable name for interval-censored data.
#' @param Status The event status variable name.
#' @param variable.list A vector of candidate variable names to be considered in the model.
#' @param in.variable A vector of variable names to always include in the model. Defaults to "NULL".
#' @param data The data frame containing the variables.
#' @param sle The significance level for entry in the model during stepwise selection. Defaults to 0.15.
#' @param sls The significance level for staying in the model during stepwise selection. Defaults to 0.15.
#' @param vif.threshold The threshold for the variance inflation factor (VIF) to assess multicollinearity. Defaults to 999.
#' @return A vector of the names of the selected variables.
#' @importFrom survival coxph
#' @importFrom car vif
#' @importFrom MASS stepAIC
#' @export

stepwise_coxph <- function(Time = NULL, T1 = NULL, T2 = NULL, Status = NULL,
                           variable.list, in.variable = "NULL", data,
                           sle = 0.15, sls = 0.15, vif.threshold = 999) {

  # Determine the appropriate survival object based on input
  surv_object <- if (!is.null(Time) && !is.null(Status)) {
    Surv(data[[Time]], data[[Status]])
  } else if (!is.null(T1) && !is.null(T2) && !is.null(Status)) {
    Surv(data[[T1]], data[[T2]], data[[Status]])
  } else {
    stop("Appropriate time and status variables not provided.")
  }

  # Initialize the variable list for the model
  model_vars <- if (in.variable != "NULL") in.variable else c()

  # Stepwise selection
  final_model <- stepAIC(coxph(surv_object ~ ., data = data[, c(model_vars, variable.list), drop = FALSE]),
                         scope = list(upper = ~., lower = ~1),
                         direction = "both", trace = FALSE)

  # Extract selected variable names
  selected_vars <- all.vars(formula(final_model))

  # Removing 'Surv' and status variables from the list
  features_stepwise <- setdiff(selected_vars, c("Surv", Time, T1, T2, Status))

  return(features_stepwise)
}
