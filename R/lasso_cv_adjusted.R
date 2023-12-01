#LASSO

#' Repeated Cross-Validation for LASSO using glmnet
#'
#' This function performs repeated cross-validation for LASSO models using the glmnet package.
#' It aims to find the optimal lambda with the smallest error. The function returns the lambda
#' values at minimum error and 1 standard error from the minimum.
#'
#' @param x Matrix of predictors.
#' @param y Response variable.
#' @param cv_iter Number of repetitions for cross-validation.
#' @param alpha Alpha parameter for glmnet.
#' @param measure Type of measure to use for cross-validation (default is 'deviance').
#' @param family Model family for glmnet.
#' @param pen Optional penalty factor for glmnet.
#' @importFrom glmnet cv.glmnet

#' @return A list containing the optimal lambda.min, lambda.1se, and the fitted glmnet model.
#' @examples
#' # Assuming x is your predictor matrix and surv is your response variable
#' cv_results <- glmnet.repcv(x, surv, nreps = 200, alpha = 1, nvar = 50)
## cross-validation, which is intended to stabilize the method and to find the optimal lambda with the smallest error.
#' @export
lasso_cv_adjusted <- function(x, y, cv_iter, alpha, pen=NULL, measure = c("mse", "deviance", "class", "auc", "mae"),
                         family=c("gaussian", "binomial", "poisson", "multinomial", "cox", "mgaussian")) {
  count = 0
  best = 500
  for (j in 1:50){
    print(paste("dfmax:",j,sep=""))
    CVMs <- NULL
    CVSDs <- NULL
    dfmax = j -1#set the maximum number of degrees-of-freedom
    #measure = match.arg(measure)
    lambda.mins <- NULL #list of minimum error lambdas
    lambda.1ses <- NULL #list of 1se lambdas

    for (i in 1:cv_iter) {
      if (is.null(pen))
        #perform 1 cross-validation (default: 10 folds)
        cv <- cv.glmnet(x, y, alpha=alpha, family=family, type.measure=measure, dfmax=dfmax, type.multinomial="grouped")
        #plot(cv, xvar = "dev")
      else
        ##perform 1 cross-validation (default: 10 folds)
        cv <- cv.glmnet(x, y, alpha=alpha, family=family, type.measure=measure, penalty.factor=pen, dfmax=dfmax, type.multinomial="grouped")
      CVMs  <- qpcR:::cbind.na(CVMs,  cv$cvm) #store partial likelihood deviances
      CVSDs <- qpcR:::cbind.na(CVSDs, cv$cvsd)#store standard errors
      lambda.mins = c(lambda.mins,cv$lambda.min) #store lambda.min for each run
      lambda.1ses = c(lambda.1ses,cv$lambda.1se) #store lambda.1se for each run
    }
    rownames(CVMs) <- cv$glmnet.fit$lambda[1:nrow(CVMs)]
    colnames(CVMs) <- NULL
    mean.CVMs = rowMeans(CVMs,na.rm=TRUE)
    mean.CVSDs = rowMeans(CVSDs,na.rm=TRUE)

    if (measure=="auc") {
      lambda.min.index = which.max(mean.CVMs) #max mean AUC
    } else {
      lambda.min.index = which.min(mean.CVMs)
    } #min mean error
    lambda.min <- as.numeric(names(lambda.min.index))

    if (measure=="auc") {
      ## which.max returns the first TRUE
      lambda.1se.index = which.max(mean.CVMs > max(mean.CVMs)-mean.CVSDs[lambda.min.index])
    } else {
      ## which.max returns the first TRUE
      lambda.1se.index = which.max(mean.CVMs < min(mean.CVMs)+mean.CVSDs[lambda.min.index])
    }
    lambda.1se <- as.numeric(names(lambda.1se.index))
    if (best > min(mean.CVMs)) {
      print(best)
      best = min(mean.CVMs)
      count = 0
    } else {
      count = count +1
    }

    if (count == 5)
    {
      print(dfmax-2)
      return(list(lambda.min = lambda.min,lambda.1se = lambda.1se, glmnet.fit = cv$glmnet.fit))
    }
  }

}

