#' R function for optimism-adjusted AUC (internal validation)
#'
#' The function allows to calculate the AUC of a (binary) Logistic Regression model, adjusted for optimism.\cr
#'
#' The function performs an internal validation of a model via a bootstrap procedure (devised by Harrell and colleagues), which enables to estimate the degree of optimism of a fitted model and the extent to which the model will be able to generalize outside the training dataset.
#' If you want more info, you can refer to this website (http://thestatsgeek.com/2014/10/04/adjusting-for-optimismoverfitting-in-measures-of-predictive-ability-using-bootstrapping/), and/or read the following interesting article (in which the bootstrap procedure is described at page 776):
#' http://thestatsgeek.com/2014/10/04/adjusting-for-optimismoverfitting-in-measures-of-predictive-ability-using-bootstrapping/\cr
#'
#' The returned boxplots represent:\cr
#' -the distribution of the AUC value in the bootstrap sample (auc.boot), which represents "an estimation of the apparent performance" (according to the aforementioned reference);\cr
#' -the distribution of the AUC value deriving from the model fitted to the bootstrap samples and evaluated on the original sample (auc.orig), which represents the model performance on independent data.\cr
#' At the bottom of the chart, the apparent AUC (i.e., the value deriving from the model fitted to the original dataset) and the AUC adjusted for optimism are reported.
#' @param data: dataframe containing the dataset (note: the Dependent Variable must be stored in the first column to the left).
#' @param fit: object returned from glm() function.
#' @param B: desired number of bootstrap resamples (suggested values: 100 or 200).
#' @keywords logistic regression model validation AUC optimism
#' @export
#' @examples
#' data(log_regr_data) # load the sample dataset
#' model <- glm(admit ~ gre + gpa + rank, data = log_regr_data, family = "binomial") # fit a logistic regression model, storing the results into an object called 'model'
#' aucadj(data=log_regr_data, fit=model, B=200)
#' @seealso \code{\link{logregr}} , \code{\link{modelvalid}}
#'
aucadj <- function(data, fit, B){
  fit.model <- fit
  data$pred.prob <- fitted(fit.model)
  auc.app <- roc(data[,1], data$pred.prob, data=data)$auc                              # require 'pROC'
  auc.boot <- vector (mode = "numeric", length = B)
  auc.orig <- vector (mode = "numeric", length = B)
  o <- vector (mode = "numeric", length = B)
  pb <- txtProgressBar(min = 0, max = B, style = 3)                                    #set the progress bar to be used inside the loop
  for(i in 1:B){
    boot.sample <- sample.rows(data, nrow(data), replace=TRUE)                         # require 'kimisc'
    fit.boot <- glm(formula(fit.model), data = boot.sample, family = "binomial")
    boot.sample$pred.prob <- fitted(fit.boot)
    auc.boot[i] <- roc(boot.sample[,1], boot.sample$pred.prob, data=boot.sample)$auc
    data$pred.prob.back <- predict.glm(fit.boot, newdata=data, type="response")
    auc.orig[i] <- roc(data[,1], data$pred.prob.back, data=data)$auc
    o[i] <- auc.boot[i] - auc.orig[i]
    setTxtProgressBar(pb, i)
    }
  auc.adj <- auc.app - (sum(o)/B)
  boxplot(auc.boot, auc.orig, names=c("auc.boot", "auc.orig"))
  title(main=paste("Optimism-adjusted AUC", "\nn of bootstrap resamples:", B), sub=paste("auc.app (blue line)=", round(auc.app, digits=4),"\nadj.auc (red line)=", round(auc.adj, digits=4)), cex.sub=0.8)
  abline(h=auc.app, col="blue", lty=2)
  abline(h=auc.adj, col="red", lty=3)
}
