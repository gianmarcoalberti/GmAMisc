#' R function for binary Logistic Regression internal validation
#'
#' The function allows to perform internal validation of a binary Logistic Regression model, implementing part of the procedure described in
#' Arboretti Giancristofaro R, Salmaso L. "Model performance analysis and model validation in logistic regression". Statistica 2003(63): 375–396.\cr
#' The procedure consists of the following steps:\cr
#' (1) the whole dataset is split into two random parts, a fitting (75 percent) and a validation (25 percent) portion;\cr
#' (2) the model is fitted on the fitting portion (i.e., its coefficients are computed considering only the observations in that portion) and its performance is evaluated on both the fitting and the validation portion, using AUC as performance measure;\cr
#' (3) steps 1-2 are repeated B times, eventually getting a fitting and validation distribution of the AUC values. The former provides an estimate of the performance of the model in the population of all the theoretical training samples; the latter represents an estimate of the model’s performance on new and independent data.\cr
#' The function returns two boxplots that represent the training and the testing (i.e., validation) distribution of the AUC value across the 1000 iterations. For an example of the interpretation of the chart, see the aforementioned article, especially page 390-91.
#' @param data: dataframe containing the dataset (Dependent Variable must be stored in the first column to the left).
#' @param fit: object returned from glm() function.
#' @param B: desired number of iterations (see description of the procedure above).
#' @keywords logistic regression validation bootstrap
#' @export
#' @examples
#' data(log_regr_data) # load the sample dataset
#' model <- glm(admit ~ gre + gpa + rank, data = log_regr_data, family = "binomial") # fit a logistic regression model, storing the results into an object called 'model'
#' modelvalid(data=log_regr_data, fit=model, B=1000) # run the function, using 1000 iterations
#'
modelvalid <- function(data, fit, B){
  auc.train <- vector(mode = "numeric", length = B)
  auc.test <- vector(mode = "numeric", length = B)
  data$pred.prob.full <- fitted(fit)
  auc.full <- roc(data[,1], data$pred.prob.full, data=data)$auc
  for(i in 1:B){
    sample <- sample.split(data[,1], SplitRatio = .75) #require caTools
    train <- subset(data, sample == TRUE)
    test <- subset(data, sample == FALSE)
    fit.train <- glm(formula(fit), data = train, family = "binomial")
    train$pred.prob <- fitted(fit.train)
    auc.train[i] <- roc(train[,1], train$pred.prob, data=train)$auc #require pROC
    test$pred.prob.back <- predict.glm(fit.train, newdata=test, type="response")
    auc.test[i] <- roc(test[,1], test$pred.prob.back, data=test)$auc
  }
  auc.train.min <- round(min(auc.train), digits=4)
  auc.train.max <- round(max(auc.train), digits=4)
  auc.train.median <- round(median(auc.train), digits=4)
  auc.test.min <- round(min(auc.test), digits=4)
  auc.test.max <- round(max(auc.test), digits=4)
  auc.test.median <- round(median(auc.test), digits=4)
  boxplot(auc.train, auc.test, names=c("AUC training", "AUC testing"))
  title(main=paste("Cross-validated AUC", "\nn of iterations:", B), sub=paste("AUC full sample=", round(auc.full, digits=4), "\nAUC training min, median, max:", auc.train.min, auc.train.median, auc.train.max, "\nAUC testing min, median, max:", auc.test.min, auc.test.median, auc.test.max), cex.sub=0.8)
}
