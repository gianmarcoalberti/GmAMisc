#' R function for binary Logistic Regression internal validation
#'
#' The function allows to perform internal validation of a binary Logistic Regression model implementing most of the procedure described in:\cr
#' Arboretti Giancristofaro R, Salmaso L. "Model performance analysis and model validation in logistic regression". Statistica 2003(63): 375–396.\cr
#'
#' The procedure consists of the following steps:\cr
#' (1) the whole dataset is split into two random parts, a fitting (75 percent) and a validation (25 percent) portion;\cr
#' (2) the model is fitted on the fitting portion (i.e., its coefficients are computed considering only the observations in that portion) and its performance is evaluated on both the fitting and the validation portion, using AUC as performance measure;\cr
#' (3) the model's estimated coefficients and p-values are stored;\cr
#' (4) steps 1-3 are repeated B times, eventually getting a fitting and validation distribution of the AUC values, and the fitting distribution of the coefficients and of the associated p-values.
#' The AUC fitting distribution provides an estimate of the performance of the model in the population of all the theoretical fitting samples; the AUC validation distribution represents an estimate of the model’s performance on new and independent data.\cr
#'
#' The function returns:\cr
#' -a chart with boxplots representing the fitting distribution of the estimated model's coefficients;
#' coefficients' labels are flagged with an asterisk when the proportion of p-values smaller than 0.05 across the selected iterations is at least 95 percent;\cr
#' - a chart with boxplots representing the fitting and the validation distribution of the AUC value across the selected iterations.
#' For an example of the interpretation of the chart, see the aforementioned article, especially page 390-91.\cr
#' -a list containing:\cr
#' $overall.model.significance: statistics related to the overall model p-value and to its distribution across the selected iterations;\cr
#' $parameters.stability: statistics related to the stability of the estimated coefficients across the selected iterations;\cr
#' $p.values.stability: statistics related to the stability of the estimated p-values across the selected iterations;\cr
#' $AUCstatistics: statistics about the fitting and validation AUC distribution.\cr
#'
#' As for the abovementioned statistics:\cr
#' -full: statistic estimated on the full dataset;\cr
#' -median: median of the statistic across the selected iterations;\cr
#' -QRNG: interquartile range across the selected iterations;\cr
#' -QRNGoverMedian: ratio between the QRNG and the median, expressed as percentage;\cr
#' -min: minimum of the statistic across the selected iterations;\cr
#' -max: maximum of the statistic across the selected iterations;\cr
#' -percent_smaller_0.05: (only for $overall.model.significance and $p.values.stability): proportion of times in which the p-values are smaller than 0.05;\cr
#' -significant (only for $p.values.stability): asterisk indicating that the p-values of the corresponding coefficient resulted smaller than 0.05
#' in at least 95percent of the iterations.
#'
#' @param data: dataframe containing the dataset (Dependent Variable must be stored in the first column to the left).
#' @param fit: object returned from glm() function.
#' @param B: desired number of iterations (200 by default).
#' @param oneplot: TRUE (default) is the user wants the charts returned in a single visualization.
#' @param excludeInterc: if set to TRUE, the chart showing the boxplots of the parameters distribution across the selected iteration will have y-axis
#' limits corresponding to the min and max of the parameters value; this allows better displaying the boxplots of the model parameters when they end up
#' showing up too much squeezed due to comparatively higher/lower values of the intercept. FALSE is default.
#' @keywords validation
#' @export
#' @examples
#' data(log_regr_data) # load the sample dataset
#' model <- glm(admit ~ gre + gpa + rank, data = log_regr_data, family = "binomial") # fit a logistic regression model, storing the results into an object called 'model'
#' res <- modelvalid(data=log_regr_data, fit=model, B=1000) # run the function, using 1000 iterations, and store the result in the 'res' object
#'
modelvalid <- function(data, fit, B=200, oneplot=TRUE, excludeInterc=FALSE){

  #disable scientific notation
  options(scipen=999)

  #set up empty containers for the AUC training (=fitting) values and for the AUC testing (=validation) values
  auc.train <- vector(mode = "numeric", length = B)
  auc.test <- vector(mode = "numeric", length = B)

  #set up empy container for the overall p-value of the model fitted to the training (=fitting) partition across the B iterations
  pvalue.full <- vector(mode = "numeric", length=B)

  #store the predicted probabilities based on the fitted model
  #this will be used to calculate the AUC of the model based on the full sample
  data$pred.prob.full <- fitted(fit)

  #calculate the number of coefficients in the model
  n.of.coeff <- length(fit$coefficients)

  #create a matrix to store the model coefficients that will be estimated on B training (fitting) samples
  coeff.matrix <- matrix(nrow=B, ncol=n.of.coeff)

  #give names to the columns of the matrix
  colnames(coeff.matrix) <- names(coefficients(fit))

  #create a matrix to store the coefficients' p values that will be estimated on B training (fitting) samples
  pvalues.matrix <- matrix(nrow=B, ncol=n.of.coeff)

  #give names to the columns of the matrix
  colnames(pvalues.matrix) <- names(coefficients(fit))

  #calculate the auc of the full model
  auc.full <- roc(data[,1], data$pred.prob.full, data=data)$auc

  #set the progress bar to be used inside the loop
  pb <- txtProgressBar(min = 0, max = B, style = 3)

  for(i in 1:B){
    #sample.split requires caTools
    #split the sample in a random part, corresponding to the 75percent of the full sample
    #observations in the input dataset are given TRUE or FALSE according to whether they are (randomly) selected as belonging to the
    #training or to the testing partition
    sample <- sample.split(data[,1], SplitRatio = .75)

    #assign to the training partition the observations flagged are TRUE
    train <- subset(data, sample == TRUE)

    #assign to the testing partition the observations flagged as FALSE
    test <- subset(data, sample == FALSE)

    #fit the model to the fitting partition
    fit.train <- glm(formula(fit), data = train, family = "binomial")

    #calculate the overall p-value for the model fitted to the fitting partition
    pvalue.full[i] <- anova(fit.train, update(fit.train, ~1), test="Chisq")$`Pr(>Chi)`[2]

    #save the estimated coefficients of the model fitted to the fitting partition
    coeff.matrix[i,] <- fit.train$coefficients

    #save the p-values of the coefficients of the model fitted to the fitting partition
    pvalues.matrix[i,] <- summary(fit.train)$coefficients[,4]

    #save the predicted probabilities of the model fitted to the fitting partition
    train$pred.prob <- fitted(fit.train)

    #save the probabilities predicted by applying to the validation partition
    #the parameters estimated on the fitting partition
    test$pred.prob.back <- predict.glm(fit.train, newdata=test, type="response")

    #roc requires pROC; calculate the AUC fitting values
    #the performance of the model fitted on the fitting partition is evaluated on the fitting portion itself
    #AND (see below) on the testing portion
    auc.train[i] <- roc(train[,1], train$pred.prob, data=train)$auc

    #calculate the AUC  validation values
    #the performance of the model fitted to the fitting portion (see above) is also evaluated on validation portion
    auc.test[i] <- roc(test[,1], test$pred.prob.back, data=test)$auc

    setTxtProgressBar(pb, i)
  }

  #unlist the values of the overall p-value calculated within the preceding loop
  pvalue.full <- unlist(pvalue.full)

  #work out some statistics for the overall p-value of the model fitted to the train (=fitting) partition
  pvalue.full.df <- data.frame(full=anova(fit, update(fit, ~1), test="Chisq")$`Pr(>Chi)`[2],
                               median=round(median(pvalue.full),6),
                               QRNG=round(quantile(pvalue.full, 0.75) - quantile(pvalue.full, 0.25),6),
                               QRNGoverMedian=round(((quantile(pvalue.full, 0.75) - quantile(pvalue.full, 0.25)) / median(pvalue.full)) * 100,1),
                               min=round(min(pvalue.full),6),
                               max=round(max(pvalue.full),6),
                               percent_smaller_0.05=sum(pvalue.full < 0.05) / B * 100)

  #attach name to the preceding dataframe
  rownames(pvalue.full.df) <- "overall p-value"

  #work out some statistics from the pvalues matrix
  pvalues.stab <- data.frame(full=round(summary(fit)$coefficients[,4],6),
                             median=round(apply(pvalues.matrix, 2, median),6),
                             QRNG=round(apply(pvalues.matrix, 2, function(x) quantile(x, 0.75)-quantile(x,0.25)),6),
                             QRNGoverMedian=round((abs(apply(pvalues.matrix, 2, function(x) quantile(x, 0.75) - quantile(x,0.25))) / abs(apply(pvalues.matrix, 2, median))) * 100,1),
                             min=round(apply(pvalues.matrix,2,min),6),
                             max=round(apply(pvalues.matrix,2,max),6),
                             percent_smaller_0.05=apply(pvalues.matrix,2, function(x) sum(x < 0.05)) / B * 100)

  #work out some statistics from the coefficients matrix
  par.estim.stab <- data.frame(full=round(coefficients(fit),3),
                               median=round(apply(coeff.matrix, 2, median),3),
                               QRNG=round(apply(coeff.matrix, 2, function(x) quantile(x, 0.75)-quantile(x,0.25)),3),
                               QRNGoverMedian= round((abs(apply(coeff.matrix, 2, function(x) quantile(x, 0.75) - quantile(x,0.25))) / abs(apply(coeff.matrix, 2, median))) * 100,1),
                               min=round(apply(coeff.matrix,2,min),3),
                               max=round(apply(coeff.matrix,2,max),3))

  #create a dataframe to store the parameters' labels and to calculate the proportion of p-values smaller than 0.05;
  #in other words, calculate the proportion of how many times each coefficients' p-value proves significant across the selected iterations
  p.valued.labls.df <- data.frame(labels=names(coefficients(fit)), prop=apply(pvalues.matrix,2, function(x) sum(x < 0.05)) / B)

  #if the proportion calculated in the above step is smaller than 0.05, flag the parameter with an asteriscs;
  #the starred labels will be used later on in the boxplot chart
  p.valued.labls.df$indicator <- ifelse(p.valued.labls.df$prop > 0.95, "*", "")

  #attach the asterisks to the dataframe with statistics about pvalues stability
  pvalues.stab$significant <- p.valued.labls.df$indicator

  #work out some statistics for the training (=fitting) AUC
  AUCtraining.df <- data.frame(full=round(auc.full,3),
                               median=round(median(auc.train),3),
                               QRNG=round(quantile(auc.train, 0.75) - quantile(auc.train, 0.25),3),
                               QRNGoverMedian=round(((quantile(auc.train, 0.75) - quantile(auc.train, 0.25)) / median(auc.train)) * 100,1),
                               min=round(min(auc.train),3),
                               max=round(max(auc.train),3))

  #work out some statistics for the testing (=validation) AUC
  AUCtesting.df <- data.frame(full="-",
                              median=round(median(auc.test),3),
                              QRNG=round(quantile(auc.test, 0.75) - quantile(auc.test, 0.25),3),
                              QRNGoverMedian=round(((quantile(auc.test, 0.75) - quantile(auc.test, 0.25)) / median(auc.test)) * 100,1),
                              min=round(min(auc.test),3),
                              max=round(max(auc.test),3))

  #bind the two above dataframes togheter
  AUCglobal.df <- rbind(AUCtraining.df,AUCtesting.df)

  #give names to the binded dataframe's rows
  rownames(AUCglobal.df) <- c("AUCfitting", "AUCvalidation")

  #set the layout of the plot visualization according to whether or not the parameter 'oneplot' is set to TRUE
  if(oneplot==TRUE){
    m <- rbind(c(1,2))
    layout(m)
  }

  #if exludeInterc is TRUE, set the y-axis limits to the min and max of the coefficient matrix values, excluding the first column,
  #which stores the Intercept values
  if(excludeInterc==TRUE){
    ylimlower <- min(coeff.matrix[,-c(1)])
    ylimupper <- max(coeff.matrix[,-c(1)])
  } else {
    ylimlower <- NULL
    ylimupper <- NULL
  }

  #boxplot of the randomized coefficients distribution
  boxplot(coeff.matrix,
          names=paste0(names(coefficients(fit)), p.valued.labls.df$indicator),
          main=paste0("Boxplots of parameters distribution across ", B, " iterations", "\nn= ", nrow(log_regr_data) * 0.75, " (75% of the full sample)"),
          sub="Starred parameters have p-value < 0.05 in at least 95% of the iterations",
          cex.main=0.95,
          cex.sub=0.75,
          cex.axis=0.7,
          ylim=c(ylimlower, ylimupper),
          las=2)

  abline(h=0, lty=2, col="grey")

  #boxplot of the training and testing AUC
  boxplot(auc.train, auc.test, names=c("AUC fitting", "AUC validation"), cex.axis=0.7)

  title(main=paste("Cross-validated AUC", "\nn of iterations:", B),
        sub=paste("AUC full sample=", AUCglobal.df$full[1], "\nAUC fitting min, median, max:", AUCglobal.df$min[1], AUCglobal.df$median[1], AUCglobal.df$max[1], "\nAUC validation min, median, max:", AUCglobal.df$min[2], AUCglobal.df$median[2], AUCglobal.df$max[2]),
        cex.main=0.95,
        cex.sub=0.75)


  results <- list("overall.model.significance"=pvalue.full.df,
                 "parameters.stability"=par.estim.stab,
                 "p.values.stability"=pvalues.stab,
                 "AUCstatistics"=AUCglobal.df)

  # restore the original graphical device's settings
  par(mfrow = c(1,1))

  return(results)
}
