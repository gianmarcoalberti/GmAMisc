#' R function for binary Logistic Regression internal validation
#'
#' The function allows to perform internal validation of a binary Logistic Regression model implementing most of the procedure described in:\cr
#' Arboretti Giancristofaro R, Salmaso L. "Model performance analysis and model validation in logistic regression". Statistica 2003(63): 375–396.\cr
#'
#' The procedure consists of the following steps:\cr
#' (1) the whole dataset is split into two random parts, a fitting (75 percent) and a validation (25 percent) portion;\cr
#' (2) the model is fitted on the fitting portion (i.e., its coefficients are computed considering only the observations in that portion) and its performance is evaluated on both the fitting and the validation portion, using AUC as performance measure;\cr
#' (3) the model's estimated coefficients, p-values, and the p-value of the Hosmer and Lemeshow test are stored;\cr
#' (4) steps 1-3 are repeated B times, eventually getting a fitting and validation distribution of the AUC values and of the HL test p-values, as well as a fitting distribution of the coefficients and of the associated p-values.
#' The AUC fitting distribution provides an estimate of the performance of the model in the population of all the theoretical fitting samples; the AUC validation distribution represents an estimate of the model’s performance on new and independent data.\cr
#'
#' The function returns:\cr
#' -a chart with boxplots representing the fitting distribution of the estimated model's coefficients;
#' coefficients' labels are flagged with an asterisk when the proportion of p-values smaller than 0.05 across the selected iterations is at least 95 percent;\cr
#' - a chart with boxplots representing the fitting and the validation distribution of the AUC value across the selected iterations.
#' for an example of the interpretation of the chart, see the aforementioned article, especially page 390-91;\cr
#' -a chart of the levels of the dependent variable plotted against the predicted probabilities (if the model has a high discriminatory power,
#' the two stripes of points will tend to be well separated, i.e. the positive outcome of the dependent variable will tend to cluster
#' around high values of the predicted probability, while the opposite will hold true for the negative outcome of the dependent variable);\cr
#' -a list containing:\cr
#' $overall.model.significance: statistics related to the overall model p-value and to its distribution across the selected iterations;\cr
#' $parameters.stability: statistics related to the stability of the estimated coefficients across the selected iterations;\cr
#' $p.values.stability: statistics related to the stability of the estimated p-values across the selected iterations;\cr
#' $AUCstatistics: statistics about the fitting and validation AUC distribution;\cr
#' $Hosmer-Lemeshow statistics: statistics about the fitting and validation distribution of the HL test p-values.\cr
#'
#' As for the abovementioned statistics:\cr
#' -full: statistic estimated on the full dataset;\cr
#' -median: median of the statistic across the selected iterations;\cr
#' -QRNG: interquartile range across the selected iterations;\cr
#' -QRNGoverMedian: ratio between the QRNG and the median, expressed as percentage;\cr
#' -min: minimum of the statistic across the selected iterations;\cr
#' -max: maximum of the statistic across the selected iterations;\cr
#' -percent_smaller_0.05: (only for $overall.model.significance, $p.values.stability, and $Hosmer-Lemeshow statistics):
#' proportion of times in which the p-values are smaller than 0.05; please notice that for the overall model significance and for the p-values stability it is desirable that
#' the percentage is at least 95percent, whereas for the HL test p-values it is indeed desirable that the proportion is not larger than 5percent (in line with the interpetation of the test p-value
#' which has to be NOT significant in order to hint at a good fit);\cr
#' -significant (only for $p.values.stability): asterisk indicating that the p-values of the corresponding coefficient resulted smaller than 0.05
#' in at least 95percent of the iterations.
#'
#' @param data: dataframe containing the dataset (Dependent Variable must be stored in the first column to the left).
#' @param fit: object returned from glm() function.
#' @param B: desired number of iterations (200 by default).
#' @param g: number of groups to be used for the Hosmer-Lemeshow test (10 by default).
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
#' @seealso \code{\link{logregr}} , \code{\link{aucadj}}
#'
modelvalid <- function(data, fit, B=200, g=10, oneplot=TRUE, excludeInterc=FALSE){

  #disable scientific notation
  options(scipen=999)

  #set up empty containers for the AUC training (=fitting) values and for the AUC testing (=validation) values
  auc.train <- vector(mode = "numeric", length = B)
  auc.test <- vector(mode = "numeric", length = B)

  #set up empy container for the overall p-value of the model fitted to the training (=fitting) partition across the B iterations
  pvalue.full <- vector(mode = "numeric", length=B)

  #set up empty containers for p-values of the HL test perofrmed on the training (=fitting) partition and on the testing (=validation) partition
  hl.train <- vector(mode = "numeric", length = B)
  hl.test <- vector(mode = "numeric", length = B)

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

  #calculate the p-value of the HL test on the full dataset
  hl.full <- HosmerLemeshowTest(fitted(fit), fit$y)$H$p.value

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

    #calculate the p-value of the HL test on the training (=fitting) portion
    hl.train[i] <- HosmerLemeshowTest(fit=train$pred.prob, obs=train[,1], ngr=g)$H$p.value

    #calculate the p-value of the HL test on the testing (=validation) portion
    hl.test[i] <- HosmerLemeshowTest(fit=test$pred.prob.back, obs=test[,1], ngr=g)$H$p.value

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

  #work out some statistics for the p-values of the HL test on the training (=fitting) partition
  HLtraining.df <- data.frame(full=round(hl.full,3),
                              median=round(median(hl.train),3),
                              QRNG=round(quantile(hl.train, 0.75) - quantile(hl.train, 0.25),3),
                              QRNGoverMedian=round(((quantile(hl.train, 0.75) - quantile(hl.train, 0.25)) / median(hl.train)) * 100,1),
                              min=round(min(hl.train),3),
                              max=round(max(hl.train),3),
                              percent_smaller_0.05=sum(hl.train < 0.05) / B)

  #work out some statistics for the p-values of the HL test on the testing (=validation) partition
  HLtesting.df <- data.frame(full="-",
                             median=round(median(hl.test),3),
                             QRNG=round(quantile(hl.test, 0.75) - quantile(hl.test, 0.25),3),
                             QRNGoverMedian=round(((quantile(hl.test, 0.75) - quantile(hl.test, 0.25)) / median(hl.test)) * 100,1),
                             min=round(min(hl.test),3),
                             max=round(max(hl.test),3),
                             percent_smaller_0.05=sum(hl.test < 0.05) / B)

  #bind the two above dataframes togheter
  HLglobal.df <- rbind(HLtraining.df, HLtesting.df)

  #give names to the binded dataframe's rows
  rownames(HLglobal.df) <- c("HLfitting", "HLvalidation")

  #set the layout of the plot visualization according to whether or not the parameter 'oneplot' is set to TRUE
  if(oneplot==TRUE){
    m <- rbind(c(1,2), c(3,3))
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

  #stripchart of the discriminatory power of the full model
  stripchart(model$fitted.values ~ model$y,
             method = "jitter",
             pch = 20,
             col="#00000088", cex = 0.7,
             xlim=c(0,1),
             xlab="Predicted probability",
             ylab="Levels of the Dependent Variable (y)",
             sub=paste0("AUC: ", round(auc.full,3)),
             main="Discriminatory power of model \n(fitted to the full sample)",
             cex.main=0.95,
             cex.sub=0.75,
             cex.axis=0.7)

  results <- list("overall.model.significance"=pvalue.full.df,
                  "parameters.stability"=par.estim.stab,
                  "p.values.stability"=pvalues.stab,
                  "AUCstatistics"=AUCglobal.df,
                  "Hosmer-Lemeshow statistics"=HLglobal.df)

  # restore the original graphical device's settings
  par(mfrow = c(1,1))

  return(results)
}
