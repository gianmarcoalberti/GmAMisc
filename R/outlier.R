#' R function for univariate outliers detection
#'
#' The function allows to perform univariate outliers detection using three different methods. These methods are those described in:\cr
#' Wilcox R R, "Fundamentals of Modern Statistical Methods: Substantially Improving Power and Accuracy", Springer 2010 (2nd edition), pages 31-35.\cr
#' Two of the three methods are robust, and are therefore less prone to the masking effect.\cr
#' (1) With the mean-based method, an observation is considered outlier if the absolute difference between that observation and the sample mean is more than 2 Standard Deviations away (in either direction) from the mean. In the plot returned by the function, the central reference line is indicating the mean value, while the other two are set at mean-2*SD and mean+2*SD.\cr
#' (2) The median-based method considers an observation as being outlier if the absolute difference between the observation and the sample median is larger than the Median Absolute Deviation divided by 0.6745. In this case, the central reference line is set at the median, while the other two are set at median-2*MAD/0.6745 and median+2*MAD/0.6745.\cr
#' (3) The boxplot-based method considers an observation as being an outlier if it is either smaller than the 1st Quartile minus 1.5 times the InterQuartile Range, or larger than the 3rd Quartile minus 1.5 times the InterQuartile Range. In the plot, the central reference line is set at the median, while the other two are set at 1Q-1.5*IQR and 3Q+1.5*IQR.\cr
#' The function also returns a list containing information about the choosen method, the mid-point, lower and upper boundaries where non-outlying observations are expected to fall, total number of outlying observations, and a dataframe listing the observations and indicating which is considered outlier. In the charts, the outlying observations are flagged with their ID number.
#' @param x: vector storing the data.
#' @param method: outliers identification method, either "mean" (default), "median", or "boxplot".
#' @param addthres: takes FALSE or TRUE (default) if user does not want or does want some threshold lines be added to the returned chart.
#' @keywords outlier identification univariate robust statistics
#' @export
#' @examples
#' mydata <- c(2,3,4,5,6,7,8,9,50,50) # create a toy dataset
#' outlier(mydata, method="median", addthres=TRUE) # locate outlier(s) using the median-based method
#'
outlier <- function (x,method="mean",addthres=TRUE){
  if (method=="mean") {
    avrg <- mean(x)
    stdev <-sd(x)
    dtf <- data.frame(ID=seq.int(length(x)), obs=x, outlier=abs(x-avrg)>2*stdev)
    midp <- avrg
    lower <- avrg-2*stdev
    upper <- avrg+2*stdev
    outliern <- length(which(dtf=="TRUE"))
  } else {}
  if (method=="median") {
    med <- median(x)
    MAD <-median(abs(med-x))
    dtf <- data.frame(ID=seq.int(length(x)), obs=x, outlier=abs(x-med)>2*(MAD/0.6745))
    midp <- med
    lower <- med-2*(MAD/0.6745)
    upper <- med+2*(MAD/0.6745)
    outliern <- length(which(dtf=="TRUE"))
    } else {}
  if (method=="boxplot") {
    Q1 <- quantile(x, 0.25)
    Q3 <- quantile(x, 0.75)
    IntQ <-Q3-Q1
    dtf <- data.frame(ID=seq.int(length(x)), obs=x, outlier=x<Q1-1.5*IntQ | x>Q3+1.5*IntQ)
    midp <- median(x)
    lower <- Q1-1.5*IntQ
    upper <- Q3+1.5*IntQ
    outliern <- length(which(dtf=="TRUE"))
    } else {}
  if (addthres==TRUE) {
    p <- ggplot(dtf, aes(x=ID, y=obs, label=ID)) + geom_point(aes(colour=outlier)) + geom_text_repel(data = subset(dtf, outlier=="TRUE"), aes(label = ID), size = 2.7, colour="black", box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines")) + labs(x=paste0("observation ID number\n number of outliers detected=", outliern, "\n(outlier detection method=", method, ")"), y="observation value") + geom_hline(yintercept = midp, colour="black", linetype = "longdash") + geom_hline(yintercept = lower, colour="black", linetype = "longdash") + geom_hline(yintercept = upper, colour="black", linetype = "longdash")
    } else {
  p <- ggplot(dtf, aes(x=ID, y=obs, label=ID)) + geom_point(aes(colour=outlier)) + geom_text_repel(data = subset(dtf, outlier=="TRUE"), aes(label = ID), size = 2.7, colour="black", box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines")) + labs(x=paste0("observation ID number\n number of outliers detected=", outliern, "\n(outlier detection method=", method, ")"), y="observation value")
  }
  print(p)
  results <- (list(method=method, midpoint=midp, lowerbound=lower, upperbound=upper, outlN=outliern, flaggedData=dtf))
  return(results)
}
