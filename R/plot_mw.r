#' R function for visually displaying Mann-Whitney test's results
#'
#' The function  allows to perform Mann-Whitney test, and to display the test's results in a plot along with two boxplots. For information about the test, and on what it is actually testing, see for instance the interesting article by R M Conroy, "What hypotheses do "nonparametric" two-group tests actually test?", in The Stata Journal 12 (2012): 1-9.\cr
#' The returned boxplots display the distribution of the values of the two samples, and jittered points represent the individual observations. At the bottom of the chart, a subtitle arranged on three lines reports relevant statistics:\cr
#' -test statistic (namely, U) and the associated z and p value;\cr
#' -Probability of Superiority value (which can be interpreted as an effect-size measure, as discussed in: https://nickredfern.wordpress.com/2011/05/12/the-mann-whitney-u-test/);\cr
#' -another measure of effect size, namely r (see https://stats.stackexchange.com/questions/124501/mann-whitney-u-test-confidence-interval-for-effect-size), whose thresholds are indicated in the last line of the plot's subtitle.\cr
#' The function may also return a density plot (coupled with a rug plot at the botton of the same chart) that displays the distribution of the pairwise differences between the values of the two samples being compared. The median of this distribution (which is represented by a blue reference line in the same chart) corresponds to the Hodges-Lehmann estimator.
#' @param x: object storing the values of the first group being compared.
#' @param y: object storing either the values of the second group being compared or a grouping variable with 2 levels.
#' @param xlabl: if y is not a grouping variable, user may want to specify here the name of the x group that will show up in the returned boxplots (default is "x").
#' @param ylabl: if y is not a grouping variable, user may want to specify here the name of the y group that will show up in the returned boxplots (default is "y").
#' @param strip: logical value which takes FALSE (by default) or TRUE if the user wants jittered points to represent individual values.
#' @param notch: logical value which takes FALSE (by default) or TRUE if user does not or do want to have notched boxplots in the final display, respectively; it is worth noting that overlapping of notches indicates a not significant difference at about 95 percent confidence.
#' @param omm: (which stands for overall mean and median) takes FALSE (by default) or TRUE if user wants the mean and median of the overall sample plotted in the chart (as a dashed RED line and dotted BLUE line respectively).
#' @param outl: logical value which takes FALSE or TRUE (by default) if users want the boxplots to display outlying values.
#' @param HL: logical value that takes TRUE or FALSE (default) if the user wants to diplay the distribution of the pairwise differences between the values of the two samples being compared; the median of that distribution is the Hodges-Lehmann estimator.
#' @keywords Mann Whitney test nonparametric
#' @export
#' @examples
#' mydata <- data.frame(values=c(rnorm(30, 100,10),rnorm(30, 80,10)), group=c(rep("A", 30),rep("B", 30))) #create a toy dataset
#' mwPlot(x=mydata$values, y=mydata$group, strip=TRUE, omm=TRUE, notch=TRUE, HL=TRUE) # performs the test, displays the test's result, including jittered points, notches, overall median and mean, and the Hodges-Lehmann estimator
#' mwPlot(x=rnorm(30,80,10), y=rnorm(30, 60,10), xlabl="A", ylabl="B", strip=TRUE) # performs the test when y is not a grouping variable but a vector of values
#'
mwPlot <- function (x,y,xlabl="x",ylabl="y", strip=FALSE,notch=FALSE,omm=FALSE, outl=TRUE, HL=FALSE){
  options(scipen=999)
  if (is.numeric(y)==FALSE) {
    data <- data.frame(value=x, group=y)
  } else {data <- data.frame(value=c(x,y), group=c(rep(xlabl, length(x)), rep(ylabl, length(y))))
  }
  res <- wilcox.test(data[,1] ~ data[,2], conf.int=TRUE)
  U <- wilcox.test(data[,1] ~ data[,2])$statistic
  p <- ifelse(res$p.value < 0.001, "< 0.001", ifelse(res$p.value < 0.01, "< 0.01", ifelse(res$p.value < 0.05, "< 0.05",round(res$p.value, 3))))
  print(paste("p-value=",res$p.value))
  samples.size <- count (data[,2]) #requires the plyr package
  PS <- round(U/(samples.size[1,2] * samples.size[2,2]), 3)
  z <- round(wilcox_test(data[,1] ~ data[,2])@statistic@teststatistic,3) #requires the coin package
  r <- round(abs(z/sqrt(samples.size[1,2] + samples.size[2,2])), 3)
  boxplot(data[,1] ~ data[,2], data = data, notch = notch, outline=outl)
  chart.title="Box Plots"
  if (strip==TRUE) {
    stripchart(data[,1] ~ data[,2], vertical = TRUE, data = data, method = "jitter", add = TRUE, pch = 16, col="#00000088", cex = 0.5)
    chart.title="Jittered Box Plots"
  } else {
  }
  title(main=chart.title, sub=paste("Mann-Whitney U=", U, ", z=",z, ", p=", p, "; Probability of Superiority=", PS, "; r=", r, "\nP{value(group to the left) > value(group to the right)}=", PS,"\nEffect size thresholds [r]: small (0.10), medium (0.30), large (0.50)"), cex.sub=0.8)
  if (omm==TRUE) {
    abline(h=mean(data[,1]), lty=2, col="red")
    abline(h=median(data[,1]), lty=3, col="blue")
  } else {
  }
  if (HL==TRUE) {
    unstacked.data <- unstack(data)
    diff <- outer(unstacked.data[[1]], unstacked.data[[2]],"-")
    m <- round(median(diff), 3)
    plot(density(diff), main="Pairwise differences distribution", xlab="", sub=paste("difference in location", "\nmedian (Hodges-Lehmann estimator):", m))
    polygon(density(diff), col="grey")
    rug(diff, col="red")
    abline(v=m, lty=2, col="blue")
  } else {
  }
}
