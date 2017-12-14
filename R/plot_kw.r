#' R function for visually displaying Kruskal-Wallis test's results
#'
#' The function  allows allows to perform Kruskal-Wallis test, and to display the test's results in a plot along with boxplots.
#' The boxplots display the distribution of the values of the two samples, and jittered points represent the individual observations. At the bottom of the chart, the test statistics (H) is reported, along with the degrees of freedom and the associated p value.\cr
#' Setting the parameter 'posthoc' to TRUE, the Dunn's test is performed (with Bonferroni adjustment by default): a dot chart is returned, as well as a list of p-values (2-sided). In the dot chart, a RED line indicates the 0.05 threshold. The groups compared on a pairwise basis are indicated on the left-hand side of the chart.
#' @param x: object storing the values to be analysed.
#' @param y: object storing a grouping variable with 3 or more levels.
#' @param strip: logical value which takes FALSE (by default) or TRUE if the user wants jittered points to represent individual values.
#' @param notch: logical value which takes FALSE (by default) or TRUE if user does not or do want to have notched boxplots in the final display, respectively; it is worth noting that overlapping of notches indicates a not significant difference at about 95 percent confidence.
#' @param omm: (which stands for overall mean and median) takes FALSE (by default) or TRUE if user wants the mean and median of the overall sample plotted in the chart (as a dashed RED line and dotted BLUE line respectively).
#' @param outl: logical value which takes FALSE or TRUE (by default) if users want the boxplots to display outlying values.
#' @param posthoc: logical value which takes FALSE (default) or TRUE if user does not or does want to perform a follow-up test (namely, the Dunn's test) in order to locate which group significantly differs from the others.
#' @param adjust: sets the desidered method for p-values adjustment in the context of the Dunn's test; the list of methods is the following: Bonferroni ("bonferroni"; default); Holm ("holm"), Hochberg ("hochberg"), Hommel ("hommel"), Benjamini & Hochberg ("BH" or its alias "fdr"), Benjamini & Yekutieli ("BY"),  none ("none"). For more info, see the p.adjust's help documentation in R (?p.adjust).
#' @keywords Kruskal Wallis test nonparametric
#' @export
#' @examples
#' mydata <- data.frame(values=c(rnorm(30, 100,10),rnorm(30, 80,10),rnorm(30, 98,10)), group=c(rep("A", 30),rep("B", 30),rep("C", 30))) #create a toy dataset
#' kwPlot(x=mydata$values, y=mydata$group, strip=TRUE, notch=TRUE, posthoc=TRUE) # performs the test, displays the test's result, including jittered points and notches. It also performs the Dunn's posthoc test using Bonferroni p-value correction.
#'
kwPlot <- function (x, y, strip=FALSE, notch=FALSE, omm=FALSE, outl=TRUE, posthoc=FALSE, adjust="bonferroni"){
  options(scipen=999)
  data <- data.frame(value=x, group=y)
  H <- round(kruskal.test(data[,1] ~ data[,2])$statistic, 3)
  p.value <- round(kruskal.test(data[,1] ~ data[,2])$p.value, 6)
  p <- ifelse(p.value < 0.001, "< 0.001", ifelse(p.value < 0.01, "< 0.01", ifelse(p.value < 0.05, "< 0.05",round(p.value, 3))))
  degree.freed <- nlevels(data[,2])-1
  boxplot(data[,1] ~ data[,2], data = data, notch = notch, outline=outl)
  chart.title="Box Plots"
  if (strip==TRUE) {
  stripchart(data[,1] ~ data[,2], vertical = TRUE, data = data, method = "jitter", add = TRUE, pch = 16, col="#00000088", cex = 0.5)
    chart.title="Jittered Box Plots"
    } else {
    }
  title(main=chart.title, sub=paste("Kruskal-Wallis H=", H, ", df=", degree.freed ,", p=", p), cex.sub=0.8)
  if (omm==TRUE) {
    abline(h=mean(data[,1]), lty=2, col="red")
    abline(h=median(data[,1]), lty=3, col="blue")
  } else {
  }
  if (posthoc==TRUE) {
    res <- DunnTest(data[,1], data[,2], method=adjust) #requires 'DescTools' package; note: since version 0.99.17, the DunnTest function returns 2-sided p values by default
    res1 <- res[[1]]
    data.f <- data.frame(pair=row.names(res1), res1, row.names = NULL)
    dotchart2(data.f$pval, labels = data.f$pair, sort = FALSE, lty = 2, xlim = c(0, 1), xlab = paste0("2-sided p-value with ", adjust, " correction (red reference line set at 0.05)\n", "(Kruskal-Wallis H=", H, ", df=", degree.freed ,", p=", p,")")) #requires 'Hmisc' package
    title(main="Dunn's Post-Hoc Test")
    abline(v = 0.05, lty = 2, col = "RED")
    return(data.f)
    } else {
  }
  }
