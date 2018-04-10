#' R function for permutation-based chi-square test of independence
#'
#' The function performs the chi-square test of independence on the basis of permuted tables, whose number is selected by user.\cr
#'
#' For the rationale of this approach, see for instance the nice description provided by Beh E.J., Lombardo R. 2014, Correspondence Analysis: Theory, Practice and New Strategies, Chichester, Wiley, pages 62-64.\cr
#'
#' The function produces:\cr
#' (1) a chart that displays the permuted distribution of the chi square statistic based on B permuted tables. The selected number of permuted tables, the observed chi square, the 95th percentile of the permuted distribution, and the associated p value are reported at the bottom of the chart;\cr
#'
#' (2) a chart that displays the bootstrap distribution of Cramer's V coefficient, based on a number of bootstrap replicates which is equal to the value of the function's parameter B;\cr
#'
#' (3) a chart that the Pearson's Standardized Residuals: a colour scale allows to easily understand which residual is smaller (BLUE) or larger (RED) than expected under the hypothesis of independence. Should the user want to only display residuals larger than a given threshold, it suffices to set the filter parameter to TRUE, and to specify the desidered threshold by means of the thresh parameter, which is set at 1.96 by default.
#' @param data: dataframe containing the input contingency table.
#' @param B: desired number of permuted tables (999 by default).
#' @param resid: TRUE or FALSE (default) if the user does or doesn't want to plot the table of Pearson's standardized residuals.
#' @param filter: takes TRUE or FALSE (default) if the user does or does't want to filter the Pearson's standardized residuals according to the threshold provided by the thresh parameter; by default, the threshold is set at 1.96, which corresponds to an alpha level of 0.05.
#' @param thresh: value of the standardized residuals below which the residuals will be not displayed (by default, the threshold is set at 1.96, which corresponds to an alpha level of 0.05).
#' @param cramer: takes TRUE or FALSE (default) if the user does or doesn't want to calculate and plot the bootstrap confidence interval for Cramer's V.
#' @keywords chi-square permutation independence cramer residuals categorical association
#' @export
#' @examples
#' data(assemblage)
#' chiperm(data=assemblage, resid=TRUE, cramer=TRUE)
#'
chiperm <- function(data, B=999, resid=FALSE, filter=FALSE, thresh=1.96, cramer=FALSE){
  options(warn=-1)

  rowTotals <- rowSums(data)
  colTotals <- colSums(data)

  obs.chi.value <- chisq.test(data)$statistic

  chistat.perm <- vector(mode = "numeric", length = B)
  chi.statistic <- function(x)  chisq.test(x)$statistic
  chistat.perm <- sapply(r2dtable(B, rowTotals, colTotals), chi.statistic)

  p.lowertail <- (1 + sum (chistat.perm < obs.chi.value)) / (1 + B)
  p.uppertail <- (1 + sum (chistat.perm > obs.chi.value)) / (1 + B)
  two.sided.p <- 2 * min(p.lowertail, p.uppertail)

  p.to.report <- ifelse(two.sided.p < 0.001, "< 0.001",
                        ifelse(two.sided.p < 0.01, "< 0.01",
                               ifelse(two.sided.p < 0.05, "< 0.05",
                                      round(two.sided.p, 3))))

  d <- density(chistat.perm)

  plot(d, main="Chi-square statistic Permuted Distribution",
       sub=paste0("\nSolid line: observed chi-sq (", round(obs.chi.value, 3), ")","\nDashed line: 95th percentile of the permuted chi-sq (=alpha 0.05 threshold) (", round(quantile(chistat.perm, c(0.95)),3), ")", "\np value: ", p.to.report, " (n. of permutations: ", B,")"),
       xlab = "",
       cex.main=0.85,
       cex.sub=0.70)
  polygon(d, col = "#BCD2EE88", border = "blue")
  rug(chistat.perm, col = "#0000FF")
  abline(v = obs.chi.value)
  abline(v = round(quantile(chistat.perm, c(0.95)), 5), lty = 2, col = "blue")

  if (resid==TRUE) {
    res <- chisq.test(data)
    col1 <- colorRampPalette(c("blue", "white", "red"))
    ifelse(filter==FALSE, res$stdres,res$stdres[abs(res$stdres) <= thresh] <- 0)
    corrplot(res$stdres, method="circle", addCoef.col="black", is.corr=FALSE, cl.lim = c(min(res$stdres), max(res$stdres)),tl.col="black", tl.cex=0.8, col = col1(100)) #requires 'corrplot'
  } else {}

  if (cramer==TRUE) {
    cramerv <- vector (mode = "numeric", length = B)
    obs.cramerv[1] <- cramersV(data) #requires 'lsr'
    for(i in 2:B){
      cramerv[i] <- cramersV(contingency.data.break(data, boot=TRUE))
    }

    d.cramer <- density(cramerv)
    perc5 <- quantile(cramerv, c(0.05))
    perc95 <- quantile(cramerv, c(0.95))
    plot(d.cramer, xlab="", main="Cramer's V Bootstrap Distribution", sub=paste0("\nCramer's V: ", round(cramerv[1], 3), "\nDashed lines: 5th and 95th percentile", "\n95% Condifence Interval: ", round(perc5, 3), "-", round(perc95, 3), " (n. of bootstrap replicates: ",B,")"), cex.sub=0.8)
    polygon(d.cramer, col = "blue", border = "red")
    rug(cramerv, col = "#0000FF")
    abline(v = perc5, lty = 2, col = "red")
    abline(v = perc95, lty = 2, col = "red")
  } else {}
}
