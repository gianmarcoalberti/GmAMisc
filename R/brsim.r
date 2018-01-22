#' R function for Brainerd-Robinson simiarity coefficient
#'
#' The function allows to calculate the Brainerd-Robinson similarity coefficient, taking as input a cross-tabulation (dataframe).\cr
#'
#' The function returns:\cr
#' a) a correlation matrix in tabular form;\cr
#' b) a heat-map representing, in a graphical form, the aforementioned correlation matrix.\cr
#'
#' In the heat-map (which is built using the 'corrplot' package), the size and the color of the squares are proportional to the Brainerd-Robinson coefficients, which are also reported by numbers.\cr
#'
#' In order to "penalize" BR similarity coefficient(s) arising from assemblages with unshared categories, the function does what follows: it divides the BR coefficient(s) by the number of unshared categories plus 0.5. The latter addition is simply a means to be still able to penalize coefficient(s) arising from assemblages having just one unshared category. Also note that joint absences will have no weight on the penalization of the coefficient(s). In case of assemblages sharing all their categories, the corrected coefficient(s) turns out to be equal to the uncorrected one.
#' @param data: dataframe containing the dataset (note: assemblages in rows, variables in columns).
#' @param which: takes "rows" (default) if the user wants the coefficients be calculated for the row categories, "cols" if the users wants the coefficients be calculated for the column categories.
#' @param correction: takes FALSE (default) if the user does not want the coefficients to be corrected, while TRUE will provide corrected coefficients.
#' @param rescale: takes FALSE if the user does NOT want the coefficients to be rescaled between 0.0 and 1.0 (i.e., the user will get the original version of the Brainerd-Robinson coefficient (spanning from 0 [maximum dissimilarity] to 200 [maximum similarity]), while TRUE (default) will return rescaled coefficient.
#' @keywords logistic regression model validation AUC optimism
#' @export
#' @examples
#' data(assemblage)
#' coeff <- BRsim(data=assemblage, correction=FALSE, rescale=TRUE)
#'
BRsim <- function(data, which="rows", correction=FALSE, rescale=TRUE) {
  x <- data
  ifelse(which=="rows", x <- x, x <- t(x))
  rd <- dim(x)[1]
  results <- matrix(0,rd,rd)
  if (correction == T){
    for (s1 in 1:rd) {
      for (s2 in 1:rd) {
        zero.categ.a <-length(which(x[s1,]==0))
        zero.categ.b <-length(which(x[s2,]==0))
        joint.absence <-sum(colSums(rbind(x[s1,], x[s2,])) == 0)
        if(zero.categ.a==zero.categ.b) {
          divisor.final <- 1
        } else {
          divisor.final <- max(zero.categ.a, zero.categ.b)-joint.absence+0.5
        }
        results[s1,s2] <- round((1 - (sum(abs(x[s1, ] / sum(x[s1,]) - x[s2, ] / sum(x[s2,]))))/2)/divisor.final, digits=3)
      }
    }
  } else {
    for (s1 in 1:rd) {
      for (s2 in 1:rd) {
        results[s1,s2] <- round(1 - (sum(abs(x[s1, ] / sum(x[s1,]) - x[s2, ] / sum(x[s2,]))))/2, digits=3)
      }
    }
  }
  rownames(results) <- rownames(x)
  colnames(results) <- rownames(x)
  col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white", "cyan", "#007FFF", "blue", "#00007F"))
  if (rescale == F) {
    upper <- 200
    results <- results * 200
  } else {
    upper <- 1.0
  }
  corrplot(results, method="square", addCoef.col="red", is.corr=FALSE, cl.lim = c(0, upper), col = col1(100), tl.col="black", tl.cex=0.8)
  return(results)
}
