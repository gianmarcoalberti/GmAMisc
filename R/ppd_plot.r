#' R function for plotting Posterior Probability Densities for Bayesian modeled 14C dates/parameters
#'
#' The function allows plot Posterior Probability Densities with a nice outlook thanks to 'ggplot2'.\cr
#' It takes as input a dataframe that must be organized as follows (it is rather easy to do that once the data have been exported from OxCal):\cr
#' -calendar dates (first column to the left);\cr
#' -posterior probabilities (second column);\cr
#' -grouping variables (third column), which could contain the names of the events of interest (e.g., phase 1 start, phase 1 end, phase 2 start, phase 2 end, etc).
#' @param data: dataframe containing the data as returned by the OxCal program.
#' @param lower: lower limit of the calendar date axis.
#' @param upper: upper limit of the calendar date axis; if the lower and upper parameters are not provided, the default values will be the earliest and latest calendar dates.
#' @param type: type of plot the user wishes to plot (a: curves outlined by a line; b: curves plotted as solid areas; c: combination of a and b).
#' @keywords Bayesian parameters dates radiocarbon posterior probability density
#' @export
#' @examples
#' data(radioc_data) #load a toy dataset
#' ppdPlot(radioc_data, type="a") #plot the Posterior Probability Densities for the phases' parameters
#' ppdPlot(radioc_data, -1000, 100, type="b") #plot the Posterior Probability Densities for the phases' parameters, setting different boundaries for the x-axis and using filled curves instead of simple outlines
#'
ppdPlot <- function(data,lower=min(data[,1]),upper=max(data[,1]),type) {
  if (type == "a") {
    a <- ggplot(data=data) + geom_line(aes(x=data[,1],y=data[,2],color=data[,3])) + scale_x_continuous(limits =c(lower,upper)) + xlab("Calibrated date") + ylab("Probability density") + labs(color="Events")
    print(a)
  } else {
    if (type == "b") {
      b <- ggplot(data=data) + geom_area(position="identity", aes(x=data[,1],y=data[,2],fill=data[,3]), alpha=0.5) + scale_x_continuous(limits =c(lower,upper)) + xlab("Calibrated date") + ylab("Probability density") + labs(fill="Events")
      print(b)
    } else {
      if (type == "c") {
        a <- ggplot(data=data) + geom_line(aes(x=data[,1],y=data[,2],color=data[,3])) + scale_x_continuous(limits =c(lower,upper)) + xlab("Calibrated date")
        c <- a + geom_area(position="identity",aes(x=data[,1],y=data[,2],fill=data[,3]), alpha=0.5) + xlab("Calibrated date") + ylab("Probability density") + labs(fill="Events") + guides(color=FALSE)
        print(c)
      }
    }
  }
}
