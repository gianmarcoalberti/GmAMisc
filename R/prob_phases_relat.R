#' R function to calculate the Posterior Probability for different chronological relations between two Bayesian radiocarbon phases
#'
#' The function allows to calculate the posterior probability for different chronological relations between two phases defined via Bayesian radiocarbon modeling. For the results to make sense, the phases have to be defined as independent if one wishes to assess what is the posterior probability for different relative chronological relations between them.\cr
#' The rationale for this approach is made clear in an article by Buck et al 1992 (https://doi.org/10.1016/0305-4403(92)90025-X), and it runs as follows: "if we do not make any assumption about the relationship between the phases, can we test how likely they are to be in any given order"?\cr
#' The function takes as input the table provided by the 'OxCal' program  as result of the 'Order' query. Once the table as been saved from 'OxCal' in .csv format, you have to feed it in R. A .csv file can be imported into R using (for instance): mydata <- read.table(file.choose(), header=TRUE, sep=",", dec=".", as.is=T)\cr
#' Make sure to insert the phases' parameters (i.e., the starting and ending boundaries of the two phases) in the OxCal's Order query in the following order: StartA, EndA, StartB, EndB. That is, first the start and end of your first phase, then the start and end of the second one. You can give any name to your phases, as long as the order is like the one described.\cr
#' Given two phases A and B, the function allows to calculate the posterior probability for:\cr
#' -A being within B\cr
#' -B being within A\cr
#' -A starting earlier but overlapping with B\cr
#' -B starting earlier but overlapping with A\cr
#' -A being entirely before B\cr
#' -B being entirely before A\cr
#' -sA being within B\cr
#' -eA being within B\cr
#' -sB being within A\cr
#' -eB being within A\cr
#' where 's' and 'e' refer to the starting and ending boundaries of a phase. The function will return a table and a dot plot.\cr
#' Thanks are due to Dr. Andrew Millard (Durham University) for the help provided in working out the operations on probabilities.
#' @param data: data matrix containing the posterior probability of the chronological relation between the Starting and Ending boundaries of two independent phases, as returned by the OxCal's 'Order' query. Given the .csv file returnd by Oxcal, it can be imported in R using: mydata <- read.table(file.choose(), header=TRUE, sep=",", dec=".", as.is=T).
#' @param sort: logical which takes TRUE of FALSE (default) if the user does or does not want the returned posterior probabilities sorted in descending order.
#' @keywords Bayesian model dates radiocarbon posterior probability chronological relation phases
#' @export
#' @examples
#' data(phases) #load a toy dataset
#' res <- prob.phases.relat(phases) #calculate the Posterior Probability for the chronological relation between two phases, stores the results in the 'res' object, and produce a dot chart.
#'
prob.phases.relat <- function(data, sort=FALSE) {
  df <- as.data.frame(matrix(nrow = 10, ncol = 2))
  colnames(df)<-c("relation","posterior_prob")
  rel.labs <-c("A within B","B within A","A overlapping earlier than B","B overlapping earlier than A", "A entirely before B","B entirely before A", "sA within B", "eA within B", "sB within A", "eB within A")
  df[,1] <- rel.labs
  df[1,2] <- data[3,2]*data[2,5]
  df[2,2] <- data[1,4]*data[4,3]
  df[3,2] <- data[1,4]*data[3,3]*data[2,5]
  df[4,2] <- data[3,2]*data[1,5]*data[4,3]
  df[5,2] <- data[2,4]
  df[6,2] <- data[4,2]
  df[7,2] <- 1-(data[1,4]+(1-data[1,5]))
  df[8,2] <- 1-(data[2,4]+(1-data[2,5]))
  df[9,2] <- 1-(data[3,2]+(1-data[3,3]))
  df[10,2] <- 1-(data[4,2]+(1-data[4,3]))
  ifelse(sort == TRUE, df <- df[order(-df$posterior_prob), ], df <- df)
  dotchart2(as.numeric(df$posterior_prob), labels=df$relation, sort=FALSE, lty=2, xlim=c(0.0, 1.0), xlab="posterior probability")
  return(df)
  }
