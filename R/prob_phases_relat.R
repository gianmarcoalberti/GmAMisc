#' R function to calculate the Posterior Probability for different chronological relations between two Bayesian radiocarbon phases
#'
#' The function allows to calculate the posterior probability for different chronological relations between two phases defined via Bayesian radiocarbon modeling.
#' For the results to make sense, the phases have to be defined as independent if one wishes to assess what is the posterior probability for different relative chronological relations between them.\cr
#'
#' The rationale for this approach is made clear in an article by Buck et al 1992 (https://doi.org/10.1016/0305-4403(92)90025-X), and it runs as follows:
#' "if we do not make any assumption about the relationship between the phases, can we test how likely they are to be in any given order"?\cr
#'
#' Data can be fed into the function in two ways:\cr
#' -the function takes as input the table provided by the 'OxCal' program  as result of the 'Order' query.
#' Once the table as been saved from 'OxCal' in .csv format, you have to feed it in R. A .csv file can be imported into R using (for instance):
#' \eqn{mydata <- read.table(file.choose(), header=TRUE, sep=",", dec=".", as.is=T)};\cr
#'
#' be sure to insert the phases' parameters (i.e., the starting and ending boundaries of the two phases) in the OxCal's Order query in the following order:
#' StartA, EndA, StartB, EndB; that is, first the start and end of your first phase, then the start and end of the second one;
#' you can give any name to your phases, as long as the order is like the one described.\cr
#'
#' -alternatively, 8 relevant parameters (which can be read off from the Oxcal's Order query output) can be manually fed into the function (see the Arguments section above).\cr
#'
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
#' @param data: matrix containing the posterior probability of the chronological relation between the Starting and Ending boundaries of two independent phases, as returned by the OxCal's 'Order' query (see Details).
#' @param sAoldersB: probability of startA being older than startB.
#' @param sAoldereB: probability of startA being older than endB.
#' @param eAoldersB: probability of endA being older than startB.
#' @param eAoldereB: probability of endA being older than endB.
#' @param sBoldersA: probability of startB being older than startA.
#' @param sBoldereA: probability of startB being older than endA.
#' @param eBoldersA: probability of endB being older than startA.
#' @param eBoldereA: probability of endB being older than endA.
#' @param sort: logical which takes TRUE or FALSE (default) if the user does or does not want the returned posterior probabilities sorted in descending order.
#' @keywords radiocarbon
#' @export
#' @examples
#' data(phases) #load a toy dataset
#' res <- prob.phases.relat(phases) #calculate the Posterior Probability for the chronological relation between two phases, stores the results in the 'res' object, and produce a dot chart.
#'
#' res <- prob.phases.relat(data=NULL, sAoldersB=0.613, sAoldereB=1, eAoldersB=0.0010, eAoldereB=0.666, sBoldersA= 0.386, sBoldereA=0.999, eBoldersA=0.000039, eBoldereA=0.3334) # same as above, but manually feeding relevant parameters
#'
prob.phases.relat <- function(data=NULL, sAoldersB=NULL, sAoldereB=NULL, eAoldersB=NULL, eAoldereB=NULL, sBoldersA=NULL, sBoldereA=NULL, eBoldersA=NULL, eBoldereA=NULL, sort=FALSE) {

  if(is.null(data)==TRUE & is.null(sAoldersB)==FALSE){
    prob.data <- matrix(nrow=4, ncol=4)
    diag(prob.data) <- 0
    prob.data[1,2] <- 1
    prob.data[2,1] <- 0
    prob.data[3,4] <- 1
    prob.data[4,3] <- 0
    prob.data[1,3] <- sAoldersB
    prob.data[1,4] <- sAoldereB
    prob.data[2,3] <- eAoldersB
    prob.data[2,4] <- eAoldereB
    prob.data[3,1] <- sBoldersA
    prob.data[3,2] <- sBoldereA
    prob.data[4,1] <- eBoldersA
    prob.data[4,2] <- eBoldereA
    prob.data <- as.data.frame(prob.data)
    new.column <- c("sA", "eA", "sB", "eB")
    prob.data <- cbind(new.column, prob.data)
    data <- prob.data
  }

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
