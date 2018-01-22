#' R function for Nearest Neighbor analysis of point patterns
#'
#' The function allows to perform the Nearest Neighbor analysis of point patterns to formally test for the presence of clustering, overdispersion, or random spatial arrangement. Significance assessed via a randomized approach\cr
#'
#' The function uses a randomized approach to test the significance of the Nearest Neighbor distance: the observed average NN distance is compared against the distribution of average NN distances computed across B iterations.
#' In each iteration, a set of random points (with a sample size equal to the number of points of the input feature) is drawn.\cr
#'
#' The function produces a density chart of the randomized average NN distances, with a reference line indicating the observed average NN and a black dot representing the average of the randomized NN distances.
#' P-values are reported at the bottom of the same chart. The two tails of the randomized distribution are given a red (left tail, indicating the area of significant clustering) and a blue color (right tail, indicating the area of significant overdispersion).\cr
#'
#' The function also returns a list storing the following:\cr
#' -$obs.aver.NN.dist\cr
#' -$rnd.aver.NN.dist\cr
#' -$Prob. of obs. aver. NN dist. <= random. aver. NN dist.\cr
#' -$Prob. of obs. aver. NN dist. >= random. aver. NN dist.\cr
#' @param feature: feature dataset (of point type).
#' @param studyplot: shapefile (of polygon type) representing the study area; if not provided, the study area is internally worked out as the convex hull enclosing the input feature dataset.
#' @param buffer: add a buffer to the studyplot (0 by default); the unit depends upon the units of the input data.
#' @param B: number of randomizations to be used (1000 by default).
#' @param addmap: TRUE (default) or FALSE if the user wants or does not want a map of the study area and of feature dataset to be also displayed.
#' @keywords NNA
#' @export
#' @examples
#' data(springs)
#' res <- NNa(springs) #perform the analysis using all default values; the result points to a significant clustering
#'
#' data(rndpoints)
#' res <- NNa(rndpoints, buffer=100, B=500) #perform the analysis using the 'rndpoints' dataset, add a 100m buffer arounf the points' convexhull, and use 500 iterations instead of 1000; the result points to a random arrangement
#'
NNa <- function(feature, studyplot=NULL, buffer=0, B=1000, addmap=TRUE){
  if(is.null(studyplot)==TRUE){
    ch <- rgeos::gConvexHull(feature)
    region <- rgeos::gBuffer(ch, width=buffer)
  } else {
    region <- studyplot
  }
  matr <- rgeos::gDistance(feature, byid=TRUE)
  diag(matr) <- NA
  NNdist <- numeric(B+1)
  NNdist[1] <- mean(apply(matr, 2, min, na.rm=TRUE))
  pb <- txtProgressBar(min = 0, max = B, style = 3)                             #set the progress bar to be used inside the loop
  for(i in 2:length(NNdist)){
    rnd <- sp::spsample(region, n=length(feature), type='random')
    matr <- rgeos::gDistance(rnd, byid=TRUE)
    diag(matr) <- NA
    NNdist[i] <- mean(apply(matr, 2, min, na.rm=TRUE))
    setTxtProgressBar(pb, i)
  }
  perm.p.value <- length(which(NNdist <= NNdist[1]))/B
  p.to.report <- ifelse(perm.p.value < 0.001, "< 0.001", ifelse(perm.p.value < 0.01, "< 0.01", ifelse(perm.p.value < 0.05, "< 0.05", round(perm.p.value, 3))))
  p.equalORlarger.to.report <- ifelse(1-perm.p.value < 0.001, "< 0.001", ifelse(1-perm.p.value < 0.01, "< 0.01", ifelse(1-perm.p.value < 0.05, "< 0.05", round(1-perm.p.value, 3))))
  if(addmap==TRUE){
    par(mfrow=c(1,2))
    plot(region, main="Map of the point dataset plus study area", cex.main=0.9, col=NA, border="red", lty=2)
    plot(feature, add=TRUE, pch=20, col="#00000088")
  } else {}
  dens <- density(NNdist)
  plot(dens, main=paste0("Nearest Neighbor Analysis: \ndensity of randomized average nearest neighbor distances across ", B, " iterations"), sub=paste0("Obs. aver. NN distance: ", round(NNdist[1],3), "; Random. aver. NN distance: ", round(mean(NNdist[-1]),3), "\nProb. of obs. aver. NN dist. <= random. aver. NN dist.: ",  p.to.report, " (Clustering [red tail])", "\nProb. of obs. aver. NN dist. >= random. aver. NN dist.: ", p.equalORlarger.to.report, " (Overdispersion [blue tail])"), xlab="", cex.main=0.9, cex.sub=0.70)
  polygon(dens, col = "#BCD2EE88", border = "blue")
  rug(NNdist, col = "#00000088")

  lower.lower <- min(NNdist[NNdist<=quantile(NNdist,0.05)])
  lower.upper <- max(NNdist[NNdist<=quantile(NNdist,0.05)])
  upper.lower <- min(NNdist[NNdist>=quantile(NNdist,0.95)])
  upper.upper <- max(NNdist[NNdist>=quantile(NNdist,0.95)])

  x1 <- min(which(dens$x >= lower.lower))
  x2 <- max(which(dens$x <=  lower.upper))
  x3 <- min(which(dens$x >= upper.lower))
  x4 <- max(which(dens$x <=  upper.upper))

  with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col="red"))
  with(dens, polygon(x=c(x[c(x3,x3:x4,x4)]), y= c(0, y[x3:x4], 0), col="blue"))

  abline(v=NNdist[1], lty = 2, col = "blue")
  points(x=mean(NNdist[-1]), y=0, pch=20, col="black")
  results <- list("obs.aver.NN.dist"=NNdist[1] , "rnd.aver.NN.dist" = mean(NNdist[-1]), "Prob. of obs. aver. NN dist. <= random. aver. NN dist."=round(perm.p.value,3), "Prob. of obs. aver. NN dist. >= random. aver. NN dist."=round(1-perm.p.value,3))
  return(results)
}
