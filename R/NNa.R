#' R function for Nearest Neighbor analysis of point patterns
#'
#' The function allows to perform the Nearest Neighbor analysis of point patterns to formally test for the presence of a clustered, dispersed, or random spatial arrangement (second-order effect).
#' It also allows to controll for a first-order effect (i.e., influence of an underlaying numerical covariate) while performing the analysis. The covariate must be of RasterLayer class.
#' Significance is assessed via a randomized approach.\cr
#'
#' The function uses a randomized approach to test the significance of the Nearest Neighbor distance:
#' the observed average NN distance is compared against the distribution of average NN distances computed across B iterations.
#' In each iteration, a set of random points (with a sample size equal to the number of points of the input feature) is drawn.\cr
#'
#' The function produces a density chart of the randomized average NN distances, with a black dot indicating the observed average NN and a hollow dot representing the average of the randomized NN distances.
#' P-values (computed following Baddeley et al., "Spatial Point Patterns. Methodology and Applications with R", CRC Press 2016, p. 387) are reported at the bottom of the same chart.
#' Two reference lines represent the two tails of the randomized distribution (left tail, indicating a significant clustered pattern; right tail, indicating a significant dispersed pattern).\cr
#'
#' The function also returns a list storing the following:\cr
#' -$obs.aver.NN.dist\cr
#' -$rnd.aver.NN.dist\cr
#' -$Prob. of obs. aver. NN dist. < random. aver. NN dist.\cr
#' -$Prob. of obs. aver. NN dist. > random. aver. NN dist.\cr
#' @param feature: feature dataset (of point type; SpatialPointsDataFrame class).
#' @param studyplot: shapefile (of polygon type; SpatialPolygonsDataFrame class) representing the study area; if not provided, the study area is internally worked out as the convex hull enclosing the input feature dataset.
#' @param buffer: add a buffer to the studyplot (0 by default); the unit depends upon the units of the input data.
#' @param cov.var: numeric covariate (of RasterLayer class).
#' @param B: number of randomizations to be used (199 by default).
#' @param addmap: TRUE (default) or FALSE if the user wants or does not want a map of the study area and of feature dataset to be also displayed.
#' @keywords NNA
#' @export
#' @examples
#' data(springs)
#' res <- NNa(springs) #perform the analysis using all default values; the result points to a significant clustering
#'
#' data(springs)
#' data(malta_polyg)
#' res <- NNa(springs, studyplot=malta_polyg) #same as above but using a polygon (SpatialPolygonsDataFrame) as studyplot
#'
#' data(rndpoints)
#' res <- NNa(rndpoints, buffer=100, B=499) #perform the analysis using the 'rndpoints' dataset, add a 100m buffer arounf the points' convexhull, and use 499 iterations; the result points to a random arrangement
#'
#' data(Starbucks)
#' data(popdensity)
#' res <- NNa(Starbucks, cov.var=popdensity)  #perform the analysis, while controlling for the effect of the population density covariate
#' @seealso \code{\link{refNNa}}
#'
NNa <- function(feature, studyplot=NULL, buffer=0, B=199, cov.var=NULL, addmap=TRUE){
  #define an ojbect in which we can store the observed average NN distance (first slot)
  #and the randomized average NN distance (from slot 2 on)
  NNdist <- numeric(B+1)

  #calculate all the pair-wise distances among points
  matr <- rgeos::gDistance(feature, byid=TRUE)

  #set to NA the values along the diagonal of the matrix,
  #which represent the distance between anypoint and itself
  diag(matr) <- NA

  #calculate the observed average NN distance, and store it in the first slot of the NNdist object
  NNdist[1] <- mean(apply(matr, 2, min, na.rm=TRUE))

  #set the progress bar to be used later inside the loops
  pb <- txtProgressBar(min = 0, max = B, style = 3)

  #if there is no covariate data, workout the studyplot according to whether or not the studyplot is
  #entered by the user; if it is not, the studyplot is the convex hull based on the points themselves;
  #either way, the studyplot is eventually stored into the region object
  if(is.null(cov.var)==TRUE){

    if(is.null(studyplot)==TRUE){
      ch <- rgeos::gConvexHull(feature)
      region <- rgeos::gBuffer(ch, width=buffer)
    } else {
      region <- studyplot
    }

    #create a loop that draws random points within the region and calculate the randomized average NN distances storing
    #in the NNdist object from the second slot on
    for(i in 2:B){
      rnd <- sp::spsample(region, n=length(feature), type='random')
      matr <- rgeos::gDistance(rnd, byid=TRUE)
      diag(matr) <- NA
      NNdist[i] <- mean(apply(matr, 2, min, na.rm=TRUE))
      setTxtProgressBar(pb, i)
    }

    #(if the covariate raster is provided; see above), then...
  } else {

    #tranform the cov.var from a RasterLayer to an object of class im, which is needed by spatstat
    cov.var.im <- as.im(cov.var)

    #draw random points via the spatstat's rpoint function,
    #using the covariate dataset as spatial covariate
    for (i in 2:B){
      rnd   <- spatstat::rpoint(n=length(feature), f=cov.var.im)
      NNdist[i] <- mean(nndist(rnd, k=1))
      setTxtProgressBar(pb, i)
    }
  }

  #calculate the p-value for a clustered pattern
  pclus <- (1 + sum (NNdist[-1] < NNdist[1])) / (1 + B)
  pclus.to.report <- ifelse(pclus < 0.001, "< 0.001",
                            ifelse(pclus < 0.01, "< 0.01",
                                   ifelse(pclus < 0.05, "< 0.05",
                                          round(pclus, 3))))

  #calculate the p-value for a regular pattern
  preg <- (1 + sum (NNdist[-1] > NNdist[1])) / (1 + B)
  preg.to.report <- ifelse(preg < 0.001, "< 0.001",
                           ifelse(preg < 0.01, "< 0.01",
                                  ifelse(preg < 0.05, "< 0.05",
                                         round(preg, 3))))

  if(addmap==TRUE){
    par(mfrow=c(1,2))
    if(is.null(cov.var)==FALSE){
      plot(cov.var,
           main="Map of the point dataset against covariate",
           cex.main=0.9)
    } else {
      plot(region,
           main="Map of the point dataset plus study area",
           cex.main=0.9,
           col=NA,
           border="red",
           lty=2)}
    plot(feature,
         add=TRUE,
         pch=20,
         col="#00000088")
  } else {}

  dens <- density(NNdist)

  #adjust the main plot title according to the presence of a covariate dataset
  if(is.null(cov.var)==TRUE){
    maintitle <- paste0("Nearest Neighbor Analysis: \ndensity of randomized average nearest neighbor distances across ", B, " iterations")
  } else {
    maintitle <- paste0("Nearest Neighbor Analysis (controlling for the covariate effect): \ndensity of randomized average nearest neighbor distances across ", B, " iterations")
  }

  plot(dens,
       main=maintitle,
       sub=paste0("Obs. aver. NN distance: ", round(NNdist[1],3), "; Random. aver. NN distance: ", round(mean(NNdist[-1]),3), "\nProb. of obs. aver. NN dist. < random. aver. NN dist.: ",  pclus.to.report, " (Clustered pattern [left tail])", "\nProb. of obs. aver. NN dist. > random. aver. NN dist.: ", preg.to.report, " (Dispersed pattern [right tail])"),
       xlab="",
       cex.main=0.9,
       cex.sub=0.70)
  polygon(dens, col = "#BCD2EE88", border = "blue")
  rug(NNdist, col = "#0000FF")
  abline(v=quantile(NNdist, 0.025), lty=2, col="blue")
  abline(v=quantile(NNdist, 0.975), lty=2, col="blue")
  points(x=mean(NNdist[-1]), y=0, pch=1, col="black")
  points(x=NNdist[1], y=0, pch=20, col = "black")

  results <- list("obs.aver.NN.dist"=NNdist[1] ,
                  "rnd.aver.NN.dist" = mean(NNdist[-1]),
                  "Prob. of obs. aver. NN dist. < random. aver. NN dist."=round(pclus,3),
                  "Prob. of obs. aver. NN dist. > random. aver. NN dist."=round(preg,3))

  return(results)
}
