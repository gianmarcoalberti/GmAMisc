#' R function for Nearest Neighbor analysis of point patterns
#'
#' The function allows to perform the Nearest Neighbor analysis of point patterns to formally test for the presence of a clustered, dispersed, or random spatial arrangement (second-order effect).
#' It also allows to control for a first-order effect (i.e., influence of an underlaying numerical covariate) while performing the analysis.
#' The covariate must be of RasterLayer class. Significance is assessed via a randomized approach.\cr
#'
#' The function uses a randomized approach to test the significance of the Clark-Evans R statistic:
#' the observed R value is set against the distribution of R values computed across B iterations (199 by default) in which a set of random points
#' (with a sample size equal to the number of points of the input feature) is drawn and the statistic recomputed.\cr
#'
#' The function produces a histogram of the randomized R values, with a black dot indicating the observed value and a hollow dot representing the average of the randomized R values.
#' P-values (computed following Baddeley et al., "Spatial Point Patterns. Methodology and Applications with R", CRC Press 2016, p. 387), are reported at the bottom of the same chart.
#' Two reference lines represent the two tails of the randomized distribution (left tail, indicating a significant clustered pattern; right tail, indicating a significant dispersed pattern).\cr
#'
#' The function also returns a list storing the following:\cr
#' -$obs.NN.dist: observed NN distances;\cr
#' -$obs.R: observed R value;\cr
#' -$aver.rand.R: average randomized R;\cr
#' -$p.value clustered: p-value for a clustered pattern;\cr
#' -$p.value.dispersed: p-value for a dispersed pattern;\cr
#' -$p.value.diff.from.random: p-value for a pattern different from random.\cr
#'
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
#' res <- NNa(rndpoints, buffer=100, B=999) #perform the analysis using the 'rndpoints' dataset, add a 100m buffer arounf the points' convexhull, and use 499 iterations; the result points to a random arrangement
#'
#' data(Starbucks)
#' data(popdensity)
#' res <- NNa(Starbucks, cov.var=popdensity)  #perform the analysis, while controlling for the effect of the population density covariate
#' @seealso \code{\link{refNNa}}
#'
NNa <- function(feature, studyplot=NULL, buffer=0, B=199, cov.var=NULL, addmap=TRUE){
  #define an ojbect in which we can store the observed average NN distance
  obs.NNdist <- numeric()

  #do the same for the R statistics
  obs.Rindex <- numeric()

  #define an ojbect in which we can store the randomized average NN distances
  rnd.NNdist <- numeric(B)

  #do the same for the R statistics
  rnd.Rindex <- numeric(B)

  #calculate all the pair-wise distances among points
  matr <- rgeos::gDistance(feature, byid=TRUE)

  #set to NA the values along the diagonal of the matrix,
  #which represent the distance between anypoint and itself
  diag(matr) <- NA

  #calculate the observed average NN distance, and store it
  obs.NNdist <- mean(apply(matr, 2, min, na.rm=TRUE))

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

    #calculate the density of points
    pointdensity <- length(feature) / area(region)

    #calculate the observed R value, and store it
    obs.Rindex <- obs.NNdist / (1/(2*sqrt(pointdensity)))

    #create a loop that draws random points within the region and calculate the randomized average NN distances storing
    #in the rnd.NNdist object
    for(i in 1:B){
      rnd <- sp::spsample(region, n=length(feature), type='random')
      matr <- rgeos::gDistance(rnd, byid=TRUE)
      diag(matr) <- NA
      rnd.NNdist[i] <- mean(apply(matr, 2, min, na.rm=TRUE))
      rnd.Rindex[i] <- rnd.NNdist[i] / (1/(2*sqrt(pointdensity)))
      setTxtProgressBar(pb, i)
    }

    #(if the covariate raster is provided; see above), then...
  } else {

    #tranform the cov.var from a RasterLayer to an object of class im, which is needed by spatstat
    cov.var.im <- as.im(cov.var)

    #calculate the density of points
    pointdensity <- length(feature) / area(cov.var.im)

    #calculate the observed R value, and store it
    obs.Rindex <- obs.NNdist / (1/(2*sqrt(pointdensity)))

    #draw random points via the spatstat's rpoint function,
    #using the covariate dataset as spatial covariate
    for (i in 1:B){
      rnd   <- spatstat::rpoint(n=length(feature), f=cov.var.im)
      rnd.NNdist[i] <- mean(nndist(rnd, k=1))
      rnd.Rindex[i] <- rnd.NNdist[i] / (1/(2*sqrt(pointdensity)))
      setTxtProgressBar(pb, i)
    }
  }

  #calculate the p-value for a clustered pattern
  pclus <- (1 + sum (rnd.Rindex < obs.Rindex[1])) / (1 + B)
  pclus.to.report <- ifelse(pclus < 0.001, "< 0.001",
                            ifelse(pclus < 0.01, "< 0.01",
                                   ifelse(pclus < 0.05, "< 0.05",
                                          round(pclus, 3))))

  #calculate the p-value for a regular pattern
  preg <- (1 + sum (rnd.Rindex > obs.Rindex[1])) / (1 + B)
  preg.to.report <- ifelse(preg < 0.001, "< 0.001",
                           ifelse(preg < 0.01, "< 0.01",
                                  ifelse(preg < 0.05, "< 0.05",
                                         round(preg, 3))))

  #calculate the 2-tailed p-value
  two.tailed.p <- 2 * (min(pclus, preg))
  two.tail.to.report <- ifelse(two.tailed.p < 0.001, "< 0.001",
                               ifelse(two.tailed.p < 0.01, "< 0.01",
                                      ifelse(two.tailed.p < 0.05, "< 0.05",
                                             round(two.tailed.p, 3))))

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

  #adjust the main plot title according to the presence of a covariate dataset
  if(is.null(cov.var)==TRUE){
    maintitle <- paste0("Nearest Neighbor Analysis: \nfreq. distrib. of randomized R values across ", B, " iterations")
  } else {
    maintitle <- paste0("Nearest Neighbor Analysis (controlling for the covariate effect): \nfreq. distrib. of randomized R values across ", B, " iterations")
  }

  hist(rnd.Rindex,
       main=maintitle,
       sub=paste0("Observed R: ", round(obs.Rindex,3), "; Average of randomized R: ", round(mean(rnd.Rindex),3), "\np-value clustered: ",  pclus.to.report, "; p-value dispersed: ", preg.to.report, "\np-value different from random: ", two.tail.to.report, "\nR<1 (clustered); R=1 (random); R>1 (dispersed)"),
       xlab="",
       cex.main=0.90,
       cex.sub=0.70,
       cex.axis=0.85)
  rug(rnd.Rindex, col = "#0000FF")
  abline(v=quantile(rnd.Rindex, 0.025), lty=2, col="blue")
  abline(v=quantile(rnd.Rindex, 0.975), lty=2, col="blue")
  points(x=mean(rnd.Rindex), y=0, pch=1, col="black")
  points(x=obs.Rindex, y=0, pch=20, col = "black")

  results <- list("obs.NN.dist"=obs.NNdist,
                  "aver.rand.NN.dist" = mean(rnd.NNdist),
                  "obs.R"=obs.Rindex,
                  "aver.rand.R"= mean(rnd.Rindex),
                  "p.value clustered"=round(pclus,3),
                  "p.value.dispersed"=round(preg,3),
                  "p.value.diff.from.random"=round(two.tailed.p,3))

  # restore the original graphical device's settings
  par(mfrow = c(1,1))

  return(results)
}
