#' R function to calculate the significance of the spatial relationship between two features (points-to-points, points-to-lines, points-to-polygons)
#'
#' The function allows to assess if there is a significant spatial association between two features.
#' For instance, users may want to assess if some locations tend to lie close to some features represented by polylines.
#' By the same token, users may want to know if there is a spatial association between the location of a given event and the location of another event. See the example provided further below (in the examples section), where the question to address is if there is a spatial association between springs and geological fault-lines; in other words: do springs tend to be located near the geological faults?\cr
#'
#' Given a from-feature (event for which we want to estimate the spatial association with the to-feature) and a to-feature (event in relation to which we want to estimate the spatial association for the from-feature), the assessment is performed by means of a randomized procedure:\cr
#'
#' -keeping fixed the location of the to-feature, random from-features are drawn B times (the number of randomized from-features is equal to the number of observed from-features);\cr
#' -for each draw, the average minimum distance to the to-features is calculated; if the to-feature is made up of polygons, the from-features falling within a polygon will have a distance of 0;\cr
#' -a distribution of average minimum distances is thus obtained;\cr
#' -p values are computed following Baddeley et al., "Spatial Point Patterns. Methodology and Applications with R", CRC Press 2016, p. 387.\cr
#'
#' The from-feature must be a point feature, whilst the to-feature can be a point or a polyline or a polygon feature. The rationale of the procedure is that, if there indeed is a spatial association between the two features, the from-feature should be on average closer to the to-feature than randomly generated from-features.
#' If the study plot shapefile is not provided, the random locations are drawn within a bounding polygon based on the union the convex hulls of the from- and of the to-feature.\cr
#'
#' The function produces a plot showing: the distribution of randomized average minimum distances; a black dot indicating the observed average minimum distance; a hollow dot representing the average of the randomized minimum distances; two blue reference lines correspond to the 0.025 and to the 0.975 percentile of the randomized distribution. P-values are reported at the bottom of the plot.\cr
#'
#' A list is also returned, containing what follows:\cr
#' -$from.feat.min.dist: distance of each entity of the from-feature to the nearest entity of the to-feature;\cr
#' -$avrg.obs.min.dist: observed average minimum distance;\cr
#' -$avrg.rnd.min.dist: randomized average minimum distance;\cr
#' -$Prob. of obs. aver. min. dist. < random. aver. min. dist;\cr
#' -$Prob. of obs. aver. min. dist. > random. aver. min. dist.
#' @param from.feat: feature (of point type; SpatialPointsDataFrame class) whose spatial association with the to-feature has to be assessed.
#' @param to.feat: feature (point, polyline, or polygon type; SpatialPointsDataFrame, SpatialLinesDataFrame, SpatialPolygonsDataFrame class) in relation to which the spatial association of the from-feature has to be assessed.
#' @param studyplot: feature (of polygon type; SpatialPolygonsDataFrame class) representing the study area; if not provided, the study area is internally worked out as the bounding polygon based on the union the convex hulls of the from- and of the to-feature.
#' @param buffer: add a buffer to the convex hull of the study area (0 by default); the unit depends upon the units of the input data.
#' @param B: number of randomizations to be used (199 by default).
#' @param addmap: TRUE (default) or FALSE if the user wants or does not want a map of the study area and of the from- and to-feature to be also displayed.
#' @keywords distance
#' @export
#' @examples
#' Example 1:
#' data(springs)
#' data(faults)
#' result <- distRandSign(from.feat=springs, to.feat=faults, addmap=TRUE) #calculate the significance of the spatial association between springs and geological fault-lines; it also returns a map.
#'
#' Example 2:
#' data(malta_polyg)
#' result <- distRandSign(from.feat=springs, to.feat=faults, studyplot=malta_polyg, addmap=TRUE) #same as above, but using a polygon (SpatialPolygonsDataFrame class) as studyplot
#'
#' Example 3:
#' data(points)
#' data(polygons)
#' result <- distRandSign(from.feat=points, to.feat=polygons, addmap=FALSE) #calculate the significance of the spatial association between points and polygons, but do not return the map.
#' @seealso \code{\link{distCovarModel}}
#'
distRandSign <- function(from.feat, to.feat, studyplot=NULL, buffer=0, B=199, addmap=TRUE){
  options(scipen=999)
  if(is.null(studyplot)==TRUE){
    ch <- gConvexHull(raster::union(gConvexHull(from.feat), gConvexHull(to.feat)))  #'union()' requires 'raster'; build the convex hull of the union of the convex hulls of the two features
    region <- gBuffer(ch, width=buffer)                                             #add a buffer to the convex hull, with width is set by the 'buffer' parameter; the unit depends upon the units of the input data
  } else {
    region <- studyplot
  }
  av.min.dist <- numeric(B+1)
  res <- gDistance(from.feat, to.feat, byid=TRUE)                                  #require 'rgeos' package; gDistance calculates all the pair-wise distances between each from-feature and to-feature
  obs.min.distances <- apply(res, 2, min)                                          #for each from-feature (i.e., column-wise), get the minimum distance to the to-feature
  av.min.dist[1] <- mean(obs.min.distances)                                        #average of the observed minimum distances, stored in the first slot of the object av.min.dist
  pb <- txtProgressBar(min = 0, max = B, style = 3)                                #set the progress bar to be used inside the loop
  for(i in 2:B){
    rnd <- spsample(region, n=length(from.feat), type='random')                    #require 'sp' package; draw a random sample of points, with sample size equal to the number of from-feature points
    res <- gDistance(rnd, to.feat, byid=TRUE)                                      #calculate all the pair-wise distances from the from-feature to the to-feature
    av.min.dist[i] <- mean(apply(res, 2, min))                                     #calculate the average of the minimum distances and store them in the av.min.dist object, starting from position 2
    setTxtProgressBar(pb, i)
  }
  pclus <- (1 + sum (av.min.dist[-1] < av.min.dist[1])) / (1 + B)                  #calculate the p-value for a clustered pattern
  pclus.to.report <- ifelse(pclus < 0.001, "< 0.001",
                            ifelse(pclus < 0.01, "< 0.01",
                                   ifelse(pclus < 0.05, "< 0.05",
                                          round(pclus, 3))))
  preg <- (1 + sum (av.min.dist[-1] > av.min.dist[1])) / (1 + B)                   #calculate the p-value for a regular pattern
  preg.to.report <- ifelse(preg < 0.001, "< 0.001",
                           ifelse(preg < 0.01, "< 0.01",
                                  ifelse(preg < 0.05, "< 0.05",
                                         round(preg, 3))))
  if(addmap==TRUE){
    par(mfrow=c(1,2))
    plot(region,
         main="Map of the from- and to-feature, plus convex hull",
         cex.main=0.9, col=NA,
         border="red",
         lty=2)
    plot(from.feat,
         add=TRUE,
         pch=20)
    plot(to.feat,
         add=TRUE)
  } else {}
  d <- density(av.min.dist)
  plot(d,
       main=paste0("Feature-to-feature distance significance: \ndensity of randomized average minimum distances across ", B, " iterations"),
       sub=paste0("Obs. aver. min. distance: ", round(av.min.dist[1],3), "; Random. aver. min. distance: ", round(mean(av.min.dist[-1]),3), "\nProb. of obs. aver. min. dist. < random. aver. min. dist.: ",  pclus.to.report, "\nProb. of obs. aver. min. dist. > random. aver. min. dist.: ", preg.to.report), xlab="",
       cex.main=0.9,
       cex.sub=0.70)
  polygon(d,
          col = "#BCD2EE88",
          border = "blue")
  rug(av.min.dist, col = "#0000FF")
  abline(v=quantile(av.min.dist, 0.025), lty=2, col="blue")
  abline(v=quantile(av.min.dist, 0.975), lty=2, col="blue")
  points(x=mean(av.min.dist[-1]), y=0, pch=1, col="black")
  points(x=av.min.dist[1], y=0, pch=20, col = "black")
  results <- list("from.feat.min.dist"=obs.min.distances ,
                  "avrg.obs.min.dist" = av.min.dist[1],
                  "avrg.rnd.min.dist"=mean(av.min.dist[-1]),
                  "Prob. of obs. aver. min. dist. < random. aver. min. dist"=round(pclus,5),
                  "Prob. of obs. aver. min. dist. > random. aver. min. dist"=round(preg,5))
  return(results)
}
