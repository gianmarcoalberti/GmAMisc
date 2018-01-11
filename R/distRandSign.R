#' R function to calculate the significance of the spatial association between two features (points-to-points, points-to-lines, points-to-polygons)
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
#' -the significance (let's call it p) of the observed average minimum distance is calculated by counting how many randomized average minimum distances are equal or smaller than the observed one, and dividing the count by B. The probability that the observed average minimum distance is equal or larger than the randomized average minimum distance is equal to 1-p.\cr
#'
#' The from-feature must be a point feature, whilst the to-feature can be a point or a polyline or a polygon feature. The rationale of the procedure is that, if there indeed is a spatial association between the two features, the from-feature should be on average closer to the to-feature than randomly generated from-features.
#' The random locations are drawn within a bounding polygon based on the union the convex hulls of the from- and of the to-feature.\cr
#'
#' The function produces a plot representing the distribution of randomized average minimum distances, and a reference line indicating the observed average minimum distance. The p-values are reported at the bottom of the plot.\cr
#'
#' A list is also returned, containing what follows:\cr
#' -$from.feat.min.dist: distance of each entity of the from-feature to the nearest entity of the to-feature;\cr
#' -$avrg.obs.min.dist: observed average minimum distance;\cr
#' -$avrg.rnd.min.dist: randomized average minimum distance;\cr
#' -$Prob. of obs. aver. min. dist. <= random. aver. min. dist;\cr
#' -$Prob. of obs. aver. min. dist. >= random. aver. min. dist.
#' @param from.feat: feature (of point type) whose spatial association with the to-feature has to be assessed.
#' @param to.feat: feature (point, polyline, or polygon type) in relation to which the spatial association of the from-feature has to be assessed.
#' @param buffer: add a buffer to the convex hull of the study area (0 by default); the unit depends upon the units of the input data.
#' @param B: number of randomizations to be used (1000 by default).
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
#' data(points)
#' data(polygons)
#' result <- distRandSign(from.feat=points, to.feat=polygons, addmap=FALSE) #calculate the significance of the spatial association between points and polygons, but do not return the map.
#'
#' Example 3:
#' build two sets of spatially correlated points and then run the distRandSign() function
#' n <- 100
#' m <- 50
#' set.seed(17)
#' a <- matrix(rnorm(2*n), ncol=2)
#' b <- matrix(rnorm(2*m, a[1:m, ], sd=1/5), ncol=2)
#' a.loc <- SpatialPoints(a)
#' b.loc <- SpatialPoints(b)
#' result <- distRandSign(a.loc,b.loc, 500, addmap=TRUE) #calculate the points-to-points spatial association using 500 randomizations, and also returns the map
#'
distRandSign <- function(from.feat, to.feat, buffer=0, B=1000, addmap=TRUE){
  options(scipen=999)
  ch <- gConvexHull(raster::union(gConvexHull(from.feat), gConvexHull(to.feat)))  #'union()' requires 'raster'; build the convex hull of the union of the convex hulls of the two features
  region <- gBuffer(ch, width=buffer)                                             #add a buffer to the convex hull, with width is set by the 'buffer' parameter; the unit depends upon the units of the input data
  av.min.dist <- numeric(B+1)
  res <- gDistance(from.feat, to.feat, byid=TRUE)                                 #require 'rgeos' package; gDistance calculates all the pair-wise distances between each from-feature and to-feature
  obs.min.distances <- apply(res, 2, min)                                         #for each from-feature (i.e., column-wise), get the minimum distance to the to-feature
  av.min.dist[1] <- mean(obs.min.distances)                                       #average of the observed minimum distances
  for(i in 2:length(av.min.dist)){
    rnd <- spsample(region, n=length(from.feat), type='random')                   #require 'sp' package; draw a random sample of points, with sample size equal to the number of from-feature points
    res <- gDistance(rnd, to.feat, byid=TRUE)                                     #calculate all the pair-wise distances from the from-feature to the to-feature
    av.min.dist[i] <- mean(apply(res, 2, min))                                    #calculate the average of the minimum distances
  }
  perm.p.value <- length(which(av.min.dist <= av.min.dist[1]))/B                  #get the p-value considering how many randominzed average minimum distances are equal or smaller than the observed average minimum distance
  p.to.report <- ifelse(perm.p.value < 0.001, "< 0.001", ifelse(perm.p.value < 0.01, "< 0.01", ifelse(perm.p.value < 0.05, "< 0.05", round(perm.p.value, 3))))
  p.equalORlarger.to.report <- ifelse(1-perm.p.value < 0.001, "< 0.001", ifelse(1-perm.p.value < 0.01, "< 0.01", ifelse(1-perm.p.value < 0.05, "< 0.05", round(1-perm.p.value, 3))))
  if(addmap==TRUE){
    par(mfrow=c(1,2))
    plot(region, main="Map of the from- and to-feature, plus convex hull", cex.main=0.9, col=NA, border="red", lty=2)
    plot(from.feat, add=TRUE, pch=20)
    plot(to.feat, add=TRUE)
  } else {}
  d <- density(av.min.dist)
  plot(d, main=paste0("Feature-to-feature distance significance: \ndensity of randomized average minimum distances across ", B, " iterations"), sub=paste0("Obs. aver. min. distance: ", round(av.min.dist[1],3), "; Random. aver. min. distance: ", round(mean(av.min.dist[-1]),3), "\nProb. of obs. aver. min. dist. <= random. aver. min. dist: ",  p.to.report, "\nProb. of obs. aver. min. dist. >= random. aver. min. dist: ", p.equalORlarger.to.report), xlab="", cex.main=0.9, cex.sub=0.70)
  polygon(d, col = "#BCD2EE88", border = "blue")
  rug(av.min.dist, col = "#0000FF")
  abline(v=av.min.dist[1], lty = 2, col = "blue")
  results <- list("from.feat.min.dist"=obs.min.distances , "avrg.obs.min.dist" = av.min.dist[1], "avrg.rnd.min.dist"=mean(av.min.dist[-1]), "Prob. of obs. aver. min. dist. <= random. aver. min. dist"=round(perm.p.value,5), "Prob. of obs. aver. min. dist. >= random. aver. min. dist"=round(1-perm.p.value,5))
  return(results)
}
