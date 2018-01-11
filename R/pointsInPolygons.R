#' R function to test points-in-polygons association
#'
#' The function allows to test if there is a significant spatial association between a set of points and a set of polygons, in terms of points falling within the polygons.
#' In other words, it aims at testing whether a set of points falls inside a set of polygons more often than would be expected by chance.
#' The basic assumption is that the polygons are completely contained within the study plot.
#' The calculations are based on the bounding polygon based on the union the convex hulls of the point and of the polygon feature. The bounding polygon is considered as representing the study plot itself.
#'
#' The computational bases of the function are described in the help documentation of the 'Point-Polygon Relationship' analysis facility provided by the PASSaGE software (http://www.passagesoftware.net/manual.php).\cr
#'
#'The function makes use of the 'dbinom() and 'pbinom()' functions.\cr
#'The probability of observed count within polygons is 'dbinom(x, size=n.of.points, prob=p)', where 'x' is the observed number of points within polygons, 'n.of.points' is the total number of points, and 'p' is the probability that a single point will be found within a polygon, which is equal to the ratio between the area of the polygons and the total area of the study plot.\cr
#'The probability that x or fewer points will be found within the polygons is 'pbinom(x, size=n.of.points, prob=p)'.
#'
#' The function produces a plot of the points and polygons (plus the convex hull of the study area), and relevant information are reported at the bottom of the chart itself.\cr
#'
#' A list is also returned, containing what follows:\cr
#' -$Polygons' area;\cr
#' -$Study area's area;\cr
#' -$Total # of points;\cr
#' -$Observed # of points in polygons;\cr
#' -$Expected # of points in polygons;\cr
#' -$Exact probability of observed count within polygons;\cr
#' -$Probability of <= observed count within polygons;\cr
#' -$Probability of >= observed count within polygons.
#' @param point.feat: feature (of point type) whose spatial association with the polygons has to be assessed.
#' @param polyg.feat: feature (polygon type) in relation to which the spatial association of the points has to be assessed.
#' @param buffer: add a buffer to the convex hull of the study area (0 by default); the unit depends upon the units of the input data.
#' @keywords association
#' @export
#' @examples
#' data(points)
#' data(polygons)
#' result <- pointsInPolygons(points, polygons)
#'
pointsInPolygons <- function(point.feat, polyg.feat, buffer=0){
  options(scipen=999)
  region <- gConvexHull(raster::union(rgeos::gConvexHull(point.feat), rgeos::gConvexHull(polyg.feat)))  #'union()' requires 'raster'; build the convex hull of the union of the convex hulls of the two features
  ch <- rgeos::gBuffer(region, width=buffer)                                                            #add a buffer to the convex hull, with width is set by the 'buffer' parameter; the unit depends upon the units of the input data
  ch.area <- rgeos::gArea(ch)
  bypoly.area <- rgeos::gArea(polyg.feat, byid=TRUE)                                                    #area of each single polygon
  n.of.points <- length(point.feat)                                                                     #total number of points
  p <- sum(bypoly.area)/ch.area                                                                         #probability that a single point will be found within a polygon
  q <- 1-p                                                                                              #probability that a point will be found outside of the polygons
  exp.n.points <- round(n.of.points * p, 2)                                                             #for n total points, expected number of points within the polygons
  x <- sum(colSums(rgeos::gContains(polyg.feat, point.feat, byid = TRUE)))                              #observed number of points within polygons
  p.x <- round(dbinom(x, size=n.of.points, prob=p),5)                                                   #probability of observed count within polygons
  p.lessORequal <- round(pbinom(x, size=n.of.points, prob=p),5)                                         #cumulative probability from 1 to x
  p.largerORequal <- round(1-p.lessORequal,5)
  plot(polyg.feat, main="Points-in-Polygons relationship", cex.main=0.9, sub=paste0("Polygons' area: ", round(sum(bypoly.area),3),"; Study area's area: ", round(ch.area,3), "\nTotal number of points: ", n.of.points, "; Obs. number of points in polygons: ", x, "; Expect. number of points within polygons: ", exp.n.points, "\nExact prob. of obs. count within polygons: ", p.x, "; Prob. of <= obs. count within polygons: ", p.lessORequal, "; Prob. of >= obs. count within polygons: ", p.largerORequal ), cex.sub=0.80)
  plot(point.feat, pch=20, col="#00000088", add=TRUE)
  plot(ch, add=TRUE, col=NA, border="red", lty=2)
  results <- list("Polygons' area"=round(sum(bypoly.area),3), "Study area's area"=round(ch.area,3), "Total # of points"=n.of.points, "Observed # of points in polygons"=x , "Expected # of points in polygons" = exp.n.points, "Exact probability of observed count within polygons"=p.x, "Probability of <= observed count within polygons"=p.lessORequal, "Probability of >= observed count within polygons"=p.largerORequal)
  return(results)
}
