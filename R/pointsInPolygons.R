#' R function to test points-in-polygons relationship
#'
#' The function allows to test:\cr
#' -scenario a: if there is a significant spatial association between a set of points and a set of polygons, in terms of points falling within the polygons. In other words, it aims at testing whether a set of points falls inside a set of polygons more often than would be expected by chance. The basic assumption is that the polygons are completely contained within the study plot.
#' If the shapefile (of polygon type) representing the study plot is not provided, the calculations use the bounding polygon based on the union the convex hulls of the point and of the polygon feature.\cr
#' -scenario b: if the distribution of points within a set of polygons totally covering the study area can be considered random, or if the observed points count for each polygon is larger or smaller than expected. P values are also reported.\cr
#'
#' The computations relative to scenario "a" are based on the 'dbinom() and 'pbinom()' functions. The probability of observed count within polygons is 'dbinom(x, size=n.of.points, prob=p)', where 'x' is the observed number of points within polygons, 'n.of.points' is the total number of points, and 'p' is the probability that a single point will be found within a polygon, which is equal to the ratio between the area of the polygons and the total area of the study plot. The probability that x or fewer points will be found within the polygons is 'pbinom(x, size=n.of.points, prob=p)'.
#'
#' The calculations relative to the scenario "b" are again based on the binomial distribution: the probability of the observed counts is 'dbinom(x, size=n.of.points, prob=p)', where 'x' is the observed number of points within a given polygon, 'n.of.points' is the total number of points, and 'p' is equal to the size of each polygon relative to sum of the polygons' area. The probability that x or fewer points will be found within a given polygon is 'pbinom(x, size=n.of.points, prob=p)'.\cr
#'
#' For scenario "a" the function produces a plot of the points and polygons (plus  the study area), and relevant information are reported at the bottom of the chart itself.\cr
#' A list is also returned, containing what follows:\cr
#' -$Polygons' area;\cr
#' -$Study area's area;\cr
#' -$Total # of points;\cr
#' -$Observed # of points in polygons;\cr
#' -$Expected # of points in polygons;\cr
#' -$Exact probability of observed count within polygons;\cr
#' -$Probability of <= observed count within polygons;\cr
#' -$Probability of >= observed count within polygons.
#'
#' For scenario "b" the function returns a plot showing the polygons plus the dots; in each polygon the observed and expected counts are reported, and the p-value of the observed count is indicated.\cr
#' A matrix is also returned, containing what follows:\cr
#' -polygons' area;\cr
#' -percentage area (size of each polygon relative to sum of the polygons' area; it corresponds to the probability (p) fed into the binomial distribution function);\cr
#' -observed number of points;\cr
#' -expected number of points;\cr
#' -probability of observed counts;\cr
#' -probability of observed counts <= than expected;\cr
#' -probability of observed counts >= than expected.\cr
#' @param point.feat: feature (of point type) whose spatial association with the polygons has to be assessed.
#' @param polyg.feat: feature (polygon type) in relation to which the spatial association of the points has to be assessed.
#' @param studyplot: shapefile (of polygon type) representing the study area; if not provided, the study area is internally worked out as the bounding polygon based on the union the convex hulls of the point and of the polygon feature.
#' @param scenario: select one of the two types of analysis available ("a" or "b").
#' @param buffer: add a buffer to the convex hull of the study area (0 by default); the unit depends upon the units of the input data.
#' @param cex.text: modify the size of the labels in the plot produced by the 'scenario b' option.
#' @keywords association
#' @export
#' @examples
#' Example 1
#' data(points)
#' data(polygons)
#' result <- pointsInPolygons(points, polygons, scenario="a")
#'
#'Example 2
#' data(events)
#' data(thiessenpolyg)
#' result <- pointsInPolygons(events, thiessenpolyg, scenario="b")
#'
pointsInPolygons <- function(point.feat, polyg.feat, studyplot=NULL, scenario, buffer=0, cex.text=0.7){
  options(scipen=999)
  if(scenario=="a"){
    if(is.null(studyplot)==TRUE){
      region <- gConvexHull(raster::union(rgeos::gConvexHull(point.feat), rgeos::gConvexHull(polyg.feat)))  #'union()' requires 'raster'; build the convex hull of the union of the convex hulls of the two features
    } else {
      region <- studyplot
    }
    ch <- rgeos::gBuffer(region, width=buffer)                                                              #add a buffer to the convex hull, with width is set by the 'buffer' parameter; the unit depends upon the units of the input data
    ch.area <- rgeos::gArea(ch)
    bypoly.area <- rgeos::gArea(polyg.feat, byid=TRUE)                                                      #area of each single polygon
    n.of.points <- length(point.feat)                                                                       #total number of points
    p <- sum(bypoly.area)/ch.area                                                                           #probability that a single point will be found within a polygon
    q <- 1-p                                                                                                #probability that a point will be found outside of the polygons
    exp.n.points <- round(n.of.points * p, 2)                                                               #for n total points, expected number of points within the polygons
    x <- sum(colSums(rgeos::gContains(polyg.feat, point.feat, byid = TRUE)))                                #observed number of points within polygons
    p.x <- round(dbinom(x, size=n.of.points, prob=p),5)                                                     #probability of observed count within polygons
    p.lessORequal <- round(pbinom(x, size=n.of.points, prob=p),5)                                           #cumulative probability from 1 to x
    p.largerORequal <- round(1-p.lessORequal,5)
    plot(ch, col=NA, border="red", lty=2, main="Points-in-Polygons relationship", cex.main=0.9, sub=paste0("Polygons' area: ", round(sum(bypoly.area),3),"; Study area's area: ", round(ch.area,3), "\nTotal number of points: ", n.of.points, "; Obs. number of points in polygons: ", x, "; Expect. number of points within polygons: ", exp.n.points, "\nExact prob. of obs. count within polygons: ", p.x, "; Prob. of <= obs. count within polygons: ", p.lessORequal, "; Prob. of >= obs. count within polygons: ", p.largerORequal ), cex.sub=0.80)
    plot(polyg.feat, add=TRUE)
    plot(point.feat, pch=20, col="#00000088", add=TRUE)
    results <- list("Polygons' area"=round(sum(bypoly.area),3), "Study area's area"=round(ch.area,3), "Total # of points"=n.of.points, "Observed # of points in polygons"=x , "Expected # of points in polygons" = exp.n.points, "Exact probability of observed count within polygons"=p.x, "Probability of <= observed count within polygons"=p.lessORequal, "Probability of >= observed count within polygons"=p.largerORequal)
    return(results) } else {
      bypoly.area <- rgeos::gArea(polyg.feat, byid = TRUE)
      perc.area <- bypoly.area / sum(bypoly.area)
      tot.n.points <- length(point.feat)
      bypoly.points <- colSums(rgeos::gContains(polyg.feat, point.feat, byid = TRUE))
      exp.n.points <- tot.n.points * perc.area
      p.obs <- dbinom(bypoly.points, tot.n.points, perc.area)
      pLessOrEqual <- pbinom(bypoly.points, tot.n.points, perc.area)
      pLargerOrEqual <- 1 - pLessOrEqual
      df <- matrix(nrow=length(polyg.feat), ncol=7)
      colnames(df)<- c("polygon.area", "%area", "obs.n.points", "exp.n.points", "p.obs", "p.<=exp", "p.>=exp")
      df[,1] <- bypoly.area
      df[,2] <- round(perc.area,2)
      df[,3] <- bypoly.points
      df[,4] <- round(exp.n.points,2)
      df[,5] <- round(p.obs,5)
      df[,6] <- round(pLessOrEqual,5)
      df[,7] <- round(pLargerOrEqual,5)
      plot(polyg.feat, main="Points-Polygons relationship", sub="By polygon observed/expected counts, and p-value of obs. counts", cex.main=0.9, cex.sub=0.8)
      plot(point.feat, pch=20, col="#ff000088", add=TRUE)
      plot(polyg.feat, add=TRUE)
      p.to.report <- ifelse(p.obs < 0.001, "< 0.001", ifelse(p.obs < 0.01, "< 0.01", ifelse(p.obs < 0.05, "< 0.05", round(p.obs, 3))))
      labls <- paste0(p.to.report, "\n(obs: ", df[,3],"; exp: ", round(df[,4],2), ")")
      text(x=coordinates(polyg.feat)[,1], y=coordinates(polyg.feat)[,2], labels=labls, cex=cex.text)
      return(df)
    }
}
