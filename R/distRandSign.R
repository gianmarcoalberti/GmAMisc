#' R function to test for a significant spatial association between two features (points-to-points, points-to-lines, points-to-polygons)
#'
#' The function allows to assess if there is a significant spatial association between a point pattern and the features of another pattern.
#' For instance, users may want to assess if some locations tend to lie close to some features represented by polylines.
#' By the same token, users may want to know if there is a spatial association between the location of a given event and the location of another event.
#' See the example provided further below (in the examples section), where the question to address is if there is a spatial association between springs and geological fault-lines;
#' in other words: do springs tend to be located near the geological faults?\cr
#'
#' Given a from-feature (event for which we want to estimate the spatial association with the to-feature) and a to-feature
#' (event in relation to which we want to estimate the spatial association for the from-feature), the assessment is performed by means of a randomized procedure:\cr
#'
#' -keeping fixed the location of the to-feature, random from-features are drawn B times (the number of randomized from-features is equal to the number of observed from-features);\cr
#' -for each draw, the average minimum distance to the to-features is calculated; if the to-feature is made up of polygons, the from-features falling within a polygon will have a distance of 0;\cr
#' -a distribution of average minimum distances is thus obtained;\cr
#' -p values are computed following Baddeley et al., "Spatial Point Patterns. Methodology and Applications with R", CRC Press 2016, p. 387.\cr
#'
#' The from-feature must be a point feature, whilst the to-feature can be a point or a polyline or a polygon feature.\cr
#'
#' The rationale of the procedure is that, if there indeed is a spatial association between the two features, the from-feature should be on average closer to the to-feature than randomly generated from-features.
#' If the studyplot shapefile is not provided, the random locations are drawn within a bounding polygon based on the union the convex hulls of the from- and of the to-feature.\cr
#'
#' If both the from-feature and the to-feature are of point type (SpatialPointsDataFrame class), the function also test the spatial association by means of
#' a permuted procedures. Unlike the procedure described above, whereby random points are drawn within the study area, the permutation-based routine
#' builds a distribution of averages minimum distances keeping the points location unchanged and randomly assigning the points to either of the two patterns.
#' The re-assigment is performed B times (199 by default) and each time the average minimum distance is calculated.\cr
#'
#' The function produces a histogram showing: the distribution of randomized average minimum distances;
#' a black dot indicating the observed average minimum distance;
#' a hollow dot representing the average of the randomized minimum distances;
#' two blue reference lines correspond to the 0.025 and to the 0.975 quantile of the randomized distribution.
#' P-values are reported at the bottom of the plot.\cr
#'
#' In case both the from- and the to- feature are of point type, another histogram is produced, which provides the same information of the preceding histogram, but
#' derived from the permutation-based routine that has been detailed above.
#'
#' A list is also returned, containing what follows:\cr
#' -$from.feat.min.dist: distance of each entity of the from-feature to the nearest entity of the to-feature;\cr
#' -$avrg.obs.min.dist: observed average minimum distance;\cr
#' -$avrg.rnd.min.dist: average randomized minimum distance;\cr
#' -$avrg.perm.min.dist: average permuted minimum distance (returned only when both the from- and to- features are of point type);\cr
#' -$p.value closer than expected-rnd-;\cr
#' -$p.value closer than expected-perm- (returned only when both the from- and to- features are of point type);\cr
#' -$p.value more distant than expected-rnd-;\cr
#' -$p.value more distant than expected-perm- (returned only when both the from- and to- features are of point type);\cr
#' -$p.value different from random-rnd-;\cr
#' -$p.value different from random-perm- (returned only when both the from- and to- features are of point type).\cr
#'
#' @param from.feat: feature (of point type; SpatialPointsDataFrame class) whose spatial association with the to-feature has to be assessed.
#' @param to.feat: feature (point, polyline, or polygon type; SpatialPointsDataFrame, SpatialLinesDataFrame, SpatialPolygonsDataFrame class) in relation to which the spatial association of the from-feature has to be assessed.
#' @param studyplot: feature (of polygon type; SpatialPolygonsDataFrame class) representing the study area; if not provided, the study area is internally worked out as the bounding polygon based on the union the convex hulls of the from- and of the to-feature.
#' @param buffer: add a buffer to the convex hull of the study area (0 by default); the unit depends upon the units of the input data.
#' @param B: number of randomizations to be used (199 by default).
#' @param oneplot: TRUE (default) or FALSE if the user wants or does not want the plots displayed in a single window.
#' @keywords distance
#' @export
#' @examples
#' Example 1:
#' data(springs)
#' data(faults)
#' result <- distRandSign(from.feat=springs, to.feat=faults, oneplot=TRUE) #calculate the significance of the spatial association between springs and geological fault-lines; plots displayed in a single panel
#'
#' Example 2:
#' data(malta_polyg)
#' result <- distRandSign(from.feat=springs, to.feat=faults, studyplot=malta_polyg, oneplot=FALSE) #same as above, but using a polygon (SpatialPolygonsDataFrame class) as studyplot; plots displayed separately
#'
#' Example 3:
#' data(points)
#' data(polygons)
#' result <- distRandSign(from.feat=points, to.feat=polygons, oneplot=FALSE) #calculate the significance of the spatial association between points and polygons
#' @seealso \code{\link{distRandCum}} , \code{\link{distCovarModel}} , \code{\link{Aindex}}
#'
distRandSign <- function(from.feat, to.feat, studyplot=NULL, buffer=0, B=199, oneplot=TRUE){

  options(scipen=999)

  if(is.null(studyplot)==TRUE){
    # build the convex hull of the union of the convex hulls of the two features
    ch <- gConvexHull(raster::union(gConvexHull(from.feat), gConvexHull(to.feat)))
    #add a buffer to the convex hull, with width is set by the 'buffer' parameter; the unit depends upon the units of the input data
    region <- gBuffer(ch, width=buffer)
  } else {
    region <- studyplot
  }

  #create an empty container for the observed average minimum distance
  obs.av.min.dist <- numeric()

  #create an empty container for the randomized average minimum distances
  rnd.av.min.dist <- numeric(B)

  #require 'rgeos' package; gDistance calculates all the pair-wise distances between each from-feature and to-feature
  res <- gDistance(from.feat, to.feat, byid=TRUE)

  #for each from-feature (i.e., column-wise), get the minimum distance to the to-feature
  obs.min.distances <- apply(res, 2, min)

  #average of the observed minimum distances
  obs.av.min.dist <- mean(obs.min.distances)

  #set the progress bar to be used inside the loop
  pb <- txtProgressBar(min = 0, max = B, style = 3)

  for(i in 1:B){
    #require 'sp' package; draw a random sample of points, with sample size equal to the number of from-feature points
    rnd <- spsample(region, n=length(from.feat), type='random')
    #calculate all the pair-wise distances from the from-feature to the to-feature
    res <- gDistance(rnd, to.feat, byid=TRUE)
    #calculate the average of the minimum distances and store them
    rnd.av.min.dist[i] <- mean(apply(res, 2, min))
    setTxtProgressBar(pb, i)
  }

  #calculate the p-value for a pattern featuring smaller than expected distances
  pclus <- (1 + sum (rnd.av.min.dist < obs.av.min.dist)) / (1 + B)
  pclus.to.report <- ifelse(pclus < 0.001, "< 0.001",
                            ifelse(pclus < 0.01, "< 0.01",
                                   ifelse(pclus < 0.05, "< 0.05",
                                          round(pclus, 3))))

  #calculate the p-value for a pattern featuring larger than expected distances
  preg <- (1 + sum (rnd.av.min.dist > obs.av.min.dist)) / (1 + B)
  preg.to.report <- ifelse(preg < 0.001, "< 0.001",
                           ifelse(preg < 0.01, "< 0.01",
                                  ifelse(preg < 0.05, "< 0.\05",
                                         round(preg, 3))))

  #calculate the 2-tailes p-value
  two.tailed.p <- 2 * (min(pclus, preg))
  two.tailed.to.report <- ifelse(two.tailed.p < 0.001, "< 0.001",
                                 ifelse(two.tailed.p < 0.01, "< 0.01",
                                        ifelse(two.tailed.p < 0.05, "< 0.05",
                                               round(two.tailed.p, 3))))


  #if both the from and to features are a point pattern, perform a randomized test
  #which, unlike the above section of the code, does not create random points within the study area but randomly assigns the
  #points membership to either the from or the to feature
  if((class(from.feat)[1]=="SpatialPointsDataFrame" | class(from.feat)[1]=="SpatialPoints") & (class(to.feat)[1]=="SpatialPointsDataFrame" | class(to.feat)[1]=="SpatialPoints")){
    #size of from.feat dataset
    nfrom <- length(from.feat)

    #size of to.feat dataset
    nto <- length(to.feat)

    #pool the coordinates of points in pattern from.feat and to.feat to be used later within the pointDistance function
    pooledData <- as.matrix(rbind(coordinates(from.feat), coordinates(to.feat)))

    #set the progress bar to used inside the loop
    pb <- txtProgressBar(min = 0, max = B, style = 3)

    #create an empty sloth for the permuted average minimum distances
    perm.av.min.dist <- numeric(length=B)

    for(i in 1:B){
      #create a random set of number to be used as index to extract the values from the pooledData matrix
      #the set of number has size equal to the number of points in pattern x
      index <- sample(1: (nfrom + nto), size=nfrom, replace=F)
      #extract from the pooledData matrix the rows corresponding to the generated random index
      from.perm <- pooledData[index,]
      #extract from the pooledData matrix all the rows that does not correspond to the random index
      #the procedure eventually creates two sets of points randomly shuffling the entries of the pooledData matrix
      to.perm <- pooledData[-index,]

      #calculate the average minimum distance using the permuted point patterns
      res <- pointDistance(from.perm, to.perm, lonlat=FALSE, allpairs=TRUE)
      #calculate the average of the minimum distances and store them
      #notice that the mean is applied row-wise (unlike to what done previously) because the distance matrix returned by pointDistance is
      #different from that returned by the gDistance function (see above)
      perm.av.min.dist[i] <- mean(apply(res, 1, min))

      setTxtProgressBar(pb, i)
    }

    #calculate the single tailed and the two-tailed permuted p values
    p.lower.perm <- (1 + sum (perm.av.min.dist < obs.av.min.dist)) / (1 + B)
    pclus.perm.to.report <- ifelse(p.lower.perm < 0.001, "< 0.001",
                                   ifelse(p.lower.perm < 0.01, "< 0.01",
                                          ifelse(p.lower.perm < 0.05, "< 0.05",
                                                 round(p.lower.perm, 3))))

    p.upper.perm <- (1 + sum (perm.av.min.dist > obs.av.min.dist)) / (1 + B)
    preg.perm.to.report <- ifelse(p.upper.perm < 0.001, "< 0.001",
                                  ifelse(p.upper.perm < 0.01, "< 0.01",
                                         ifelse(p.upper.perm < 0.05, "< 0.05",
                                                round(p.upper.perm, 3))))

    two.tailed.p.perm <- 2 * (min(p.lower.perm, p.upper.perm))
    p.two.t.perm.to.report <- ifelse(two.tailed.p.perm < 0.001, "< 0.001",
                                     ifelse(two.tailed.p.perm < 0.01, "< 0.01",
                                            ifelse(two.tailed.p.perm < 0.05, "< 0.05",
                                                   round(two.tailed.p.perm, 3))))
  } else {}

  #if 'oneplot' is TRUE
  if(oneplot==TRUE){

    #if both the input datasets are of point type
    if((class(from.feat)[1]=="SpatialPointsDataFrame" | class(from.feat)[1]=="SpatialPoints") & (class(to.feat)[1]=="SpatialPointsDataFrame" | class(to.feat)[1]=="SpatialPoints")){
      #set a layout featuring the map on top, and the two histograms side-by-side below the top panel
      m <- rbind(c(1,1), c(2,3))
      layout(m)

    } else {
      #if the two datasets are of different type (and oneplot is TRUE), set a layout featuting the map to the left and the histogram to the right
      par(mfrow=c(1,2))

    }
  }

  plot(region,
       main="Map of the from- and to-feature, plus study area",
       cex.main=0.9, col=NA,
       border="red",
       lty=2,
       axes=TRUE)
  plot(from.feat,
       add=TRUE,
       pch=20)
  plot(to.feat,
       add=TRUE)

  #plot the histogram of the randomized average minimum distances
  hist(rnd.av.min.dist,
       main=paste0("Feature-to-feature distance analysis: \nfreq. distribution of randomized average minimum distances across ", B, " iterations"),
       sub= paste0("Obs. aver. min. distance: ", round(obs.av.min.dist,3), "; Average of randomized min. distances: ", round(mean(rnd.av.min.dist),3), "\np-value closer than expected: ",  pclus.to.report, "; p-value more distant than expected: ", preg.to.report, "\np-value different from random: ", two.tailed.to.report),
       xlab="",
       cex.main=0.9,
       cex.sub=0.70)
  rug(rnd.av.min.dist, col = "#0000FF")
  abline(v=quantile(rnd.av.min.dist, 0.025), lty=2, col="blue")
  abline(v=quantile(rnd.av.min.dist, 0.975), lty=2, col="blue")
  points(x=mean(rnd.av.min.dist), y=0, pch=1, col="black")
  points(x=obs.av.min.dist, y=0, pch=20, col = "black")

  #plot an additional plot according to whether or not the input features are both a point pattern,
  #in which case the histogram shows the distribution of the permuted average minimum distances
  if((class(from.feat)[1]=="SpatialPointsDataFrame" | class(from.feat)[1]=="SpatialPoints") & (class(to.feat)[1]=="SpatialPointsDataFrame" | class(to.feat)[1]=="SpatialPoints")){

    hist(perm.av.min.dist,
         main=paste0("Feature-to-feature distance analysis: \nfreq. distribution of permuted average minimum distances across ", B, " iterations"),
         sub= paste0("Obs. aver. min. distance: ", round(obs.av.min.dist,3), "; Average of permuted min. distances: ", round(mean(perm.av.min.dist)), "\np-value closer than expected: ", pclus.perm.to.report,"; p-value more distant than expected: ", preg.perm.to.report,"\np-value different from random: ", p.two.t.perm.to.report),
         xlab="",
         cex.main=0.9,
         cex.sub=0.70)
    rug(perm.av.min.dist, col = "#0000FF")
    abline(v=quantile(perm.av.min.dist, 0.025), lty=2, col="blue")
    abline(v=quantile(perm.av.min.dist, 0.975), lty=2, col="blue")
    points(x=mean(perm.av.min.dist), y=0, pch=1, col="black")
    points(x=obs.av.min.dist, y=0, pch=20, col = "black")
  } else {}

  #if there are two inout point datasets, create a list also containing the permuted statistics,
  #otherwise create a list with only the randomized statistics
  if((class(from.feat)[1]=="SpatialPointsDataFrame" | class(from.feat)[1]=="SpatialPoints") & (class(to.feat)[1]=="SpatialPointsDataFrame" | class(to.feat)[1]=="SpatialPoints")){

    results <- list("from.feat.min.dist"=obs.min.distances,
                    "avrg.obs.min.dist" = obs.av.min.dist,
                    "avrg.rnd.min.dist"=mean(rnd.av.min.dist),
                    "avrg.perm.min.dist"=mean(perm.av.min.dist),
                    "p.value closer than expected-rnd-"=round(pclus,3),
                    "p.value closer than expected-perm-"=round(p.lower.perm,3),
                    "p.value more distant than expected-rnd-"=round(preg,3),
                    "p.value more distant than expected-perm-"=round(p.upper.perm,3),
                    "p.value different from random-rnd-"=round(two.tailed.p, 3),
                    "p.value different from random-perm-"=round(two.tailed.p.perm, 3))

  } else {

    results <- list("from.feat.min.dist"=obs.min.distances,
                    "avrg.obs.min.dist" = obs.av.min.dist,
                    "avrg.rnd.min.dist"=mean(rnd.av.min.dist),
                    "p.value closer than expected"=round(pclus,3),
                    "p.value more distant than expected"=round(preg,3),
                    "p.value different from random"=round(two.tailed.p, 3))
  }

  # restore the original graphical device's settings
  par(mfrow = c(1,1))

  return(results)
}
