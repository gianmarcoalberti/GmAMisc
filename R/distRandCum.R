#' R function to test the significance of the spatial relationship between two features in terms of the cumulative distribution of minimum distances
#'
#' The function allows to assess if there is a significant spatial association between a point pattern and the features of another pattern.
#' For instance, users may want to assess if the features of a point pattern tend to lie close to some features represented by polylines.\cr
#'
#' Given a from-feature (event for which we want to estimate the spatial association with the to-feature) and a to-feature
#' (event in relation to which we want to estimate the spatial association for the from-feature), the assessment is performed by means of q randomized procedure:\cr
#'
#' -keeping fixed the location of the to-feature, random from-features are drawn B times (the number of randomized from-features is equal to the number of observed from-features);\cr
#' -for each draw, the minimum distance to the to-features is calculated; if the to-feature is made up of polygons, the from-features falling within a polygon will have a distance of 0;\cr
#' -a cumulative distribution of random minimum distances is thus obtained;\cr
#' -the cumulative random minimum distances are used to work out a 95percent confidence envelope that allows to assess the statistical significance of the
#' cumulative distribution of the observed minimum distances.\cr
#'
#' The from-feature must be a point feature, whilst the to-feature can be a point or a polyline or a polygon feature.\cr
#'
#' The rationale of the procedure is that, if there indeed is a spatial association between the two features, the from-feature should be closer to the to-feature than randomly generated from-features.
#' If the studyplot shapefile is not provided, the random locations are drawn within a bounding polygon based on the union the convex hulls of the from- and of the to-feature.\cr
#'
#' If both the from-feature and the to-feature are of point type (SpatialPointsDataFrame class), the user may opt for the randomized procedure described above (parameter 'type' set to 'rand'),
#' or for a permutation-based procedure (parameter 'type' set to 'perm'). Unlike the procedure described above, whereby random points are drawn within the study area, the permutation-based routine
#' builds a cumulative distribution of minimum distances keeping the points location unchanged and randomly assigning the points to either of the two patterns.
#' The re-assigment is performed B times (200 by default) and each time the minimum distance is calculated.\cr
#'
#' The function produces a cumulative distribution chart in which the distribution of the observed minimum distances is represented by a black line,
#' and the 95percent confidence envelope is represented in grey. The number of iteration used and the type of analysis
#' (whether randomization-based or permutation-based) are reported in the chart's title.\cr
#'
#' @param from.feat: feature (of point type; SpatialPointsDataFrame class) whose spatial association with the to-feature has to be assessed.
#' @param to.feat: feature (point, polyline, or polygon type; SpatialPointsDataFrame, SpatialLinesDataFrame, SpatialPolygonsDataFrame class) in relation to which the spatial association of the from-feature has to be assessed.
#' @param studyplot: feature (of polygon type; SpatialPolygonsDataFrame class) representing the study area; if not provided, the study area is internally worked out as the bounding polygon based on the union the convex hulls of the from- and of the to-feature.
#' @param buffer: add a buffer to the convex hull of the study area (0 by default); the unit depends upon the units of the input data.
#' @param B: number of randomizations to be used (200 by default).
#' @param type: by default is set to "rand", which performs the randomization-based analysis; if both the from.feature and the to.feature dataset are of
#' point type, setting the parameter to "perm" allows to opt for the permutation-based approach.
#' @keywords distance
#' @export
#' @examples
#' data(springs)
#' data(faults)
#' distRandCum(from.feat=springs, to.feat=faults) #perform the analysis using the default 200 iterations and the default randomization-based approach
#'
#' data(malta_polyg)
#' distRandCum(from.feat=springs, to.feat=faults, studyplot=malta_polyg, B=100) #same as above, but using a polygon (SpatialPolygonsDataFrame class) as studyplot and using 100 iteration
#'
#' data("malta_polyg") # load a sample polygon
#' pA <- spsample(malta_polyg, n=30, type='random')  #create a set of 30 random points within the polygon
#' pB <- spsample(malta_polyg, n=40, type='random')  #create a set of 40 random points within the polygon
#' distRandCum(pA, pB, studyplot=malta_polyg) #perform the analysis; since both patterns are of point type but the 'type' parameter is left in its default value ('rand'), the
#' randomization-based approach is used
#'
#' distRandCum(pA, pB, studyplot=malta_polyg, type="perm") #same as above, but using the permutation-based approach
#'
#' @seealso \code{\link{distRandSign}}
#'
distRandCum <- function (from.feat, to.feat, studyplot=NULL, buffer=0, B=200, type="rand") {

  if(is.null(studyplot)==TRUE){
    # build the convex hull of the union of the convex hulls of the two features
    ch <- gConvexHull(raster::union(gConvexHull(from.feat), gConvexHull(to.feat)))
    #add a buffer to the convex hull, with width is set by the 'buffer' parameter; the unit depends upon the units of the input data
    region <- gBuffer(ch, width=buffer)
  } else {
    region <- studyplot
  }

  #require 'rgeos' package; gDistance calculates all the pair-wise distances between each from-feature and to-feature
  dst <- gDistance(from.feat, to.feat, byid=TRUE)

  #for each from-feature (i.e., column-wise), get the minimum distance to the to-feature
  obs.min.dst <- apply(dst, 2, min)

  #calculate the ECDF of the observed minimum distances
  dst.ecdf <- ecdf(obs.min.dst)

  #create a matrix to store the distance of each random point to its closest to.feature;
  #each column correspond to a random set of points
  dist.rnd.mtrx <- matrix(nrow=length(from.feat), ncol=B)

  #set the progress bar to be used later on within the loop
  pb <- txtProgressBar(min = 0, max = B, style = 3)

  #if both the from and to features are a point pattern, AND if the type parameter is set to 'perm', perform a randomized test
  #which, unlike the above section of the code, does not create random points within the study area but randomly assigns the
  #points membership to either the from or the to feature
  if ((class(from.feat)[1] == "SpatialPointsDataFrame" | class(from.feat)[1] == "SpatialPoints") &
      (class(to.feat)[1] == "SpatialPointsDataFrame" | class(to.feat)[1] == "SpatialPoints")
      & type=="perm") {

    #size of from.feat dataset
    nfrom <- length(from.feat)

    #size of to.feat dataset
    nto <- length(to.feat)

    #pool the coordinates of points in pattern from.feat and to.feat to be used later within the pointDistance function
    pooledData <- as.matrix(rbind(coordinates(from.feat), coordinates(to.feat)))

    #set the progress bar to used inside the loop
    pb <- txtProgressBar(min = 0, max = B, style = 3)

    #loop to perform permuted distance calculations
    for (i in 1:B) {
      #create a random set of number to be used as index to extract the values from the pooledData matrix
      #the set of number has size equal to the number of points in pattern x
      index <- sample(1:(nfrom + nto), size = nfrom, replace = F)
      #extract from the pooledData matrix the rows corresponding to the generated random index
      from.perm <- pooledData[index, ]
      #extract from the pooledData matrix all the rows that does not correspond to the random index
      #the procedure eventually creates two sets of points randomly shuffling the entries of the pooledData matrix
      to.perm <- pooledData[-index, ]
      #calculate the average minimum distance using the permuted point patterns
      res <- pointDistance(from.perm, to.perm, lonlat = FALSE, allpairs = TRUE)
      #calculate all the pair-wise distances between each permuted from-feature and to-feature
      dist.rnd.mtrx[,i] <- apply(res, 1, min)
      setTxtProgressBar(pb, i)
    }

    #if the preceding condition does not hold, ie if both datasets are not of point type
    #or of the type parameter is not set to perm....
    } else {

    #loop to perform randomized distance calculations
    for (i in 1:B){
      #draw a random sample of points within the study region
      rnd <- sp::spsample(region, n=length(from.feat), type='random')
      #calculate all the pair-wise distances between each random from-feature and to-feature
      rnd.dst <- gDistance(rnd, to.feat, byid=TRUE)
      #for each random from-feature (i.e., column-wise), get the minimum distance to the to-feature
      dist.rnd.mtrx[,i]  <- apply(rnd.dst, 2, min)
      setTxtProgressBar(pb, i)
    }
  }

  # Make a list for the ecdfs of the randominzed minimum distances
  rnd.ecdfs <- list()
  for(i in 1:ncol(dist.rnd.mtrx)){
    rnd.ecdfs[[i]] <- ecdf(dist.rnd.mtrx[,i])
  }

  #set the limit of the plot's a-axis
  #the lower limit is the minimum of the randomized minimum distances or of the observed minimum distances, whichever is smaller
  #the upper limit is the max of the observed minimum distances
  xlim = c(min(min(dist.rnd.mtrx), min(obs.min.dst)), max(obs.min.dst))

  # we will evaluate the ecdfs on a grid of 1000 points between
  # the x limits
  xs <- seq(xlim[1], xlim[2], length.out = 1000)

  # this actually gets those evaluations and puts them into a matrix
  out <- lapply(seq_along(rnd.ecdfs), function(i){rnd.ecdfs[[i]](xs)})
  tmp <- do.call(rbind, out)

  # get the .025 and .975 quantile for each column
  # at this point each column is a fixed 'x' and the rows
  # are the different ecdfs applied to that
  lower <- apply(tmp, 2, quantile, probs = .025)
  upper <- apply(tmp, 2, quantile, probs = .975)

  #set the plot's subtitle title
  if ((class(from.feat)[1] == "SpatialPointsDataFrame" | class(from.feat)[1] == "SpatialPoints") &
      (class(to.feat)[1] == "SpatialPointsDataFrame" | class(to.feat)[1] == "SpatialPoints")
      & type=="perm") {
  subtitle <- paste0("confidence envelope based on a permutation-based routine employing ", B, " iterations")
  } else {
    subtitle <- paste0("confidence envelope based on a randomization-based routine employing ", B, " iterations")
  }

  # plot the original data
  # ECDF of the first random dataset
  plot(dst.ecdf,
       verticals=TRUE,
       do.points=FALSE,
       col="white",
       xlab="distance to the closest to.feature (d)",
       ylab="Cumulative (d)",
       main="Feature-to-feature distance analysis: \ncumulative distribution of observ. min. distances against 95% conf. envel.",
       sub=subtitle,
       cex.main=0.90,
       cex.sub=0.70,
       xlim=xlim)

  # add in the quantiles
  polygon(c(xs,rev(xs)), c(upper, rev(lower)), col = "#DBDBDB88", border = NA)
  plot(dst.ecdf,
       verticals=TRUE,
       do.points=FALSE,
       add=TRUE)

  # restore the original graphical device's settings
  par(mfrow = c(1,1))

}
