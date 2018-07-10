#' R function for calculating the Hodder-Okell's A index of spatial association
#'
#' The function allows to calculate the Hodder-Okell's A index of spatial association between the features of two point patterns.\cr
#'
#' The functions takes as input two point patterns (SpatialPointDataframe class) and calculate the A index. Details about the latter are provided by:\cr
#' Orton C. 1980, "Mathematics in Archeology", Glasgow: William Collins Sons & Co Ltd, pp. 154-155\cr
#' Blankholm P. 1990, "Intrasite spatial Analysis in Theory and Practice", Aarhus: Aarhus University Press, pp. 130-135.\cr
#'
#' The A index is about equal to 1 when the two patterns are randomly mingled; it is smaller than 1 when the two patterns are segregrated; it is larger than 1 when the features
#' of the two point patterns tend to occur together. The computational details are provided by Blankholm's book cited above (page 132).\cr
#'
#' The significance of the A index is calculated via the randomized approach devised by:\cr
#' Kintigh K W. 1990, “Intrasite Spatial Analysis: A Commentary of Major Methids”.
#' In Voorrips A, “Mathematics and Information Science in Archaeology: A Flexible Framework”, Studies in Modern Archaeology 3: 165-200\cr
#'
#' Given two patterns A and B being analysed, the procedure keeps the points location unchanged and randomly assigns the points to either pattern.
#' The random re-assigment is performed B times (199 by default) and each time the A index is calculated. One-tailed and two-tailed p values are calculated
#' following the procedure described by Baddeley et al., "Spatial Point Patterns. Methodology and Applications with R", CRC Press 2016, p. 387.\cr
#'
#' The function produces:\cr
#' -an histogram showing the frequency distribution of the randomized A index, with vertical reference lines representing the 0.025th and 0.975th quantile
#' of the distribution. A black dot represents the observed A index. At the bottom of the chart the randomized p values are reported;\cr
#' -optionally (setting the 'addmap' parameter to TRUE), a map showing the point patterns (and the study area, if supplied).
#' @param x: point pattern (SpatialPointDataframe class).
#' @param y: point pattern (SpatialPointDataframe class).
#' @param studyplot: feature (of polygon type; SpatialPolygonsDataFrame class) representing the study area;
#' if not provided, the study area is internally worked out as the bounding polygon based on the union the convex hulls of the x and y patterns.
#' This is only used for visualization purpose, should the user want to plot the two point patterns within the actualy study area.
#' @param addmap: FALSE (default) or TRUE if the user does not want or wants a map of the study area and of the two patterns to be displayed.
#' @keywords A index
#' @export
#' @examples
#' data("malta_polyg") # load a sample polygon
#' pA <- spsample(malta_polyg, n=30, type='random')  #create a set of 30 random points within the polygon
#' pB <- spsample(malta_polyg, n=40, type='random')  #create a set of 40 random points within the polygon
#' Aindex(pA,pB) # calculate the Hodder-Okell's A index for the two patterns
#' @seealso \code{\link{distRandSign}} , \code{\link{distCovarModel}} , \code{\link{pointsCovarModel}}
#'
Aindex <- function(x, y, studyplot=NULL, B=199, addmap=FALSE){

  if(is.null(studyplot)==TRUE){
    # build the convex hull of the union of the convex hulls of the two features
    ch <- gConvexHull(raster::union(gConvexHull(x), gConvexHull(y)))
    region <- ch
  } else {
    region <- studyplot
  }

  #size of patter x
  nx <- length(x)

  #interpoint distance matrix
  distxx <- gDistance(x, x, byid=TRUE)

  #sum the interpoint distance matrix
  sum.distxx <- sum(distxx)

  #calculation of the "average" interpoint distance for pattern x
  rxx <- sum.distxx / (nx * (nx-1))

  #size of patter y
  ny <- length(y)

  #interpoint distance matrix
  distyy <- gDistance(y, y, byid=TRUE)

  #sum the interpoint distance matrix
  sum.distyy <- sum(distyy)

  #calculation of the "average" interpoint distance for pattern y
  ryy <- sum.distyy / (ny * (ny-1))

  #matrix of the distance between points in pattern x to the points in pattern y
  distxy <- gDistance(x, y, byid=TRUE)

  #sum of the interpoint distance matrix
  sum.distxy <- sum(distxy)

  #calculation of the "average" interpoint
  rxy <- sum.distxy / (nx * ny)

  #Aindex
  Aind <- (rxx * ryy) / (rxy^2)

  #pool the coordinates of points in pattern x and y to be used later within the pointDistance function
  pooledData <- as.matrix(rbind(coordinates(x), coordinates(y)))

  #set the progress bar to used inside the loop
  pb <- txtProgressBar(min = 0, max = B, style = 3)

  #create an empty sloth for the permuted A index
  Aind.perm <- numeric(length=B)

  for(i in 1:B){
    #create a random set of number to be used as index to extract the values from the pooledData matrix
    #the set of number has size equal to the number of points in pattern x
    index <- sample(1: (nx + ny), size=nx, replace=F)
    #extract from the pooledData matrix the rows corresponding to the generated random index
    x.perm <- pooledData[index,]
    #extract from the pooledData matrix all the rows that does not correspond to the random index
    #the procedure eventually creates two sets of points randomly shuffling the entries of the pooledData matrix
    y.perm <- pooledData[-index,]

    #repeat all the steps needed to calculate the A index using the random samples
    distxx.perm <- pointDistance(x.perm, lonlat=FALSE, allpairs=TRUE)
    sum.distxx.perm <- sum(distxx.perm)
    rxx.perm <- sum.distxx.perm / (nx * (nx-1))

    distyy.perm <- pointDistance(y.perm, lonlat=FALSE, allpairs=TRUE)
    sum.distyy.perm <- sum(distyy.perm)
    ryy.perm <- sum.distyy.perm / (ny * (ny-1))

    distxy.perm <- pointDistance(x.perm, y.perm, lonlat=FALSE, allpairs=TRUE)
    sum.distxy.perm <- sum(distxy.perm)
    rxy.perm <- sum.distxy / (nx * ny)
    Aind.perm[i] <- (rxx.perm * ryy.perm) / (rxy.perm^2)

    setTxtProgressBar(pb, i)
  }

  #calculate the single tailed and the two-tailed p values
  p.lower <- (1 + sum (Aind.perm < Aind)) / (1 + B)
  p.upper <- (1 + sum (Aind.perm > Aind)) / (1 + B)
  two.tailed.p <- 2 * (min(p.lower, p.upper))

  if(addmap==TRUE){
    par(mfrow=c(1,2))
    plot(region,
         main="Map of point patterns plus study area",
         cex.main=0.9, col=NA,
         border="red",
         lty=2,
         axes=TRUE)
    plot(x,
         add=TRUE,
         pch=20)
    plot(y,
         add=TRUE,
         pch=20,
         col="red")
  }

  #plot the histogram of the permuted distirbution of the A index
  hist(Aind.perm, main=paste0("Freq. distribution of Hodder-Okell's A index\n across ", B, " permutations"),
       sub=paste0("A-index: ", round(Aind, 3),"\n p-value segregated: ", p.lower, "; p-value associated: ", p.upper, "\np-value different from randomly mingled: ", two.tailed.p, "\nA<1 (segregated); A=1 (randomly minlged); A>1 (associated)"),
       xlab="",
       cex.main=0.90,
       cex.sub=0.70,
       cex.axis=0.85)
  rug(Aind.perm, col = "#0000FF")
  abline(v=quantile(Aind.perm, 0.025), lty=2, col="blue")
  abline(v=quantile(Aind.perm, 0.975), lty=2, col="blue")
  points(x=Aind, y=0, pch=20, col = "black")

  # restore the original graphical device's settings
  par(mfrow = c(1,1))
}
