#' R function for refined Nearest Neighbor analysis of point patterns (G function)
#'
#' The function allows to perform the refined Nearest Neighbor analysis of point patterns by plotting the cumulative Nearest Neighbour distance,
#' along with a 95percent confidence envelope and a curve representing the expected cumulative distribution under the assumption of complete spatial randomness.\cr
#'
#' The function uses a randomized approach to build the confidence envelope, whereby cumulative distributions of average NN distances of random points are computed across B iterations (1000 by default).
#' In each iteration, a set of random points (with sample size equal to the number of points of the input feature) is drawn.\cr
#'
#' @param feature: feature dataset (of point type).
#' @param studyplot: shapefile (of polygon type) representing the study area; if not provided, the study area is internally worked out as the convex hull enclosing the input feature dataset.
#' @param buffer: add a buffer to the studyplot (0 by default); the unit depends upon the units of the input data.
#' @param B: number of randomizations to be used (1000 by default).
#' @param order: integer indicating the kth nearest neighbour (1 by default).
#' @keywords refNNA
#' @export
#' @examples
#' data(springs)
#' refNNa(springs) #produces a plot representing the cumulative nearest neighbour distance distribution; a confidence envelope based on 1000 ranodmized simulations is also shown.
#' @seealso \code{\link{NNa}}
#'
refNNa <- function (feature, studyplot=NULL, buffer=0, B=1000, order=1) {
  if(is.null(studyplot)==TRUE){
    ch <- rgeos::gConvexHull(feature)
    region <- rgeos::gBuffer(ch, width=buffer)
  } else {
    region <- studyplot
  }
  dst <- spatstat::nndist(coordinates(feature), k=order)                                 #for each point in the input feature dataset, calculate the distance to its nearest neighbor
  dst.ecdf <- ecdf(dst)                                                                  #calculate the ECDF of the observed NN distances
  dist.rnd.mtrx <- matrix(nrow=length(feature), ncol=B)                                  #create a matrix to store the distance of each random point to its nearest neighbor; each column correspond to a random set of points
  pb <- txtProgressBar(min = 0, max = B, style = 3)

  for (i in 1:B){
    rnd <- sp::spsample(region, n=length(feature), type='random')                        #draw a random sample of points within the study region
    dist.rnd.mtrx[,i] <- spatstat::nndist(coordinates(rnd), k=order)                     #calculate the NN distances of the random points and store them in the matrix (column-wise)
    setTxtProgressBar(pb, i)
  }

  # Make a list for the ecdfs
  rnd.ecdfs <- list()
  for(i in 1:ncol(dist.rnd.mtrx)){
    rnd.ecdfs[[i]] <- ecdf(dist.rnd.mtrx[,i])
  }

  xlim = c(min(min(dist.rnd.mtrx), min(dst)), max(max(dist.rnd.mtrx), max(dst)))

  # We will evaluate the ecdfs on a grid of 1000 points between
  # the x limits
  xs <- seq(xlim[1], xlim[2], length.out = 1000)
  # This actually gets those evaluations and puts them into a matrix
  out <- lapply(seq_along(rnd.ecdfs), function(i){rnd.ecdfs[[i]](xs)})
  tmp <- do.call(rbind, out)

  # Get the .025 and .975 quantile for each column
  # at this point each column is a fixed 'x' and the rows
  # are the different ecdfs applied to that
  lower <- apply(tmp, 2, quantile, probs = .025)
  upper <- apply(tmp, 2, quantile, probs = .975)

  # Calculate the expected distribution
  d <- 0:max(dst)
  lambda <- length(feature) / gArea(region)
  E <- 1 - exp(-1 * lambda * pi * d^2)

  # Plot the original data
  plot(dst.ecdf,                                                                          #plot the ECDF of the first random dataset
       verticals=TRUE,
       do.points=FALSE,
       col="white",
       xlab="Nearest Neighbor distance (d)",
       ylab="G (d)",
       main="Refined Nearest Neighbor analysis (G function)",
       cex.main=0.95,
       xlim= xlim)

  # Add in the quantiles.
  polygon(c(xs,rev(xs)), c(upper, rev(lower)), col = "#DBDBDB88", border = NA)
  lines(d, E, col="red")
  plot(dst.ecdf,
       verticals=TRUE,
       do.points=FALSE,
       add=TRUE)
}
