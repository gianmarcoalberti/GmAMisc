#' R function for refined Nearest Neighbor analysis of point patterns (G function)
#'
#' The function allows to perform the refined Nearest Neighbor analysis of point patterns by plotting the cumulative Nearest Neighbour distance,
#' along with an acceptance interval (with significance level equal to 0.05;
#' sensu Baddeley et al., "Spatial Point Patterns. Methodology and Applications with R", CRC Press 2016, 208) based on B (set to 200 by default) realizations of a Complete Spatial Random process.
#' The function also allows to control for a first-order effect (i.e., influence of an underlaying numerical covariate) while performing the analysis.
#' The covariate must be of RasterLayer class.\cr
#'
#' The function uses a randomized approach to build the mentioned acceptance interval whereby cumulative distributions of average NN distances of random points are computed across B iterations.
#' In each iteration, a set of random points (with sample size equal to the number of points of the input feature) is drawn.\cr
#'
#' @param feature: feature dataset (of point type).
#' @param studyplot: shapefile (of polygon type) representing the study area; if not provided, the study area is internally worked out as the convex hull enclosing the input feature dataset.
#' @param buffer: add a buffer to the studyplot (0 by default); the unit depends upon the units of the input data.
#' @param B: number of randomizations to be used (200 by default).
#' @param cov.var: numeric covariate (of RasterLayer class) (NULL by default).
#' @param order: integer indicating the kth nearest neighbour (1 by default).
#' @keywords refNNA
#' @export
#' @examples
#' data(springs)
#' refNNa(springs) #produces a plot representing the cumulative nearest neighbour distance distribution; a confidence envelope based on 1000 ranodmized simulations is also shown.
#'
#' data(Starbucks) #load the Startbucks datset
#' data(popdensity) #load the raster representing the numerical covariate
#' refNNa(Starbucks, cov.var=popdensity) #perform the analysis, controlling for the 1st order effect
#' @seealso \code{\link{NNa}}
#'
refNNa <- function (feature, studyplot=NULL, buffer=0, B=200, cov.var=NULL, order=1) {

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
    #(if the covariate raster is provided; see above), then...
  } else {

    #tranform the cov.var from a RasterLayer to an object of class im, which is needed by spatstat
    cov.var.im <- as.im(cov.var)
  }

  #for each point in the input feature dataset, calculate the distance to its nearest neighbor
  dst <- spatstat::nndist(coordinates(feature), k=order)

  #calculate the ECDF of the observed NN distances
  dst.ecdf <- ecdf(dst)

  #create a matrix to store the distance of each random point to its nearest neighbor;
  #each column correspond to a random set of points
  dist.rnd.mtrx <- matrix(nrow=length(feature), ncol=B)

  #set the progress bar to be used later on within the loop
  pb <- txtProgressBar(min = 0, max = B, style = 3)

  #if there is no covariate data, draw random points within the study region
  if(is.null(cov.var)==TRUE){

    for (i in 1:B){
      #draw a random sample of points within the study region
      rnd <- sp::spsample(region, n=length(feature), type='random')
      #calculate the NN distances of the random points and store them in the matrix (column-wise)
      dist.rnd.mtrx[,i] <- spatstat::nndist(coordinates(rnd), k=order)
      setTxtProgressBar(pb, i)

    }

  } else {

    #if there is a covariate dataset
    for (i in 1:B){
      #draw random points via the spatstat's rpoint function,
      #using the covariate dataset as spatial covariate
      rnd <- spatstat::rpoint(n=length(feature), f=cov.var.im)
      #calculate the NN distances of the random points and store them in the matrix (column-wise)
      dist.rnd.mtrx[,i] <- spatstat::nndist(rnd, k=order)
      setTxtProgressBar(pb, i)
    }
  }

  # Make a list for the ecdfs
  rnd.ecdfs <- list()
  for(i in 1:ncol(dist.rnd.mtrx)){
    rnd.ecdfs[[i]] <- ecdf(dist.rnd.mtrx[,i])
  }

  xlim = c(min(min(dist.rnd.mtrx), min(dst)), max(dst))

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

  #adjust the main plot title according to the presence of a covariate dataset
  if(is.null(cov.var)==TRUE){
    maintitle <- paste0("Refined Nearest Neighbor analysis (G function) \n(acceptance interval based on ", B, " randomized iterations)")
  } else {
    maintitle <- paste0("Refined Nearest Neighbor analysis (G function, controlling for the covariate effect) \n(acceptance interval based on ", B, " randomized iterations)")
  }

  # Plot the original data
  # plot the ECDF of the first random dataset
  plot(dst.ecdf,
       verticals=TRUE,
       do.points=FALSE,
       col="white",
       xlab="Nearest Neighbor distance (d)",
       ylab="G (d)",
       main=maintitle,
       cex.main=0.90,
       xlim= xlim)

  # Add in the quantiles
  polygon(c(xs,rev(xs)), c(upper, rev(lower)), col = "#DBDBDB88", border = NA)
  plot(dst.ecdf,
       verticals=TRUE,
       do.points=FALSE,
       add=TRUE)
}
