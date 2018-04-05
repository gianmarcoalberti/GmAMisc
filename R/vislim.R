#' R function for computing the limit of visibility of an object given its height
#'
#' The function allows to plot the angular size of an object (in degrees) against the distance from the observer, and to compute at which distance
#' from the observer the angular size of the object hits the limit of human visual acuity (0.01667 degrees).
#'
#' The function returns:\cr
#' -a plot displaying the decay in angular size as function of the object's distance from the observer; a black dot represents the distance at which the angular size hits the limit of human visual acuity;\cr
#' -the value (in km) of the visibility limit.
#' @param vis.degree: limit of human visual acuity (0.01667 by default).
#' @param targ.h: target size (=height in meters).
#' @export
#' @examples
#' limit <- vislim(targ.h=6) # calculate the visibility limit of an object of size 6m, and store the result (20.62 km) in the 'limit' variable
#'
vislim <- function(vis.degree=0.01667, targ.h) {
  deg2rad <- (vis.degree*pi)/180
  a <- 1/(2*tan(deg2rad/2))
  b2 <- a * (targ.h * 0.001)                                                 # use the calculated multiplier 'a' to multiply the target's height; the latter is first converted in km to keep the output distance on the same scale
  dataf <- data.frame(size=targ.h, dist=seq(1000, (b2/0.001)+1000, by=500))
  dataf$ang.size <- 57.3 * (dataf$size/dataf$dist)
  plot(dataf$dist, dataf$ang.size,
       type="l",
       xlab=paste0("distance of a ", targ.h, "-m-high object from the observer (m)"),
       ylab=paste0("target object (height: ", targ.h, "m) angular size (degrees)"),
       main="Object's distance from the observer vs. object's angular size",
       sub=paste0("red line: limit of visual acuity (", vis.degree, " degrees); object's visibility limit: ", round(b2*1000,2), " m"),
       cex.main=0.9,
       cex.sub=0.75)
  abline(h=vis.degree, col="red", lty=2)
  points(x=b2*1000, y=vis.degree, pch=20, col="black")
  return(b2)
}

