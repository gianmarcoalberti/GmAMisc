#' R function for calculating accumulated cost of movement across the terrain and least-cost paths from an origin.
#'
#' The function provides the facility to calculate the accumulated cost of movement around a starting location and to optionally calculate least-cost paths toward
#' one or multiple destinations. The function implements different cost estimations directly or inderectly related to human movement across the landscape.
#' The function takes as input a Digital Terrain Model (RasterLayer class) and a point feature (SpatialPointsDataFrame class), the latter representing
#' the starting location, i.e. the location from which the accumulated cost is calculated. \cr
#'
#' If the parameter 'destin' is fed with a dataset representing destination location(s) (SpatialPointsDataFrame class), the function will also calculate
#' least-cost path(s) plotted on the input DTM; the length of each path will be saved under the variable 'length' stored in the 'LCPs' dataset (SpatialLines class) returned by the function.
#' The red dot(s) representing the destination location(s) will be labelled with numeric values representing
#' cost value at the location(s). The cost value will be also appended to the updated destination dataset returned by the function and
#' storing a new variable named 'cost'.\cr
#'
#' The function builds on functions out of the 'gdistance' package, and by default uses a 16-directions movement in calculating the accumulated cost-surface.
#' The number of movements can be set by the user via the 'moves' parameter. As noted in Nakoinz-Knitter (2016). "Modelling Human Behaviour in Landscapes". New York: Springer, p. 183,
#' "gdistance works with conductivity rather than the more usual approach using costs, we need inverse cost functions". For this reason, in this function the cost is calculated
#' using the inverse of the published cost functions.\cr
#'
#' The following cost functions are implemented:\cr
#'
#' \strong{Tobler's hiking function (on-path)}: reshaped as follows \cr
#'
#' \eqn{((6 * exp(3.5 * abs(slope + 0.05))) * 0.278)^-1}\cr
#'
#' if we use speed, the final accumulated values will be 1/travel time (according to the 'gdistance' help documentation; see page 16 of this PDF: https://cran.r-project.org/web/packages/gdistance/vignettes/gdistance1.pdf);
#' therefore, we use the reciprocal of speed to eventually get travel time/1. The Tobler's equation is multiplied by 0.278 (which is the ratio between 1000 -meters in 1 km- and 3600 -seconds in 1 hour-) to turn KmH into m/s,
#' and then reciprocated to turn m/s to s/m; before being reciprocating, we drop the minus before the 3.5. The same applies to the other Tobler's function-related
#' equations listed below.\cr
#'
#'
#' \strong{Tobler's hiking function (off-path)}:\cr
#'
#' \eqn{(((6 * exp(3.5 * abs(slope + 0.05))) * 1.666667) * 0.278)^-1}\cr
#'
#' as per Tobler's indication, the off-path walking speed is reduced by 0.6; we use the reciprocal (1.666667) in our case since we are dealing with pace instead of speed.\cr
#'
#'
#' \strong{Márquez-Pérez et al.'s modified Tobler hiking function}:\cr
#'
#' \eqn{((4.8 * exp(5.3 * abs((slope * 0.7) + 0.03))) * 0.278)^-1}\cr
#'
#' modified version as proposed by Joaquín Márquez-Pérez, Ismael Vallejo-Villalta & José I. Álvarez-Francoso (2017), "Estimated travel time for walking trails in natural areas",
#' Geografisk Tidsskrift-Danish Journal of Geography, 117:1, 53-62, DOI: 10.1080/00167223.2017.1316212.\cr
#'
#'
#' \strong{Relative energetic expenditure cost function}:\cr
#'
#' \eqn{(tan((atan(abs(slope) * 57.29578) * 0.0174532925)  / tan(1 * 0.0174532925)))^-1}\cr
#'
#' slope-based cost function expressing change in potential energy expenditure; in the above formula, atan(abs(slope)) * 57.29578) turns rise-over-run to degrees; multiplying by 0.0174532925 turns degrees to radians before calculating the tangent; the same applies to the degrees in the denominator. \cr
#' \strong{See} Conolly, J., & Lake, M. (2006). Geographic Information Systems in Archaeology. Cambridge: Cambridge University Press, p. 220;
#' \strong{see also} Newhard, J. M. L., Levine, N. S., & Phebus, A. D. (2014). The development of integrated terrestrial and marine pathways in the Argo-Saronic region, Greece. Cartography and Geographic Information Science, 41(4), 379–390, with references to studies that use this
#' function.\cr
#'
#'
#' \strong{Herzog's metabolic cost function in J/(kg*m)}:\cr
#'
#' \eqn{(1337.8 * slope^6 + 278.19 * slope^5 - 517.39 * slope^4 - 78.199 * slope^3 + 93.419 * slope^2 + 19.825 * slope + 1.64)^-1}\cr
#'
#' \strong{see} Herzog, I. (2016). Potential and Limits of Optimal Path Analysis. In A. Bevan & M. Lake (Eds.), Computational Approaches to Archaeological Spaces (pp. 179–211). New York: Routledge.\cr
#'
#'
#' \strong{Wheeled-vehicle critical slope cost function}:\cr
#'
#' \eqn{(1 + ((abs(slope)*100) / sl.crit))^-1}\cr
#'
#' where sl.crit (=critical slope, in percent) is "the transition where switchbacks become more effective than direct uphill or downhill paths" and typically is in the range 8-16;
#' \strong{see} Herzog, I. (2016). Potential and Limits of Optimal Path Analysis. In A. Bevan & M. Lake (Eds.), Computational Approaches to Archaeological Spaces (pp. 179–211). New York: Routledge. \cr
#'
#'
#' \strong{Pandolf et al.'s metabolic energy expenditure cost function (in Watts)}:\cr
#'
#' \eqn{(1.5 * W + 2.0 * (W + L) * (L / W)^2 + N * (W + L) * (1.5 * V^2 + 0.35 * V * abs(slope*100)))^-1}\cr
#'
#' where W is the walker's body weight (Kg), L is the carried load (in Kg), V is the velocity in m/s, N is a coefficient representing ease of movement on the terrain.\cr
#' As for the latter, suggested values available in literature are: Asphalt/blacktop=1.0; Dirt road=1.1; Grass=1.1; Light brush=1.2; Heavy brush=1.5; Swampy bog=1.8; Loose sand=2.1; Hard-packed snow=1.6; Ploughed field=1.3;
#' \strong{see} de Gruchy, M., Caswell, E., & Edwards, J. (2017). Velocity-Based Terrain Coefficients for Time-Based Models of Human Movement. Internet Archaeology, 45(45). https://doi.org/10.11141/ia.45.4.\cr
#' For this cost function, \strong{see} Pandolf, K. B., Givoni, B., & Goldman, R. F. (1977). Predicting energy expenditure with loads while standing or walking very slowly. Journal of Applied Physiology, 43(4), 577–581. https://doi.org/10.1152/jappl.1977.43.4.577.\cr
#' For the use of this cost function in a case study, \strong{see} Rademaker, K., Reid, D. A., & Bromley, G. R. M. (2012). Connecting the Dots: Least Cost Analysis, Paleogeography, and the Search for Paleoindian Sites in Southern Highland Peru. In D. A. White & S. L. Surface-Evans (Eds.), Least Cost Analysis of Social Landscapes. Archaeological Case Studies (pp. 32–45). University of Utah Press;
#' \strong{see also} Herzog, I. (2013). Least-cost Paths - Some Methodological Issues, Internet Archaeology 36 (http://intarch.ac.uk/journal/issue36/index.html) with references.\cr
#'
#'
#' \strong{Van Leusen's metabolic energy expenditure cost function (in Watts)}:\cr
#'
#' \eqn{(1.5 * W + 2.0 * (W + L) * (L / W)^2 + N * (W + L) * (1.5 * V^2 + 0.35 * V * abs(slope*100 + 10)))^-1}\cr
#'
#' which modifies the Pandolf et al.'s equation; \strong{see} Van Leusen, P. M. (2002). Pattern to process: methodological investigations into the formation and interpretation of spatial patterns in archaeological landscapes. University of Groningen.\cr
#' \strong{Note} that, as per Herzog, I. (2013). Least-cost Paths - Some Methodological Issues, Internet Archaeology 36 (http://intarch.ac.uk/journal/issue36/index.html) and
#' unlike Van Leusen (2002), in the above equation slope is expressed in percent and speed in m/s; also, in the last bit of the equantion, 10 replaces
#' the value of 6 used by Van Leusen.\cr
#'
#' When using the Tobler-related cost functions, the time unit can be selected by the user setting the 'time' parameter to 'h' (hour) or to 'm' (minutes).\cr
#'
#' In general, the user can also select which type of visualization the function has to produce; this is achieved setting the 'outp' parameter to either 'r' (=raster)
#' or to 'c' (=contours). The former will produce a raster image with a colour scale and contour lines representing the accumulated cost surface; the latter parameter will only
#' produce contour lines.\cr
#'
#' The contour lines' interval is set using the parameter 'breaks'; is not value is passed to the parameter, the interval will be set by default to
#' 1/10 of the range of values of the accumulated cost surface.\cr
#'
#' The function returns a list storing:\cr
#' -$accumulated.cost.raster: raster representing the accumualted cost (RasterLayer class);\cr
#' -$isolines: contour lines derived from the accumulated cost surface (SpatialLinesDataFrame class);\cr
#' -$LCPs: estimated least-cost paths (SpatialLines class);\cr
#' -$LCPs$length: length of each least-cost path (units depend on the unit used in the input DTM);\cr
#' -$dest.loc.w.cost: copy of the input destination location(s) dataset with a new variable ('cost') added.
#' @param dtm: digital terrain model (RasterLayer class).
#' @param origin: location from which the walking time is computed (SpatialPointsDataFrame class).
#' @param destin: location(s) to which least-cost path(s) is calculated (SpatialPointsDataFrame class).
#' @param funct: 't' (default) uses the on-path Tobler's hiking function; 'tofp' uses the off-path Tobler's hiking function; 'mt' uses the modified Tobler's function;
#' 'ree' uses the relative energetic expenditure cost function; 'hrz' uses the Herzog's metabolic cost function; 'wcs' uses the wheeled-vehicle critical slope cost function; 'p' uses the
#' Pandolf et al.'s metabolic energy expenditure cost function; 'vl' uses the Van Leusen's metabolic energy expenditure cost function (see Details).
#' @param time: time-unit expressed by the accumulated raster and by the isolines if a walking-pace cost function is used (Tobler's and Tobler-related cost functions);
#' 'h' for hour, 'm' for minutes.
#' @param outp: type of output: 'raster' or 'contours' (see Details).
#' @param sl.crit: critical slope (in percent), typically in the range 8-16 (10 by default) (used by the wheeled-vehicle cost function; see Details).
#' @param W: walker's body weight (in Kg; used by the Pandolf's or Van Leusen's cost function; see Details).
#' @param L: carried load weight (in Kg; used by the Pandolf's or Van Leusen's cost function; see Details).
#' @param N: coefficient representing ease of movement (1 by default) (used by the Pandolf's or Van Leusen's cost function; see Details).
#' @param V: speed in m/s (1.2 by default) (used by the Pandolf's or Van Leusen's cost function; see Details).
#' @param moves: number of directions used when computing the accumulated cost-surface (16 by default).
#' @param breaks: isolines interval; if no value is supplied, the interval is set by default to 1/10 of the range of values of the accumulated cost surface.
#' @param cex.breaks: set the size of the time labels used in the isochrones plot (0.6 by default).
#' @param cex.lcp.lab: set the size of the labels used in least-cost path(s) plot (0.6 by default).
#' @param oneplot: TRUE (default) or FALSE if the user wants or does not want the plots displayed in a single window.
#' @keywords accumulated cost
#' @export
#' @examples
#' data(volc) # load a sample Digital Terrain Model
#' data(volc.loc) # load a sample start location on the above DTM
#' data(destin_loc) # load the sample destination locations on the above DTM
#' res <- moveCost(dtm=volc, origin=volc.loc, destin=destin_loc, funct="t", time="h", outp="r", breaks=0.05) # calculate walking-time isochrones based on the on-path Tobler's hiking function,
#' setting the time unit to hours and the breaks (isochrones interval) to 0.05 hour; also, since destination locations are provided,
#' least-cost paths from the origin to the destination locations will be calculated and plotted
#'
#' # To compare two different sets of least-cost paths:
#' tobler <- moveCost(dtm=volc, origin=volc.loc, destin=destin_loc, funct="t", time="h", outp="r") #use the Tobler's on-path hiking cost function
#' wheeled <- moveCost(dtm=volc, origin=volc.loc, destin=destin_loc, funct="wcs", outp="r") #use the wheeled-vehicle critical slope cost function
#' plot(volc)
#' plot(volc.loc, add=TRUE, pch=20)
#' plot(destin_loc, add=TRUE, pch=20, col="red")
#' plot(tobler$LCPs, add=TRUE)
#' plot(wheeled$LCPs, add=TRUE, col="blue", lty=2)
#'
moveCost <- function (dtm=NULL, origin=NULL, destin=NULL, funct="t", time="h", outp="r", sl.crit=10, W=70, L=0, N=1, V=1.2, moves=16, breaks=NULL, cex.breaks=0.6, cex.lcp.lab=0.6, oneplot=TRUE){

  #suppress the warning messages to prevent a warning message from gdistance package to be returned;
  #as per gdistance vignette description, that specific warning message can be safely ignored
  options(warn=-1)

  #calculate the terrain slope
  altDiff <- function(x){x[2] - x[1]}
  hd <- transition(dtm, altDiff, moves, symm=FALSE)
  slope <- geoCorrection(hd)

  #create an index to identify adjacent cells
  adj <- adjacent(dtm, cells=1:ncell(dtm), pairs=TRUE, directions=moves)
  cost <- slope

  # apply the cost function to only those slope cells that are adjacents;
  # as for Tobler's function, if we use speed, the final accum values will be 1/travel time (according to gdistance vignette); therefore, we use the reciprocal of speed to eventually get travel time/1;
  # note: the Tobler's equation is multiplied by 0.278 (which is the ratio between 1000 [meters in 1 km] and 3600 [seconds in 1 hour]) to turn KmH into m/s,
  # and then reciprocated to turn m/s to s/m; before being reciprocating, we drop the minus before the 3.5 (see, e.g., Wikipedia's pace formula, at the voice Tolber's hiking function);
  # the same applies to the modified Tobler's equation; the use of the inverse applies to the other cost functions as well
  if (funct=="t") {
    #Tobler's hiking function in s/m
    cost[adj] <- ((6 * exp(3.5 * abs(slope[adj] + 0.05))) * 0.278)^-1

    #set the labels to be used within the returned plot
    main.title <- paste0("Walking-time isochrones (in ", time, ") around origin")
    sub.title <- "Walking-time based on the Tobler's on-path hiking function"
    legend.cost <- paste0("walking-time (", time,")")
    sub.title.lcp.plot <- paste0("LCP(s) and walking-time distance(s) based on the Tobler's on-path hiking function (time in ", time, ") \nblack dot=start location\n red dot(s)=destination location(s)")
  }

  if (funct=="tofp") {
    #Tobler's hiking function in s/m for off-path routes
    #note that the multiplier 1.666667 is the reciprocal of the multiplier 0.6 suggested by Tobler to reduce the off-path walking speed
    cost[adj] <- (((6 * exp(3.5 * abs(slope[adj] + 0.05))) * 1.666667) * 0.278)^-1

    #set the labels to be used within the returned plot
    main.title <- paste0("Walking-time isochrones (in ", time, ") around origin")
    sub.title <- "Walking-time based on the Tobler's off-path hiking function"
    legend.cost <- paste0("walking-time (", time,")")
    sub.title.lcp.plot <- paste0("LCP(s) and walking-time distance(s) based on the Tobler's off-path hiking function (time in ", time, ") \nblack dot=start location\n red dot(s)=destination location(s)")
  }

  if(funct=="mt") {
    #Márquez-Pérez et al.'s modified Tobler hiking function in s/m
    cost[adj] <- ((4.8 * exp(5.3 * abs((slope[adj] * 0.7) + 0.03))) * 0.278)^-1

    #set the labels to be used within the returned plot
    main.title <- paste0("Walking-time isochrones (in ", time, ") around origin")
    sub.title <- "Walking-time based on the Márquez-Pérez et al.'s modified Tobler hiking function"
    legend.cost <- paste0("walking-time (", time,")")
    sub.title.lcp.plot <- paste0("LCP(s) and walking-time distance(s) based on the Márquez-Pérez et al.'s modified Tobler hiking function (time in ", time, ") \nblack dot=start location\n red dot(s)=destination location(s)")
  }

  if(funct=="ree") {
    #relative energetic expenditure
    # atan(abs(slope[adj])) * 57.29578) turns rise-over-run to degrees
    # multiplying by 0.0174532925 turns degrees to radians before calculating tan; the same applies to the degrees in the denominator
    cost[adj] <- (tan((atan(abs(slope[adj]) * 57.29578) * 0.0174532925)  / tan(1 * 0.0174532925)))^-1

    #set the labels to be used within the returned plot
    main.title <- "Accumulated cost isolines around origin"
    sub.title <- "Cost based on the slope-based relative energetic expenditure cost function"
    legend.cost <- "cost"
    sub.title.lcp.plot <- paste0("LCP(s) and cost distance(s) based on the slope-based relative energetic expenditure cost function \nblack dot=start location\n red dot(s)=destination location(s)")
  }

  if(funct=="hrz") {
    #Herzog metabolic cost function in J/(kg*m)
    cost[adj] <- (1337.8 * slope[adj]^6 + 278.19 * slope[adj]^5 - 517.39 * slope[adj]^4 - 78.199 * slope[adj]^3 + 93.419 * slope[adj]^2 + 19.825 * slope[adj] + 1.64)^-1

    #set the labels to be used within the returned plot
    main.title <- "Accumulated cost isolines around origin"
    sub.title <- "Cost based on the Herzog's metabolic cost function \n cost in J / (Kg*m)"
    legend.cost <- "metabolic cost J / (Kg*m)"
    sub.title.lcp.plot <- paste0("LCP(s) and cost distance(s) based on the Herzog's metabolic cost function \ncost in J / (Kg*m) \nblack dot=start location\n red dot(s)=destination location(s)")
  }

  if(funct=="wcs") {
    #wheel critical slope cost function; slope is multiplid by 100 to make it in percent
    cost[adj] <- (1 + ((abs(slope[adj])*100) / sl.crit))^-1

    #set the labels to be used within the returned plot
    main.title <- "Accumulated cost isolines around origin"
    sub.title <- paste0("Cost based on the wheeled-vehicle critical slope cost function \ncritical slope set to ", sl.crit, " percent")
    legend.cost <- "cost"
    sub.title.lcp.plot <- paste0("LCP(s) and cost distance(s) based on the wheeled-vehicle critical slope cost function \ncritical slope set to ", sl.crit, " percent \nblack dot=start location\n red dot(s)=destination location(s)")
  }

  if(funct=="vl") {
    #Van Leusen's metabolic energy expenditure cost function
    #note: V is velocity in m/s
    cost[adj] <- (1.5 * W + 2.0 * (W + L) * (L / W)^2 + N * (W + L) * (1.5 * V^2 + 0.35 * V * abs(slope[adj]*100 + 10)))^-1

    #set the labels to be used within the returned plot
    main.title <- "Accumulated cost isolines around origin"
    sub.title <- paste0("Cost based on the Van Leusen's metabolic energy expenditure cost function \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V: ", V)
    legend.cost <- "energy expenditure cost (Megawatts)"
    sub.title.lcp.plot <- paste0("LCP(s) and cost distance(s) based on the Van Leusen's metabolic energy expenditure cost function \n cost in Megawatts; parameters: W: ", W, "; L: ", L, "; N: ", N, "; V: ", V, "\nblack dot=start location\n red dot(s)=destination location(s)")
  }

  if(funct=="p") {
    #Pandolf et al.'s metabolic energy expenditure cost function
    #note: V is velocity in m/s
    cost[adj] <- (1.5 * W + 2.0 * (W + L) * (L / W)^2 + N * (W + L) * (1.5 * V^2 + 0.35 * V * abs(slope[adj]*100)))^-1

    #set the labels to be used within the returned plot
    main.title <- "Accumulated cost isolines around origin"
    sub.title <- paste0("Cost based on the Pandolf et al.'s metabolic energy expenditure cost function \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V: ", V)
    legend.cost <- "energy expenditure cost (Megawatts)"
    sub.title.lcp.plot <- paste0("LCP(s) and cost distance(s) based on the Pandolf et al.'s metabolic energy expenditure cost function \n cost in Megawatts; parameters: W: ", W, "; L: ", L, "; N: ", N, "; V: ", V, "\nblack dot=start location\n red dot(s)=destination location(s)")
  }

  #geocorrection of the pace values
  Conductance <- geoCorrection(cost)

  #accumulate the pace outwards from the origin
  accum_final <- accCost(Conductance, coordinates(origin))

  #if user select the Tobler's or the modified Tobler's function, turn seconds into the user-defined time-scale
  if (funct=="t" | funct=="tofp" | funct=="mt") {
    if (time=="h") {
      #turn seconds into hours
      accum_final <- accum_final / 3600
    } else {
      #turn seconds into minutes
      accum_final <- accum_final / 60
    }
  }

  #if user select the Val Leusen's or the Pandolf et al.'s function, turn the cost from Watts to Megawatts
  if (funct=="vl" | funct=="p") {
    accum_final <- accum_final / 1000000
  }

  #if no break value is entered, set the breaks to one tenth of the range of the values of the final accumulated cost surface
  if(is.null(breaks)==TRUE){
    breaks <- round((max(values(accum_final)) - min(values(accum_final))) /10, 2)
  }

  #set the break values for the isolines
  levels <- seq(min(values(accum_final)), max(values(accum_final)), breaks)

  #conditionally set the layout in just one visualization
  if(is.null(destin)==FALSE & oneplot==TRUE){
    m <- rbind(c(1,2))
    layout(m)
  }

  #produce the output
  if (outp=="r") {
    #produce a raster with contours
    plot(accum_final,
         main=main.title,
         sub=sub.title,
         cex.main=0.95,
         cex.sub=0.75,
         legend.lab=legend.cost)
    contour(accum_final,
            add=TRUE,
            levels=levels,
            labcex=cex.breaks)
    plot(origin,
         pch=20,
         add=TRUE)

  } else {
    #only produce contours
    contour(accum_final,
            levels=levels,
            main=main.title,
            sub=sub.title,
            cex.main=0.95,
            cex.sub=0.75,
            labcex=cex.breaks)
    plot(origin,
         pch=20,
         add=TRUE)
  }

  #calculate and store the contours as a SpatialLinesDataFrame
  isolines <- rasterToContour(accum_final, levels=levels)

  #if 'destin' is NOT NULL, calculate the least-cost path(s) from the origin to the destination(s);
  #the 'Conductance' transitional layer is used
  if(is.null(destin)==FALSE){
    #calculate the least-cost path(s)
    sPath <- shortestPath(Conductance, coordinates(origin), coordinates(destin), output="SpatialLines")

    #plot the dtm
    plot(dtm, main="Digital Terrain Model with Least-cost Path(s)",
         sub=sub.title.lcp.plot,
         cex.main=0.90,
         cex.sub=0.7,
         legend.lab="Elevation (masl)")

     #add the origin
    plot(origin, add=TRUE, pch=20)

    #add the destination(s)
    plot(destin,
         add=TRUE,
         pch=20,
         col="red")

    #add the LCPs
    lines(sPath)

    #calculate the length of the least-cost paths and store the values by appending them to a new variable of the sPath object
    sPath$length <- gLength(sPath, byid=TRUE)

    #extract the cost from the accum_final to the destination location(s), appending the data to a new column
    destin$cost <- extract(accum_final, destin)

    #add the point(s)'s labels
    text(coordinates(destin),
         labels=round(destin$cost,2),
         pos = 4,
         cex=cex.lcp.lab)

  } else {
    sPath=NULL
    dest.loc.w.walking.time=NULL
  }

  #restore the original graphical device's settings if previously modified
  if(is.null(destin)==FALSE & oneplot==TRUE){
    par(mfrow = c(1,1))
  }

  #restore the warning messages
  options(warn=0)

  results <- list("accumulated.cost.raster"=accum_final,
                  "isolines" = isolines,
                  "LCPs"=sPath,
                  "dest.loc.w.cost"=destin)
}
