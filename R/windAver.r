#' R function for averaging wind speed and direction
#'
#' The function provides the facility to average wind speed and direction data stored in two or more rasters ('RasterLayer' class). Of course, the input
#' rasters must have the same extent, resolution, and coordinate system.\cr
#'
#' The wind and direction rasters must be fed into the function as a vector storing the rasters' name.
#' The wind direction data are averaged by first computing the u and v components, then averaging them separately, and eventually converting them back to
#' degrees.\cr
#'
#' The u component is derived using the following formula:\cr
#'
#'  \eqn{ -x * sin(2 * pi * y / 360) }\cr
#'
#'  where x and y are the wind speed and directions respectively;\cr
#'
#'  the v component is derived using the following formula:\cr
#'
#'   \eqn{ -x * cos(2 * pi * y / 360)  }\cr
#'
#'  After being averaged, the averaged u and v components are converted back to degrees using the formula below:\cr
#'
#'  \eqn{ (atan2(mean.u.comp, mean.v.comp) * 360/2/pi) + 180 }\cr
#'
#'  The whole described procedure follows Grange, S. K. (2014). Technical note: Averaging wind speeds and directions.
#'  Auckland. https://doi.org/10.13140/RG.2.1.3349.2006\cr
#'
#' The function produces:\cr
#' -two plots representing the average wind speed and direction;\cr
#' -two GTiff files (saved in the R working directory) for the average wind speed and direction.\cr
#'
#' The function also returns a list containing the following data:\cr
#' -$wind.speed.avrg: wind speed data ('RasterLayer' class);\cr
#' -$wind.dir.avrg: wind direction data ('RasterLayer' class).\cr
#'
#' @param speed.data: vector containing the name of the raster datasets ('RasterLayer' class) storing the wind speed data.
#' @param dir.data: vector containing the name of the raster datasets ('RasterLayer' class) storing the wind direction data.
#' @keywords wind averaging
#' @export
#' @examples
#' # simulate toy data
#' set.seed(12345)
#'
#' monthA_speed <- raster(xmn=0,xmx=4,ymn=0,ymx=4,res=1)
#' monthA_speed[] <- runif(16, 5,15)
#' plot(monthA_speed, main="month A speed")
#' text(monthA_speed)
#'
#' monthB_speed <- raster(xmn=0,xmx=4,ymn=0,ymx=4,res=1)
#' monthB_speed[] <- runif(16, 10,20)
#' plot(monthB_speed, main="month B speed")
#' text(monthB_speed)
#'
#' monthA_dir <- raster(xmn=0,xmx=4,ymn=0,ymx=4,res=1)
#' monthA_dir[] <- runif(16, 0,369)
#'plot(monthA_dir, main="month A dir")
#' text(monthA_dir)
#'
#' monthB_dir <- raster(xmn=0,xmx=4,ymn=0,ymx=4,res=1)
#' monthB_dir[] <- runif(16, 0,369)
#' plot(monthB_dir, main="month B dir")
#' text(monthB_dir)
#'
#' # calculate the average wind speed and direction for the two toy datasets
#' res <- windAver(speed.data= c(monthA_speed, monthB_speed), dir.data= c(monthA_dir, monthB_dir))
#'
#' @seealso \code{\link{monthlyWind}}
#'
windAver <- function (speed.data, dir.data){

  #create a raster stack from the input wind speed rasters
  speed.stack <- raster::stack(speed.data)

  #create a raster stack from the input wind direction rasters
  dir.stack <- raster::stack(dir.data)

  #build two empty list where we can store the data created by the following loop
  u.comp <- list()
  v.comp <- list()

  #calculate the u and v components that will be averaged in a subsequent step;
  #in each list, the first element of the list will store the component derived from the first wind speed raster and from the first wind direction raster;
  #the second element of the list will store the component derived form the second wind speed and from the second wind direction rasters, and so on
  for(i in 1:length(speed.data)){
    u.comp[[i]] <- raster::overlay(speed.stack[[i]], dir.stack[[i]], fun=function(x, y) -x * sin(2 * pi * y / 360) )
    v.comp[[i]] <- raster::overlay(speed.stack[[i]], dir.stack[[i]], fun=function(x, y) -x * cos(2 * pi * y / 360) )
  }

  #create raster stacks from the u and v component lists
  u.comp <- raster::stack(u.comp)
  v.comp <- raster::stack(v.comp)

  #calculate the average for the u and v components using the raster stacks
  mean.u.comp <- mean(u.comp)
  mean.v.comp <- mean(v.comp)

  #uses the averaged u and v components to get an average in degrees
  dir.average <- (atan2(mean.u.comp, mean.v.comp) * 360/2/pi) + 180

  #calculate the average wind speed using the wind speed raster stack
  speed.average <- mean(speed.stack)

  #plot the wind speed and direction rasters
  raster::plot(speed.average,
       main="Wind speed average",
       cex.main=0.90)

  raster::plot(dir.average,
       main="Wind direction average",
       cex.main=0.90)

  #create a list for the results
  results <- list("wind.speed.avrg"=speed.average,
                  "wind.dir.avrg"=dir.average)

  #export the raster as GTiff files
  raster::writeRaster(results$wind.speed.avrg, "wind_speed_avrg", format="GTiff", overwrite=TRUE)
  raster::writeRaster(results$wind.dir.avrg, "wind_direct_avrg", format="GTiff", overwrite=TRUE)

  return(results)

}
