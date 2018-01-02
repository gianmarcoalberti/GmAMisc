#' R function to download month-averaged wind data (from Global Forecast System (GFS) of the USA's National Weather Service (NWS))
#'
#' The function allows to download wind data from NOAA/NCEP Global Forecast System (GFS) Atmospheric Model colection, creating monthly averages.
#' The extent of the study area can be specified: (1) by entering the geographic coordinates; (2) on the basis of an input raster dataset representing the study area itself (for instance, a Digital Elevation Model); (3) using a country code (for a list of country codes see http://kirste.userpage.fu-berlin.de/diverse/doc/ISO_3166.html).
#' The function saves two .geotiff files in the computer's working directory, one representing the wind speed, the other the wind direction. In both cases, the values are the average of the wind speed and direction values in the study area across the days of the selected month, in the selected year. A plot is also returned in the R console.\cr
#' The function returns a list containing the following data:\cr
#' -$windMonth: stores the U and V components for each output grid cell (spatial resolution 0.5 degrees=50 Km)\cr
#' -$windMonthFit: stores the wind speed and direction for each output grid cell.\cr
#' The function builds upon the 'wind.dl()' function from Javier Fernández-López's package 'rWind'. The help provided by Dr Fernández-López in creating an earlier version of the 'monthlyWind()' function is gratefully acknowledged.
#' @param raster: raster dataset representing the study area.
#' @param country: code of the country for which monthly wind average has to be calculated.
#' @param year: selected year (from 2011 to current); 2015 by default.
#' @param month: selected month; 01 by default.
#' @param days: number of days featuring the selected month; 31 by default.
#' @param lon1: western longitude. If it is set to NULL (default), it will be derived from the input raster.
#' @param lon2: eastern longitude. If it is set to NULL (default), it will be derived from the input raster.
#' @param lat1: southern latitude. If it is set to NULL (default), it will be derived from the input raster.
#' @param lat2: northern latitude. If it is set to NULL (default), it will be derived from the input raster.
#' @keywords wind
#' @export
#' @examples
#' res <- monthlyWind(year=2014, month=12, days=31, lon1=-10, lon2=5, lat1=35, lat2=45) #download wind data for Spain region, averaging the values across the 31 days of December 2014
#' res <- monthlyWind(country="ESP", year=2014, month=12, days=31) #same as above, but using the country code ESP (=Spain).
monthlyWind <- function(raster, country=NULL, year=2015, month=01, days=31, lon1=NULL, lon2=NULL, lat1=NULL, lat2=NULL){
  if (is.null(lon1)) {
    if (is.null(country)) {
      #reproject the raster representing the study area using a geographic coordinate system
      reprojected_raster <- projectRaster(raster, crs = "+proj=longlat +datum=WGS84")
      #retrieve the geographic coordinates of the reprojected raster's extent
      lon1 <- reprojected_raster@extent[1]
      lon2 <- reprojected_raster@extent[2]
      lat1 <- reprojected_raster@extent[3]
      lat2 <- reprojected_raster@extent[4]
    } else {
      #if a raster is not fed and a country code is input instead, get the extent of the area by first dowloading a DEM using 'getData()'
      #for a list of country codes see http://kirste.userpage.fu-berlin.de/diverse/doc/ISO_3166.html
      raster <- getData("alt", country = country)
      lon1 <- raster@extent[1]
      lon2 <- raster@extent[2]
      lat1 <- raster@extent[3]
      lat2 <- raster@extent[4]
    }
  }
  #create an empty list to store each day data
  wind_serie <- list()

  #create an object with times for which data is available. It will be used later
  times <- c(0,3,6,9,12,15,18,21)

  #create a loop into a loop; the internal loop will store daily data (one per "time") in a "wind_day" list, and will perform the daily mean.
  #the external loop will store each daily mean in a "wind_serie" list.
  for (d in 1:days){
    wind_day <- list()
    for (t in 1:8){
      w <- wind.dl(year,month,d,times[t],lon1,lon2,lat1,lat2)
      wind_day[[t]] <- w
    }
    wd <- wind.mean(wind_day)
    wind_serie[[d]] <- wd
  }

  #finally, you can compute month average with wind.mean:
  wind_month <- wind.mean(wind_serie)
  wind_month_fit <- wind.fit(wind_month)

  r_dir <- wind2raster(wind_month_fit, type="dir")
  r_speed <- wind2raster(wind_month_fit, type="speed")

  writeRaster(r_dir, paste0(month,"-", year,"_dir"), format="GTiff")
  writeRaster(r_speed, paste0(month, "-", year,"_sp"), format="GTiff")

  alpha <- arrowDir(wind_month_fit)
  plot(r_speed, main=paste0("wind speed and direction (", days,"-days average; month: ", month, "; year: ", year, ")"))
  newmap <- getMap(resolution = "low")
  lines(newmap, lwd=4)
  Arrowhead(wind_month_fit$lon, wind_month_fit$lat, angle=alpha, arr.length = 0.12, arr.type="curved")

  results <- list("windMonth"=wind_month, "windMonthFit"=wind_month_fit)
  return(results)
}
