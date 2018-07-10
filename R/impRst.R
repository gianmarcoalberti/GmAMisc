#' R function to easily import a raster dataset into R
#'
#' The function is a wrapper for the 'raster()' function out of the 'raster' package. It provides the facility to import
#' a raster dataset ('RasterLayer' class) by means of a window that allows the user to navigate through the computer's folders and to select the appropriate file.
#' @param none: no parameter is used within the function.
#' @keywords raster
#' @export
#' @examples
#' my.raster <- impRst() #a window will pop up allowing the user to select the raster dataset
#'
impRst <- function (){
  my.raster <- raster(file.choose())
  return(my.raster)
}
