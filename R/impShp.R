#' R function to easily import a vectorial dataset (shapefile) into R
#'
#' The function is a wrapper for the 'shapefile()' function out of the 'raster' package. It provides the facility to import
#' a vectorial dataset (of shapefile type) by means of a window that allows the user to navigate through the computer's folders and to select the appropriate file.
#' @param none: no parameter is used within the function.
#' @keywords shapefile
#' @export
#' @examples
#' my.shapefile <- impShp() #a window will pop up allowing the user to select the shapefile
#'
impShp <- function (){
  my.shapef <- shapefile(file.choose())
  return(my.shapef)
}
