#' R function to test the relationship of a set of points with the Thiessen tessellation built around points belonging to another feature dataset
#'
#' The function can be considered as a special case of the scenario "b" tested by the 'pointsInPolygons()' function provided by this same package, with the exception that in this case the polygons are not entered by the user but are internally created by the function.\cr
#'
#' The question this function may allow to address is: do the points belonging to a feature dataset tend to occur close to any of the points in another feature dataset than expected if the points would be randomly scattered across the study area?
#' To help addressing this question, the function creates Thiessen polygons around the input 'to.feature' and then runs the 'pointsInPolygons()' function using its 'scenario b'.\cr
#' For further details, see the help documentation of the 'pointsInPolygons()' function.
#' @param from.feat: feature (of point type) whose spatial association with to-feature has to be assessed.
#' @param to.feat: feature (of point type) in relation to which the spatial association of the from-feature has to be assessed.
#' @param cex.text: modify the size of the labels in the returned plot.
#' @keywords association
#' @export
#' @examples
#' Example 1
#' data(locations)
#' data(events)
#' result <- pointsToPointsTess(events, locations)
#'
#' Example 2:
#' data(deaths)
#' data(pumps)
#' result <- pointsToPointsTess(deaths, pumps)
#'
pointsToPointsTess <- function (from.feat, to.feat, cex.text=0.7) {
  vor <- dismo::voronoi(to.feat)
  pointsInPolygons(point.feat=from.feat, polyg.feat=vor, scenario="b", cex.text=cex.text)
  plot(to.feat, add=TRUE)
}
