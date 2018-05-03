#' R function to model (and test) the dependence of a point pattern on the distance to another pattern
#'
#' The function is a wrapper for a number of functions out of the extremely useful 'spatstat' package (specifically, ppm(), cdf.test(), auc(), roc(), effectfun()).
#' It allows to test if there is a significant dependence of the input point pattern on a spatial covariate (first-order effect), the latter being the distance to another feature (of either point or line type).\cr
#' The function takes as input two datasets: a point patter (SpatialPointsDataFrame class) and a feature (either SpatialPointsDataFrame or SpatialLinesDataFrame class)
#' the distance to which is used as spatial covariate for the input point pattern.\cr
#'
#' The function fits a inhomogeneous Poisson point process (Alternative Model-H1) with the distance to the second feature entered by the user (cov.var parameter) used as spatial covariate.
#' In other words, the fitted alternative model is a Poisson point process with intensity of the point pattern as a loglinear function of the distance to the second pattern entered by the user
#' (see Baddeley et al., "Spatial Point Patterns. Methodology and Applications with R", CRC Press 2016, 309-313). The distance to the second feature is internally calculated
#' via the spatstat's 'distfun()' function.
#'
#' Also, the function fits a homogeneous Poisson point model (Null Model-H0, equivalent to Complete Spatial Randomness:
#' Baddeley et al., "Spatial Point Patterns. Methodology and Applications with R", CRC Press 2016, 305-306), that is used as comparison for the inhomogeneous point process model
#' in a Likelihood Ratio test (Baddeley et al., "Spatial Point Patterns. Methodology and Applications with R", CRC Press 2016, 334-335). A significant result, i.e. a low p-value,
#' suggests rejecting the Null Hypothesis of CSR in favour of the Alternative Hypothesis of a Poisson point process affected by a covariate effect (i.e., inhomogeneous intensity due to the influence of the covariate)
#' (Baddeley et al., "Spatial Point Patterns. Methodology and Applications with R", CRC Press 2016, 305). \cr
#'
#' The function returns a 4 plots, which can be arranged in just one visualization setting the parameter oneplot to TRUE:\cr
#'
#' -plot of the study area along with the point pattern of interest and the second feature entered by the user (whose distance is the spatial covariate);\cr
#'
#' -plot of the fitted intensity against the spatial covariate (Baddeley et al., "Spatial Point Patterns. Methodology and Applications with R", CRC Press 2016, 308);\cr
#'
#' -plot of the cumulative distribution of the covariate at the data points against the cumulative distribution of the covariate at all the spatial location within the
#' study area (rationale: Baddeley et al., "Spatial Point Patterns. Methodology and Applications with R", CRC Press 2016, 184-185);\cr
#'
#' -plot of the ROC curve, which help assessing the strenght of the dependence on the covariate (Baddeley et al., "Spatial Point Patterns. Methodology and Applications with R", CRC Press 2016, 187-188).\cr
#'
#' Setting the parameter Foxall to TRUE, the third plot will be replaced by the chart of the Foxall's J function, which is
#' another "useful statistic" when the covariate is the distance to a spatial pattern (Baddeley et al., "Spatial Point Patterns. Methodology and Applications with R", CRC Press 2016, 187, 282-284).
#' Values of J are uqual to 1 when the two patterns are independent random patterns; values <1 indicate that the input point pattern tends to be closer to the cov.var pattern than expected
#' for random points; values >1 indicate that the input point pattern avoid the cov.var pattern, i.e. the point pattern is more likely than random points to lie far away from the cov.var pattern
#' (see Baddeley et al., "Spatial Point Patterns. Methodology and Applications with R", CRC Press 2016, 284).
#'
#' A list is also returned, containing what follows:\cr
#'
#' -$H0-model: info and relevant statistics regarding the Null Model;\cr
#' -$H1-model: info and relevant statistics regarding the Alternative Model;\cr
#' -$Model comparison (LRT): results of the Likelihood Ratio test;\cr
#' -$AIC-H0: AIC of the Null Model;\cr
#' -$AIC-H1: AIC of the Atlernative Model;\cr
#' -$KS test: information regarding the cumulative distribution comparison via Kolmogorov-Smirnov test;\cr
#' -$AUC: the AUC statistics.
#' @param feature: feature (of point type; SpatialPointsDataFrame class) representing the spatial point pattern of interest.
#' @param cov.var: feature (of either point or line type; SpatialPointsDataFrame or SpatialLinesDataFrame class) the distance to which represents the spatial covariate.
#' @param studyplot: feature (of polygon type; SpatialPolygonsDataFrame) representing the study area; if not provided, the study area is internally worked out as the bounding polygon based on the union the convex hulls of the feature and of the cov.var data.
#' @param buffer: add a buffer to the convex hull of the study area (0 by default); the unit depends upon the units of the input data.
#' @param oneplot: set to TRUE (default), will plot the charts into a single visualization.
#' @param Foxall: set to TRUE, will plot the Foxall's J function.
#' @keywords distance
#' @export
#' @examples
#' data(faults)
#' data(springs)
#' data(malta_polyg)
#' results <- distCovarModel(springs, faults, malta_polyg)
#' @seealso \code{\link{distRandSign}} , \code{\link{Aindex}} , \code{\link{pointsCovarModel}}
#'
distCovarModel <- function(feature, cov.var, studyplot=NULL, buffer=0,Foxall=FALSE, oneplot=FALSE){
  options(scipen=999)

  if(is.null(studyplot)==TRUE){
    #union() requires raster; build the convex hull of the union of the convex hulls of the two features
    ch <- gConvexHull(raster::union(gConvexHull(feature), gConvexHull(cov.var)))

    #add a buffer to the convex hull, with width set by the buffer parameter;
    #the unit depends upon the units of the input data
    region <- gBuffer(ch, width=buffer)

  } else {
    #if the studyplot is entered by the user, the region of interest is studyplot parameter
    region <- studyplot
  }

  #set the analytical window (needed by spatstat) to the extent of the region variable
  W  <- as(region, "owin")

  #transform the feature variable to a ppp class (needed by spatstat) and remove appended data by using ummark()
  feature.ppp <- unmark(as.ppp(feature))

  #give the feature dataset the same analytical window of the study region
  Window(feature.ppp) <- W

  #tranform the cov.var dataset to a ppp class or to a psp class (both needed by spatstat)
  #according to whether the cov.var is a SpatialPointsDataframe (or SpatialPoints) class, or a SpatialLinesDataFrame
  ifelse(class(cov.var)[1]=="SpatialPointsDataFrame" | class(cov.var)[1]=="SpatialPoints",
         cov.var.ppp <- unmark(as.ppp(cov.var)),
         cov.var.ppp <- unmark(as.psp(cov.var)))

  #give the cov.var dataset the same analytical window of the study region
  Window(cov.var.ppp) <- W

  #fit the Null Model (homegeneous point process)
  PPM0 <- ppm(feature.ppp ~ 1)

  #calculate the distance from the cov.var
  distance.fun <- distfun(cov.var.ppp)

  #fit the Alternative Model (inhomogeneous point process with intensity as function of the cavariate)
  PPM1 <- ppm(feature.ppp ~ distance.fun)

  #perform the KS test and calculate the Foxall's J function
  kolmsmirn <- cdf.test(feature.ppp, distance.fun)
  Jfunction <- Jfox(feature.ppp, cov.var.ppp)

  #calculate the AUC
  areaundercurve <- auc(feature.ppp, distance.fun, high=FALSE)

  #compare the models via likelihood ratio test
  model.comp <- anova(PPM0, PPM1, test="LRT")

  #set the output of the graphic device according to the oneplot parameter
  if(oneplot==TRUE){
    par(mfrow=c(2,2))
  } else {}

  #extract the p-value of the likelihod ratio test to be used in the plot subtitle
  anova.p <- model.comp$"Pr(>Chi)"[2]

  #classify the anova p value to be used in the plot subtitle
  anova.p.to.report <- ifelse(anova.p < 0.001, "< 0.001",
                              ifelse(anova.p < 0.01, "< 0.01",
                                     ifelse(anova.p < 0.05, "< 0.05",
                                            round(anova.p, 3))))

  #plot the study region along with the point pattern and the feature whose distance is used as covariate
  plot(cov.var, main="Study region, point pattern, and target pattern", cex.main=0.75,
       sub=paste0("Null Hypothesis (H0): Homogeneous Poisson process model\nAlternative Hyphotesis (H1): Inhomogeneous Poisson process model (intensity as loglinear function of the distance to the target pattern)\nH1 ", ifelse(anova.p > 0.05, "is not", "is"), " a significant improvement over H0 (Likelihood Ratio p-value: ", anova.p.to.report, "; AUC: ", round(areaundercurve,3), ")"), cex.sub=0.70)

  plot(feature, add=TRUE, pch=20, col="red")
  plot(region, add=TRUE)

  #plot the fitted Alternative Model
  #i.e., modelled intensity against the covariate
  plot(effectfun(PPM1, names(PPM1$covariates), se.fit=TRUE),
       main="Fitted intensity of the point pattern \nas (loglinear) function of the 'distance' covariate",
       cex.main=0.8,
       cex.axis=0.7,
       cex.lab=0.8, legend=TRUE)

  #if Foxall=FALSE plot the cumulative distr chart, otherwise plot the Foxall's J function
  if(Foxall==FALSE){
    plot(kolmsmirn, cex.main=0.8)
  } else {
    plot(Jfunction, cex.main=0.8)
  }

  #plot the ROC curve
  plot(roc(PPM1),
       main=paste0("ROC curve of the fitted intensity of point patter \nas (loglinear) function of the 'distance' cavariate \nAUC: ", round(areaundercurve,3)),
       cex.main=0.8)

  #create a list to store relevant results
  results <- list("H0-model"=PPM0,
                  "H1-model"=PPM1,
                  "Model comparison (LRT)"=model.comp,
                  "AIC-H0"=AIC(PPM0),
                  "AIC-H1"=AIC(PPM1),
                  "KS test"=kolmsmirn,
                  "AUC"=areaundercurve)
  return(results)
}
