# GmAMisc (Gianmarco Alberti Miscellaneous)
vers 0.19

`GmAMisc` is a collection of functions that I have built in different points in time. The functions' aim spans from univariate outlier detection, to permutation t test, permutation chi-square test, calculation of Brainerd-Robinson similarity coefficient, validation of logistic regression models, point pattern analysis, and more. 

The package comes with some toy datasets:

`assemblage`: distribution of 7 archaeological objects across 9 assemblages.

`deaths`: SpatialPointsDataFrame representing the location of cholera deaths in London (after Dr Snow's mid-1800s study of cholera outbreak in Soho).

`destin_loc`: SpatialPointsDataFrame representing spots on the volcano Maunga Whau (Auckland, New Zealand), to be used as destination locations for least-cost paths calculation.

`events`: SpatialPointsDataFrame representing fictional events.

`faults`: SpatialLinesDataFrame representing the geological fault-lines in Malta.

`locations`: SpatialPointsDataFrame representing fictional locations.

`log_regr_data`: admission to graduate school.

`malta_dtm_40`: A RasterLayer representing a Digital Terrain Model for Malta (40m resolution).

`malta_polyg`: A SpatialPolygonsDataFrame representing Malta.

`Massachusetts`: A SpatialPolygonsDataFrame representing the limits of Massachusetts.

`phases`: posterior probabilities for the chronological relation of the Starting and Ending boundaries of two Bayesian independent 14C phases.

`radioc_data`: posterior probabilities for the Starting and Ending boundaries of two 14C phases.

`points`: SpatialPointsDataFrame representing fictional locations.

`polygons`: SpatialPolygonsDataFrame representing fictional polygons.

`popdensity`: A RasterLayer representing the population density in Massachusetts.

`pumps`: SpatialPointsDataFrame representing the location of public water pumps in London (after Dr Snow's mid-1800s study of cholera outbreak in Soho).

`rndpoints`: SpatialPointsDataFrame representing the random locations.

`springs`: SpatialPointsDataFrame representing the location of springs in Malta.

`Starbucks`: SpatialPointsDataFrame representing the location of Starbucks in Massachusetts.

`thiessenpolyg`: SpatialPolygonsDataFrame representing Thiessen polygons around the points represented in the 'locations' dataset.

`volc`: Digital Terrain Model representing the volcano Maunga Whau (Auckland, New Zealand).

`volc.loc`: A SpatialPointsDataFrame representing a spot on the volcano Maunga Whau (Auckland, New Zealand).

<br>

## List of implemented functions
* `Aindex()`: function for calculating the Hodder-Okell's A index of spatial association.
* `aucadj()`: function for optimism-adjusted AUC (Logistic Regression internal validation).
* `BRsim()`: function for Brainerd-Robinson simiarity coefficient.
* `chiperm()`: function for permutation-based chi-square test of independence.
* `distCovarModel()`: function to model (and test) the dependence of a point pattern on the distance to another pattern.
* `distRandCum()`: function to test the significance of the spatial relationship between two features in terms of the cumulative distribution of minimum distances.
* `distRandSign()`: function to test the significance of the spatial relationship between two features (points-to-points, points-to-lines, points-to-polygons).
* `featClust()`: function for points clustering on the basis of planar distance.
* `kwPlot()`: function for visually displaying Kruskal-Wallis test's results.
* `landfClass()`: function for landform classification on the basis od Topographic Position Index.
* `logregr()`: function easy binary Logistic Regression and model diagnostics.
* `modelvalid()`: function for binary Logistic Regression internal validation.
* `monthlyWind()`: function that allows to download wind data from NOAA/NCEP Global Forecast System (GFS) Atmospheric Model colection, creating monthly averages.
* `moveCost()`: function for calculating accumulated cost of movement across the terrain and least-cost paths from an origin.
* `mwPlot()`: function for visually displaying Mann-Whitney test's results.
* `NNa()`: function for Nearest Neighbor analysis of point patterns.
* `outlier()`: function for univariate outliers detection.
* `perm.t.test()`: function for permutation-based t-test.
* `plotJenks()`: function for plotting univariate classification using Jenks' natural break method.
* `pointsCovarModel()`: function to model (and test) the dependence of a point pattern on a spatial numeric covariate.
* `pointsInPolygons()`: function to test points-in-polygons relationship.
* `pointsToPointsTess()`: function to test the relationship of a set of points with the Thiessen tessellation built around points belonging to another feature dataset.
* `ppdPlot()`: function for plotting Posterior Probability Densities for Bayesian modeled 14C dates/parameters.
* `prob.phases.relat()`: function to calculate the Posterior Probability for different chronological relations between two Bayesian radiocarbon phases.
* `refNNa()`: function for refined Nearest Neighbor analysis of point patterns (G function).
* `resc.val()`: function to rescale the values of a datset between a minimum and a maximum set by the user.
* `robustBAplot()`: function to plot a robust version of the Bland-Altman plot.
* `vislim()`: function for computing the limit of visibility of an object given its height.

<br>

## Description of implemented functions
`Aindex()`: the function allows to calculate the Hodder-Okell's A index of spatial association between the features of two point patterns. It takes as input two point patterns (SpatialPointDataframe class) and calculate the A index. Details about the latter are provided by:
* Orton C. 1980, "Mathematics in Archeology", Glasgow: William Collins Sons & Co Ltd, pp. 154-155
* Blankholm P. 1990, "Intrasite spatial Analysis in Theory and Practice", Aarhus: Aarhus University Press, pp. 130-135.

The A index is about equal to 1 when the two patterns are randomly mingled; it is smaller than 1 when the two patterns are segregrated; it is larger than 1 when the features of the two point patterns tend to occur together. The computational details are provided by Blankholm's book cited above (page 132). The significance of the A index is calculated via the randomized approach devised by:
* Kintigh K W. 1990, “Intrasite Spatial Analysis: A Commentary of Major Methids”. In Voorrips A, “Mathematics and Information Science in Archaeology: A Flexible Framework”, Studies in Modern Archaeology 3: 165-200.

Given two patterns A and B being analysed, the procedure keeps the points location unchanged and randomly assign the points to either pattern.
The random re-assigment is performed B times (199 by default) and each time the A index is calculated. One-tailed and two-tailed p values are calculated following the procedure described by Baddeley et al., "Spatial Point Patterns. Methodology and Applications with R", CRC Press 2016, p. 387.

The function produces:
* an histogram showing the frequency distribution of the randomized A index, with vertical reference lines representing the 0.025th and 97.5th quantile of the distribution. A black dot represents the observed A index. At the bottom of the chart the randomized p values are reported;
* optionally (setting the 'addmap' parameter to TRUE), a map showing the point patterns (and the study area, if supplied).

<br>

`aucadj()`: function for optimism-adjusted AUC (Logistic Regression internal validation). The function allows to calculate the AUC of a (binary) Logistic Regression model, adjusted for optimism. The function performs an internal validation of a model via a bootstrap procedure (devised by Harrell and colleagues), which enables to estimate the degree of optimism of a fitted model and the extent to which the model will be able to generalize outside the training dataset.
The returned boxplots represent:
* the distribution of the AUC value in the bootstrap sample (auc.boot), which represents "an estimation of the apparent performance" (according to the aforementioned reference);
* the distribution of the AUC value deriving from the model fitted to the bootstrap samples and evaluated on the original sample (auc.orig), which represents the model performance on independent data.

At the bottom of the chart, the apparent AUC (i.e., the value deriving from the model fitted to the original dataset) and the AUC adjusted for optimism are reported.

<br>

`BRsim()`: function for Brainerd-Robinson simiarity coefficient. The function allows to calculate the Brainerd-Robinson similarity coefficient, taking as input a cross-tabulation (dataframe).
The function will return:
* a) a correlation matrix in tabular form;
* b) a heat-map representing, in a graphical form, the aforementioned correlation matrix.
In the heat-map (which is built using the 'corrplot' package), the size and the color of the squares are proportional to the Brainerd-Robinson coefficients, which are also reported by numbers.

In order to "penalize" BR similarity coefficient(s) arising from assemblages with unshared categories, the function does what follows: it divides the BR coefficient(s) by the number of unshared categories plus 0.5. The latter addition is simply a means to be still able to penalize coefficient(s) arising from assemblages having just one unshared category. Also note that joint absences will have no weight on the penalization of the coefficient(s). In case of assemblages sharing all their categories, the corrected coefficient(s) turns out to be equal to the uncorrected one.

<br>

`chiperm()`: function for permutation-based chi-square test of independence. The function performs the chi-square test of independence on the basis of permuted tables, whose number is selected by user. For the rationale of this approach, see for instance the nice description provided by Beh E.J., Lombardo R. 2014, *Correspondence Analysis: Theory, Practice and New Strategies*, Chichester, Wiley, pages 62-64.
The function produces:
* (1) a chart that displays the permuted distribution of the chi square statistic based on B permuted tables. The selected number of permuted tables, the observed chi square, the 95th percentile of the permuted distribution, and the associated p value are reported at the bottom of the chart; 
* (2) a chart that displays the bootstrap distribution of Cramer's V coefficient, based on a number of bootstrap replicates which is equal to the value of the function's parameter B; 
* (3) a chart that the Pearson's Standardized Residuals: a colour scale allows to easily understand which residual is smaller (BLUE) or larger (RED) than expected under the hypothesis of independence. Should the user want to only display residuals larger than a given threshold, it suffices to set the filter parameter to TRUE, and to specify the desidered threshold by means of the thresh parameter, which is set at 1.96 by default.

<br>

`distCovarModel()`: is a wrapper for a number of functions out of the extremely useful `spatstat` package (specifically, `ppm()`, `cdf.test()`, `auc()`, `roc()`, `effectfun()`). It allows to test if there is a significant dependence of the input point pattern on a spatial covariate (first-order effect), the latter being the distance to another feature (of either point or line type).The function takes as input two datasets: a point patter (`SpatialPointsDataframe` class) and a feature (either `SpatialPointsDataframe` or `SpatialLinesDataframe` class) the distance to which is used as spatial covariate for the input point pattern.

The function fits a inhomogeneous Poisson point process (Alternative Model-H1) with the distance to the second feature entered by the user (`cov.var` parameter) used as spatial covariate. In other words, the fitted alternative model is a Poisson point process with intensity of the point pattern as a loglinear function of the distance to the second pattern entered by the user (see Baddeley et al., "Spatial Point Patterns. Methodology and Applications with R", CRC Press 2016, 309-313). The distance to the second feature is internally calculated via the spatstat's `distfun()` function. Also, the function fits a homogeneous Poisson point model (Null Model-H0, equivalent to Complete Spatial Randomness: Baddeley et al., "Spatial Point Patterns. Methodology and Applications with R", CRC Press 2016, 305-306), that is used as comparison for the inhomogeneous point process model in a Likelihood Ratio test (Baddeley et al., "Spatial Point Patterns. Methodology and Applications with R", CRC Press 2016, 334-335). A significant result, i.e. a low p-value, suggests rejecting the Null Hypothesis of CSR in favour of the Alternative Hypothesis of a Poisson point process affected by a covariate effect (i.e., inhomogeneous intensity due to the influence of the covariate) (Baddeley et al., "Spatial Point Patterns. Methodology and Applications with R", CRC Press 2016, 305).

The function returns a 4 plots, which can be arranged in just one visualization setting the parameter `oneplot` to TRUE:
* plot of the study area along with the point pattern of interest and the second feature entered by the user (whose distance is the spatial covariate);
* plot of the fitted intensity against the spatial covariate (Baddeley et al., "Spatial Point Patterns. Methodology and Applications with R", CRC Press 2016, 308);
* plot of the cumulative distribution of the covariate at the data points against the cumulative distribution of the covariate at all the spatial location within the study area (rationale: Baddeley et al., "Spatial Point Patterns. Methodology and Applications with R", CRC Press 2016, 184-185);
* plot of the ROC curve, which help assessing the strenght of the dependence on the covariate (Baddeley et al., "Spatial Point Patterns. Methodology and Applications with R", CRC Press 2016, 187-188).

Setting the parameter Foxall to TRUE, the third plot will be replaced by the chart of the Foxall's J function, which is another "useful statistic" when the covariate is the distance to a spatial pattern (Baddeley et al., "Spatial Point Patterns. Methodology and Applications with R", CRC Press 2016, 187, 282-284). Values of J are uqual to 1 when the two patterns are independent random patterns; values <1 indicate that the input point pattern tends to be closer to the cov.var pattern than expected for random points); values >1 indicate that the input point pattern avoid the cov.var pattern, i.e. the point pattern is more likely than random points to lie far away from the cov.var pattern (see Baddeley et al., "Spatial Point Patterns. Methodology and Applications with R", CRC Press 2016, 284).

A list is also returned, containing what follows:
* $H0-model: info and relevant statistics regarding the Null Model;
* $H1-model: info and relevant statistics regarding the Alternative Model;
* $Model comparison (LRT): results of the Likelihood Ratio test;
* $AIC-H0: AIC of the Null Model;
* $AIC-H1: AIC of the Atlernative Model;
* $KS test: information regarding the cumulative distribution comparison via Kolmogorov-Smirnov test;
* $AUC: the AUC statistics.

<br>

`distRandCum()`: allows to assess if there is a significant spatial association between a point pattern and the features of another pattern.
For instance, users may want to assess if the features of a point pattern tend to lie close to some features represented by polylines. Given a from-feature (event for which we want to estimate the spatial association with the to-feature) and a to-feature (event in relation to which we want to estimate the spatial association for the from-feature), the assessment is performed by means of a randomized procedure:
* keeping fixed the location of the to-feature, random from-features are drawn B times (the number of randomized from-features is equal to the number of observed from-features);
* for each draw, the minimum distance to the to-features is calculated; if the to-feature is made up of polygons, the from-features falling within a polygon will have a distance of 0;
* a cumulative distribution of random minimum distances is thus obtained;
* the cumulative random minimum distances are used to work out an acceptance interval (with significance level equal to 0.05; sensu Baddeley et al., "Spatial Point Patterns. Methodology and Applications with R", CRC Press 2016, 208) that allows to assess the statistical significance of the
cumulative distribution of the observed minimum distances, and that is built using the above-mentioned B realizations of a Complete Spatial Random.

The from-feature must be a point feature, whilst the to-feature can be a point or a polyline or a polygon feature.

The rationale of the procedure is that, if there indeed is a spatial association between the two features, the from-feature should be closer to the to-feature than randomly generated from-features. If the studyplot shapefile is not provided, the random locations are drawn within a bounding polygon based on the union the convex hulls of the from- and of the to-feature.

If both the from-feature and the to-feature are of point type (SpatialPointsDataFrame class), the user may opt for the randomized procedure described above (parameter 'type' set to 'rand'), or for a permutation-based procedure (parameter 'type' set to 'perm'). Unlike the procedure described above, whereby random points are drawn within the study area, the permutation-based routine builds a cumulative distribution of minimum distances keeping the points location unchanged and randomly assigning the points to either of the two patterns. The re-assigment is performed B times (200 by default) and each time the minimum distance is calculated.

The function produces a cumulative distribution chart in which the distribution of the observed minimum distances is represented by a black line,
and the 95percent confidence envelope is represented in grey. The number of iteration used and the type of analysis (whether randomization-based or permutation-based) are reported in the chart's title.

<br>

`distRandSign()`: allows to assess if there is a significant spatial association between a point pattern and the features of another pattern. For instance, users may want to assess if some locations tend to lie close to some features represented by polylines. By the same token, users may want to know if there is a spatial association between the location of a given event and the location of another event.

Given a from-feature (event for which we want to estimate the spatial association with the to-feature) and a to-feature (event in relation to which we want to estimate the spatial association for the from-feature), the assessment is performed by means of a randomized procedure:
* keeping fixed the location of the to-feature, random from-features are drawn B times (the number of randomized from-features is equal to the number of observed from-features);
* for each draw, the average minimum distance to the to-features is calculated; if the to-feature is made up of polygons, the from-features falling within a polygon will have a distance of 0;
* a distribution of average minimum distances is thus obtained;
* p values are computed following Baddeley et al., "Spatial Point Patterns. Methodology and Applications with R", CRC Press 2016, p. 387.

The from-feature must be a point feature, whilst the to-feature can be a point or a polyline or a polygon feature.

The rationale of the procedure is that, if there indeed is a spatial association between the two features, the from-feature should be on average closer to the to-feature than randomly generated from-features. If the study plot shapefile is not provided, the random locations are drawn within a bounding polygon based on the union the convex hulls of the from- and of the to-feature.

If both the from-feature and the to-feature are of point type (SpatialPointsDataFrame class), the function also test the spatial association by means of a permuted procedures. Unlike the procedure described above, whereby random points are drawn within the study area, the permutation-based routine builds a distribution of averages minimum distances keeping the points location unchanged and randomly assigning the points to either of the two patterns. The re-assigment is performed B times (199 by default) and each time the average minimum distance is calculated.

The function produces a histogram showing: the distribution of randomized average minimum distances; a black dot indicating the observed average minimum distance; a hollow dot representing the average of the randomized minimum distances; two blue reference lines correspond to the 0.025th and to the 0.975th quantile of the randomized distribution. P-values are reported at the bottom of the plot. In case both the from- and the to- feature are of point type, another histogram is produced, which provides the same information of the preceding histogram, but derived from the permutation-based routine that has been detailed above.

A list is also returned, containing what follows:
* $from.feat.min.dist: distance of each entity of the from-feature to the nearest entity of the to-feature;
* $avrg.obs.min.dist: observed average minimum distance;
* $avrg.rnd.min.dist: average randomized minimum distance;
* $avrg.perm.min.dist: average permuted minimum distance (returned only when both the from- and to- features are of point type);
* $p.value closer than expected-rnd-;
* $p.value closer than expected-perm- (returned only when both the from- and to- features are of point type);
* $p.value more distant than expected-rnd-;
* $p.value more distant than expected-perm- (returned only when both the from- and to- features are of point type);
* $p.value different from random-rnd-;
* $p.value different from random-perm- (returned only when both the from- and to- features are of point type).

<br>

`featClust()`: provides the facility to cluster the features of the input dataset on the basis of either their (projected) coordinates (for points; SpatialPointsDataFrame class) or of their area (for polygons; SpatialPolygonsDataFrame class). The function internally calculates a distance matrix (based on the Euclidean Distance) on the basis of the points' coordinates or polygons' area. A dendrogram is produced which depicts the hierarchical clustering based (by default) on the Ward's agglomeration method; rectangles identify the selected cluster partition. Besides the dendrogram, a silhouette plot is produced, which allows to measure how 'good' is the selected cluster solution.

As for the latter, if the parameter `part` is left empty (default), an optimal cluster solution is obtained. The optimal partition is selected via an iterative procedure which locates at which cluster solution the highest average silhouette width is achieved. If a user-defined partition is needed, the user can input the desired number of clusters using the parameter `part`. In either case, an additional plot is returned besides the cluster dendrogram and the silhouette plot; it displays a scatterplot in which the cluster solution (x-axis) is plotted against the average silhouette width (y-axis). A black dot represent the partition selected either by the iterative procedure or by the user.

Notice that in the silhouette plot, the labels on the left-hand side of the chart show the point ID number and the cluster to which each point is closer. Also, the function returns a plot showing the input point dataset, with points colored by cluster membership. Two new variables are added to the point shapefile's dataframe, storing a point ID number and the corresponding cluster membership.

The silhouette plot is obtained from the silhouette() function out from the 'cluster' package (https://cran.r-project.org/web/packages/cluster/index.html).

For a detailed description of the silhouette plot, its rationale, and its interpretation, see:
* Rousseeuw P J. 1987. "Silhouettes: A graphical aid to the interpretation and validation of cluster analysis", Journal of Computational and Applied Mathematics 20, 53-65 (http://www.sciencedirect.com/science/article/pii/0377042787901257)

For the hierarchical clustering of point features, see: Conolly, J., & Lake, M. (2006). Geographic Information Systems in Archaeology. Cambridge: Cambridge University Press, 168-173.

The function also returns a list storing the following:
* $dist.matrix: distance matrix;
* $avr.silh.width.by.n.of.clusters: average silhouette width by number of clusters;
* $partition.silh.data: silhouette data for the selected partition;
* $dataset: the input dataset with two variables added ($feat_ID and $clust, the latter storing the cluster membership).

<br>

`kwPlot()`: function for visually displaying Kruskal-Wallis test's results. The function allows to perform Kruskal-Wallis test, and to display the test's results in a plot along with boxplots. The boxplots display the distribution of the values of the two samples, and jittered points represent the individual observations. At the bottom of the chart, the test statistics (H) is reported, along with the degrees of freedom and the associated p value. Setting the parameter 'posthoc' to TRUE, the Dunn's test is performed (with Bonferroni adjustment by default): a dot chart is returned, as well as a list of p-values (2-sided). In the dot chart, a RED line indicates the 0.05 threshold. The groups compared on a pairwise basis are indicated on the left-hand side of the chart.

<br>

`landfClass()`: function that allows to perform landform classification on the basis of the Topographic Position Index calculated from an input Digital Terrain Model (RasterLayer class). The TPI is the difference between the elevation of a given cell and the average elevation of the surrounding cells in a user defined moving window. For landform classification, the TPI is first standardized and then thresholded; to isolate certain classes, a slope raster (which is internally worked out) is also needed.
For details about the implemented classification, see: http://www.jennessent.com/downloads/tpi_documentation_online.pdf.
Two methods are available:
-the first (devised by Weiss) produces a 6-class landform classification comprising
* valley
* lower slope
* flat slope
* middle slope
* upper clope
* ridge
-the second (devised by Jennes) produces a 10-class classification comprising
* canyons, deeply incised streams
* midslope drainages, shallow valleys
* upland drainages, headwaters
* u-shaped valleys
* plains
* open slopes
* upper slopes, mesas
* local ridges, hills in valleys
* midslope ridges, small hills
* mountain tops, high ridges.
The second classification is based on two TPI that make use of two neighborhoods (moving windows) of different size: a s(mall) n(eighborhood) and a l(arge) n(eighborhood), defined by the parameters `sn` and `ln`. Besides rasters representing the different landform classes, the function optionally returns the TPI raster, either un- or standarized.

<br>

`logregr()`: function easy binary Logistic Regression and model diagnostics. The function allows to make it easy to perform binary Logistic Regression, and to graphically display the estimated coefficients and odds ratios. It also allows to visually check model's diagnostics such as outliers, leverage, and Cook's distance. The function may take a while (just matter of few seconds) to completed all the operations, and will eventually return the following charts:
* (1) Estimated coefficients, along with each coefficient's confidence interval; a reference line is set to 0. Each bar is given a color according to the associated p-value, and the key to the color scale is reported in the chart's legend.
* (2) Odds ratios and their confidence intervals.
* (3) A chart that is helpful in visually gauging the discriminatory power of the model: the predicted probability (x axis) are plotted against the dependent variable (y axis). If the model proves to have a high discriminatory power, the two stripes of points will tend to be well separated, i.e. the positive outcome of the dependent variable (points with color corresponding to 1) would tend to cluster around high values of the predicted probability, while the opposite will hold true for the negative outcome of the dependent variable (points with color corresponding to 0). In this case, the AUC (which is reported at the bottom of the chart) points to a low discriminatory power.
* (4) Model's standardized (Pearson's) residuals against the predicted probability; the size of the points is proportional to the Cook's distance, and problematic points are flagged by a label reporting their observation number if the following two conditions happen: residual value larger than 3 (in terms of absolute value) AND Cook's distance larger than 1. Recall that an observation is an outlier if it has a response value that is very different from the predicted value based on the model. But, being an outlier doesn't automatically imply that that observation has a negative effect on the model; for this reason, it is good to also check for the Cook's distance, which quantifies how influential is an observation on the model's estimates. Cook's distance should not be larger than 1.
* (5) Predicted probability plotted against the leverage value; dots represent observations, and their size is proportional to their leverage value, and their color is coded according to whether or not the leverage is above (lever. not ok) or below (lever. ok) the critical threshold. The latter is represented by a grey reference line, and is also reported at the bottom of the chart itself. An observation has high leverage if it has a particularly unusual combination of predictor values. Observations with high leverage are flagged with their observation number, making it easy to spot them within the dataset. Remember that values with high leverage and/or with high residual may be potential influencial points and may potentially negatively impact the regression. As for the leverage threshold, it is set at `3*(k+1)/N` (following Pituch-Stevens, *Applied Multivariate Statistics for the Social Science. Analyses with SAS and IBM's SPSS*, Routledge: New York 2016), where k is the number of predictors and N is the sample size.
* (6) Predicted probability against the Cook's distance.
* (7) Standardized (Pearson's) residuals against the leverage; points representing observations with positive or negative outcome of the dependent variable are given different colors. Further, points' size is proportional to the Cook's distance. Leverage threshold is indicated by a grey reference line, and the threshold value is also reported at the bottom of the chart. Observations are flagged with their observation number if their residual is larger than 3 (in terms of absolute value) OR if leverage is larger than the critical threshold OR if Cook's distance is larger than 1. This allows to easily check which observation turns out to be an outlier or a high-leverage data point or an influential point, or a combination of the three.
* (8) Chart that is almost the same as (7) except for the way in which observations are flagged. In fact, they are flagged if the residual is larger than 3 (again, in terms of absolute value) OR if the leverage is higher than the critical threshold AND if a Cook's distance larger than 1 plainly declares them as having a high influence on the model's estimates. Since an observation may be either an outlier or a high-leverage data point, or both, and yet not being influential, the chart allows to spot observations that have an undue influence on our model, regardless of them being either outliers or high-leverage data points, or both.
* (9) Observation numbers are plotted against the standardized (Pearson's) residuals, the leverage, and the Cook's distance. Points are labelled according to the rationales explained in the preceding points. By the way, the rationale is also explained at the bottom of each plots.

The function also returns a list storing two objects: one is named 'formula' and stores the formula used for the logistic regression; the other contains the model's results.

<br>

`modelvalid()`: allows to perform internal validation of a binary Logistic Regression model implementing most of the procedure described in
Arboretti Giancristofaro R, Salmaso L. "Model performance analysis and model validation in logistic regression". Statistica 2003(63): 375–396.

The procedure consists of the following steps:
* (1) the whole dataset is split into two random parts, a fitting (75 percent) and a validation (25 percent) portion;
* (2) the model is fitted on the fitting portion (i.e., its coefficients are computed considering only the observations in that portion) and its performance is evaluated on both the fitting and the validation portion, using AUC as performance measure;
* (3) the model's estimated coefficients, p-values, and the p-value of the Hosmer and Lemeshow test are stored;
* (4) steps 1-3 are repeated B times, eventually getting a fitting and validation distribution of the AUC values and of the HL test p-values, as well as a fitting distribution of the coefficients and of the associated p-values.

The AUC fitting distribution provides an estimate of the performance of the model in the population of all the theoretical fitting samples; the AUC validation distribution represents an estimate of the model’s performance on new and independent data.

The function returns:
* a chart with boxplots representing the fitting distribution of the estimated model's coefficients; coefficients' labels are flagged with an asterisk when the proportion of p-values smaller than 0.05 across the selected iterations is at least 95 percent;
* a chart with boxplots representing the fitting and the validation distribution of the AUC value across the selected iterations; for an example of the interpretation of the chart, see the aforementioned article, especially page 390-91;
* a chart of the levels of the dependent variable plotted against the predicted probabilities (if the model has a high discriminatory power, the two stripes of points will tend to be well separated, i.e. the positive outcome of the dependent variable will tend to cluster around high values of the predicted probability, while the opposite will hold true for the negative outcome of the dependent variable);
* a list containing:
  + $overall.model.significance: statistics related to the overall model p-value and to its distribution across the selected iterations;
  + $parameters.stability: statistics related to the stability of the estimated coefficients across the selected iterations;
  + $p.values.stability: statistics related to the stability of the estimated p-values across the selected iterations;
  + $AUCstatistics: statistics about the fitting and validation AUC distribution;
  + $Hosmer-Lemeshow statistics: statistics about the fitting and validation distribution of the HL test p-values.

As for the abovementioned statistics:
* full: statistic estimated on the full dataset;
* median: median of the statistic across the selected iterations;
* QRNG: interquartile range across the selected iterations;
* QRNGoverMedian: ratio between the QRNG and the median, expressed as percentage;
* min: minimum of the statistic across the selected iterations;
* max: maximum of the statistic across the selected iterations;
* percent_smaller_0.05: (only for $overall.model.significance, $p.values.stability, and $Hosmer-Lemeshow statistics): proportion of times in which the p-values are smaller than 0.05; please notice that for the overall model significance and for the p-values stability it is desirable that
the percentage is at least 95percent, whereas for the HL test p-values it is indeed desirable that the proportion is not larger than 5percent (in line with the interpetation of the test p-value which has to be NOT significant in order to hint at a good fit);
* significant (only for $p.values.stability): asterisk indicating that the p-values of the corresponding coefficient resulted smaller than 0.05 in at least 95percent of the iterations.

<br>

`monthlyWind()`: function that allows to download wind data from NOAA/NCEP Global Forecast System (GFS) Atmospheric Model colection, creating monthly averages. The extent of the study area can be specified: (1) by entering the geographic coordinates; (2) on the basis of an input raster dataset representing the study area itself (for instance, a Digital Elevation Model); (3) using a country code (for a list of country codes see http://kirste.userpage.fu-berlin.de/diverse/doc/ISO_3166.html). 

The function saves two .geotiff files in the computer's working directory, one representing the wind speed, the other the wind direction. In both cases, the values are the average of the wind speed and direction values in the study area across the days of the selected month, in the selected year. A plot is also returned in the R console.
The function returns a list containing the following data:
* `$windMonth`: stores the U and V components for each output grid cell (spatial resolution 0.5 degrees=50 Km);
* `$windMonthFit`: stores the wind speed and direction for each output grid cell.

The function builds upon the `wind.dl()` function from Javier Fernández-López's package `rWind`. The help provided by Dr Fernández-López in creating an earlier version of the `monthlyWind()` function is gratefully acknowledged.

<br>

`moveCost()`: provides the facility to calculate the accumulated cost of movement around a starting location and to optionally calculate least-cost paths toward one or multiple destinations. The function implements different cost estimations directly or inderectly related to human movement across the landscape.

The function takes as input a Digital Terrain Model (`RasterLayer` class) and a point feature (`SpatialPointsDataFrame` class), the latter representing
the starting location, i.e. the location from which the accumulated cost is calculated. If the parameter `destin` is fed with a datset representing destination location(s) (`SpatialPointsDataFrame` class), the function will also calculate least-cost path(s) plotted on the input DTM; the length of each path will be saved under the variable `length` stored in the 'LCPs' dataset (SpatialLines class) returned by the function. The red dot(s) representing the destination location(s) will be labelled with numeric values representing cost value at the location(s). The cost value will be also appended to the updated destination dataset returned by the function and storing a new variable named `cost`.

The function builds on functions out of the `gdistance` package, and by default uses a 16-directions movement in calculating the accumulated cost-surface.
The number of movements can be set by the user via the `moves` parameter. 

As noted in Nakoinz-Knitter (2016). "Modelling Human Behaviour in Landscapes". New York: Springer, p. 183, "gdistance works with conductivity rather than the more usual approach using costs, we need inverse cost functions". For this reason, in this function the cost is calculated using the inverse of the published cost functions.

The following cost functions are implemented:

* Tobler's hiking function (on-path): reshaped as follows

`((6 * exp(3.5 * abs(slope + 0.05))) * 0.278)^-1`

if we use speed, the final accumulated values will be 1/travel time (according to the 'gdistance' help documentation; see page 16 of this PDF: https://cran.r-project.org/web/packages/gdistance/vignettes/gdistance1.pdf); therefore, we use the reciprocal of speed to eventually get travel time/1. The Tobler's equation is multiplied by 0.278 (which is the ratio between 1000 -meters in 1 km- and 3600 -seconds in 1 hour-) to turn KmH into m/s,
and then reciprocated to turn m/s to s/m; before being reciprocating, we drop the minus before the 3.5. The same applies to the other Tobler's function-related equations listed below.


* Tobler's hiking function (off-path):

`(((6 * exp(3.5 * abs(slope + 0.05))) * 1.666667) * 0.278)^-1`

as per Tobler's indication, the off-path walking speed is reduced by 0.6; we use the reciprocal (1.666667) in our case since we are dealing with pace instead of speed.


* Márquez-Pérez et al.'s modified Tobler hiking function:

`((4.8 * exp(5.3 * abs((slope[adj] * 0.7) + 0.03))) * 0.278)^-1`

modified version as proposed by Joaquín Márquez-Pérez, Ismael Vallejo-Villalta & José I. Álvarez-Francoso (2017), "Estimated travel time for walking trails in natural areas", Geografisk Tidsskrift-Danish Journal of Geography, 117:1, 53-62, DOI: 10.1080/00167223.2017.1316212.


* Relative energetic expenditure cost function:

`(tan((atan(abs(slope) * 57.29578) * 0.0174532925)  / tan(1 * 0.0174532925)))^-1`

slope-based cost function expressing change in potential energy expenditure; in the above formula, `atan(abs(slope)) * 57.29578)` turns rise-over-run to degrees; multiplying by `0.0174532925` turns degrees to radians before calculating the tangent; the same applies to the degrees in the denominator. See Conolly, J., & Lake, M. (2006). Geographic Information Systems in Archaeology. Cambridge: Cambridge University Press, p. 220; see also Newhard, J. M. L., Levine, N. S., & Phebus, A. D. (2014). The development of integrated terrestrial and marine pathways in the Argo-Saronic region, Greece. Cartography and Geographic Information Science, 41(4), 379–390, with references to studies that use this function.


* Herzog's metabolic cost function in J/(kg*m):

`(1337.8 * slope^6 + 278.19 * slope^5 - 517.39 * slope^4 - 78.199 * slope^3 + 93.419 * slope^2 + 19.825 * slope + 1.64)^-1`

see Herzog, I. (2016). Potential and Limits of Optimal Path Analysis. In A. Bevan & M. Lake (Eds.), Computational Approaches to Archaeological Spaces (pp. 179–211). New York: Routledge.


* Wheeled-vehicle critical slope cost function:

`(1 + ((abs(slope)*100) / sl.crit))^-1`

where `sl.crit` (=critical slope, in percent) is "the transition where switchbacks become more effective than direct uphill or downhill paths" and typically is in the range 8-16; see Herzog, I. (2016). Potential and Limits of Optimal Path Analysis. In A. Bevan & M. Lake (Eds.), Computational Approaches to Archaeological Spaces (pp. 179–211). New York: Routledge.


* Pandolf et al.'s metabolic energy expenditure cost function (in watts):

`(1.5 * W + 2.0 * (W + L) * (L / W)^2 + N * (W + L) * (1.5 * V^2 + 0.35 * V * abs(slope*100)))^-1`

where `W` is the walker's body weight (Kg), `L` is the carried load (in Kg), `V` is the velocity in m/s, `N` is a coefficient representing ease of movement on the terrain. See Pandolf, K. B., Givoni, B., & Goldman, R. F. (1977). Predicting energy expenditure with loads while standing or walking very slowly. Journal of Applied Physiology, 43(4), 577–581. https://doi.org/10.1152/jappl.1977.43.4.577. This cost function is used, for instance, in this case study: Rademaker, K., Reid, D. A., & Bromley, G. R. M. (2012). Connecting the Dots: Least Cost Analysis, Paleogeography, and the Search for Paleoindian Sites in Southern Highland Peru. In D. A. White & S. L. Surface-Evans (Eds.), Least Cost Analysis of Social Landscapes. Archaeological Case Studies (pp. 32–45). University of Utah Press; see also Herzog, I. (2013). Least-cost Paths - Some Methodological Issues, Internet Archaeology 36 (http://intarch.ac.uk/journal/issue36/index.html) with references.


* Van Leusen's metabolic energy expenditure cost function (in watts):

`(1.5 * W + 2.0 * (W + L) * (L / W)^2 + N * (W + L) * (1.5 * V^2 + 0.35 * V * abs(slope*100 + 10)))^-1`

which modifies the Pandolf et al.'s equation; see Van Leusen, P. M. (2002). Pattern to process: methodological investigations into the formation and interpretation of spatial patterns in archaeological landscapes. University of Groningen. Note that, as per Herzog, I. (2013). Least-cost Paths - Some Methodological Issues, Internet Archaeology 36 (http://intarch.ac.uk/journal/issue36/index.html) and unlike Van Leusen (2002), in the above equation slope is expressed in percent and speed in m/s; also, in the last bit of the equantion, 10 replaces
the value of 6 used by Van Leusen.

When using the Tobler-related cost functions, the time unit can be selected by the user setting the `time` parameter to `h` (hour) or to `m` (minutes). In general, the user can also select which type of visualization the function has to produce; this is achieved setting the 'outp' parameter to either `r` (=raster) or to `c` (=contours). The former will produce a raster image with a colour scale and contour lines representing the accumulated cost surface; the latter parameter will only produce contour lines. The contour lines' interval is set using the parameter `breaks`; is not value is passed to the parameter, the interval will be set by default to 1/10 of the range of values of the accumulated cost surface.

The function returns a list storing:
* `$accumulated.cost.raster`: raster representing the accumualted cost (`RasterLayer` class);
* `$isolines`: contour lines derived from the accumulated cost surface (`SpatialLinesDataFrame` class);
* `$LCPs`: estimated least-cost paths (`SpatialLines` class);
* `$LCPs$length`: length of each least-cost path (units depend on the unit used in the input DTM);
* `$dest.loc.w.cost`: copy of the input destination location(s) dataset with a new variable (`cost`) added.

<br>

`mwPlot()`: function for visually displaying Mann-Whitney test's results. The function allows to perform Mann-Whitney test, and to display the test's results in a plot along with two boxplots. For information about the test, and on what it is actually testing, see for instance the interesting article by R M Conroy, *What hypotheses do "nonparametric" two-group tests actually test?*, in The Stata Journal 12 (2012): 1-9. The returned boxplots display the distribution of the values of the two samples, and jittered points represent the individual observations. At the bottom of the chart, a subtitle arranged on three lines reports relevant statistics:
* test statistic (namely, U) and the associated z and p value;
* Probability of Superiority value (which can be interpreted as an effect-size measure);
* another measure of effect size, namely r, whose thresholds are indicated in the last line of the plot's subtitle.

The function may also return a density plot (coupled with a rug plot at the botton of the same chart) that displays the distribution of the pairwise differences between the values of the two samples being compared. The median of this distribution (which is represented by a blue reference line in the same chart) corresponds to the Hodges-Lehmann estimator.

<br>

`NNa()`: futhe function allows to perform the Nearest Neighbor analysis of point patterns to formally test for the presence of a clustered, dispersed, or random spatial arrangement (second-order effect). It also allows to control for a first-order effect (i.e., influence of an underlaying numerical covariate) while performing the analysis. The covariate must be of RasterLayer class. Significance is assessed via a randomized approach. 

The function uses a randomized approach to test the significance of the Clark-Evans R statistic: the observed R value is set against the distribution of R values computed across B iterations (199 by default) in which a set of random points (with a sample size equal to the number of points of the input feature) is drawn and the statistic recomputed. The function produces a histogram of the randomized R values, with a black dot indicating the observed value and a hollow dot representing the average of the randomized R values. P-values (computed following Baddeley et al., "Spatial Point Patterns. Methodology and Applications with R", CRC Press 2016, p. 387), are reported at the bottom of the same chart. Two reference lines represent the two tails of the randomized distribution (left tail, indicating a significant clustered pattern; right tail, indicating a significant dispersed pattern).

The function also returns a list storing the following:
* $obs.NN.dist: observed NN distances;
* $obs.R: observed R value;
* $aver.rand.R: average randomized R;
* $p.value clustered: p-value for a clustered pattern;
* $p.value.dispersed: p-value for a dispersed pattern;
* $p.value.diff.from.random: p-value for a pattern different from random.

<br>

`outlier()`: function for univariate outliers detection. The function allows to perform univariate outliers detection using three different methods. These methods are those described in: Wilcox R R, *Fundamentals of Modern Statistical Methods: Substantially Improving Power and Accuracy*, Springer 2010 (2nd edition), pages 31-35.
Two of the three methods are robust, and are therefore less prone to the masking effect.
* (1) With the mean-based method, an observation is considered outlier if the absolute difference between that observation and the sample mean is more than 2 Standard Deviations away (in either direction) from the mean. In the plot returned by the function, the central reference line is indicating the mean value, while the other two are set at `mean-2*SD and mean+2*SD`.
* (2) The median-based method considers an observation as being outlier if the absolute difference between the observation and the sample median is larger than the Median Absolute Deviation divided by 0.6745. In this case, the central reference line is set at the median, while the other two are set at `median-2*MAD/0.6745` and `median+2*MAD/0.6745`.
* (3) The boxplot-based method considers an observation as being an outlier if it is either smaller than the 1st Quartile minus 1.5 times the InterQuartile Range, or larger than the 3rd Quartile minus 1.5 times the InterQuartile Range. In the plot, the central reference line is set at the median, while the other two are set at `1Q-1.5*IQR` and `3Q+1.5*IQR`.

The function also returns a list containing information about the choosen method, the mid-point, lower and upper boundaries where non-outlying observations are expected to fall, total number of outlying observations, and a dataframe listing the observations and indicating which is considered outlier. In the charts, the outlying observations are flagged with their ID number.

<br>

`perm.t.test()`: function for permutation-based t-test. The function allows to perform a permutation-based t-test to compare two independent groups. The test's results are graphically displayed within the returned chart. A permutation t-test proves useful when the assumption of 'regular' t-test are not met. In particular, when the two groups being compared show a very skewed distribution, and when the sample sizes are very unbalanced.
*"The permutation test is useful even if we plan to use the two-sample t test. Rather than relying on Normal quantile plots of the two samples and the central limit theorem, we can directly check the Normality of the sampling distribution by looking at the permutation distribution. Permutation tests provide a “gold standard” for assessing two-sample t tests. If the two P-values differ considerably, it usually indicates that the conditions for the two-sample t don’t hold for these data. Because permutation tests give accurate P-values even when the sampling distribution is skewed, they are often used when accuracy is very important."* (Moore, McCabe, Craig, *Introduction to the Practice of Statistics*, New York: W. H. Freeman and Company, 2009).
The chart returned by the function diplays the distribution of the permuted mean difference between the two samples; a dashed line indicates the observed mean difference. A rug plot at the bottom of the density curve indicates the individual permuted mean differences. Under the chart, a number of information are displayed. In particular, the observed mean difference, the number of permutations used, and the permuted p-value are reported. In the last row, the result of the regular t-test (both assuming and not assuming equal variances) is reported to allow users to compare the outcome of these different versions of the test.

<br>

`plotJenks()`: function for plotting univariate classification using Jenks' natural break method. The function allows to break a dataset down into a user-defined number of breaks and to nicely plot the results, adding a number of other relevant information. Implementing the Jenks' natural breaks method, it allows to find the best arrangement of values into different classes.
The function produces a chart in which the values of the input variable are arranged on the x-axis in ascending order, while the index of the individual observations is reported on the y-axis. Vertical dotted red lines correspond to the optimal break-points which best divide the input variable into the selected classes. The break-points (and their values) are reported in the upper part of the chart, onto the corresponding break lines. Also, the chart's subtitle reports the Goodness of Fit value relative to the selected partition, and the partition which correspond to the maximum GoF value.
The function also returns a list containing the following:
* `$info`: information about whether or not the method created non-unique breaks;
* `$classif`: created classes and number of observations falling in each class;
* `$classif$brks`: classes' break-points;
* `$breaks$max.GoF`: number of classes at which the maximum GoF is achieved;
* `$class.data`: dataframe storing the values and the class in which each value actually falls into.

<br>

`pointsCovarModel()`: the function is a wrapper for a number of functions out of the extremely useful `spatstat` package (specifically, `ppm()`, `cdf.test()`, `auc()`, `roc()`, `effectfun()`). It allows to test if there is a significant dependence of the input point pattern on a underlying spatial numeric covariate (first-order effect). The function takes as input three datasets: a point patter (`SpatialPointsDataFrame` class), a covariate layer (of `RasterLayer` class), and a polygon feature (`SpatialPolygonsDataFrame` class) representing the study area and exactly matching the extent of the covariate layer. If the latter is not provided, it is internally worked out from the covariate raster and may make the whole function take a while to complete.

The function fits a inhomogeneous Poisson point process (Alternative Model-H1) with intensity of the point pattern as a loglinear function of the underlaying numerical covariate (see Baddeley et al., "Spatial Point Patterns. Methodology and Applications with R", CRC Press 2016, 307-309). Also, the function fits a homogeneous Poisson point model (Null Model-H0, equivalent to Complete Spatial Randomness: Baddeley et al., "Spatial Point Patterns. Methodology and Applications with R", CRC Press 2016, 305-306), that is used as comparison for the inhomogeneous point process model in a Likelihood Ratio test (Baddeley et al., "Spatial Point Patterns. Methodology and Applications with R", CRC Press 2016, 334-335). A significant result, i.e. a low p-value, suggests rejecting the Null Hypothesis of CSR in favour of the Alternative Hypothesis of a Poisson point process affected by a covariate effect (i.e., inhomogeneous intensity due to the influence of the covariate) (Baddeley et al., "Spatial Point Patterns. Methodology and Applications with R", CRC Press 2016, 305).

The function returns a 4 plots, which can be arranged in just one visualization setting the parameter `oneplot` to TRUE:
* plot of the point pattern along with the underlaying covariate raster;
* plot of the fitted intensity against the spatial covariate (Baddeley et al., "Spatial Point Patterns. Methodology and Applications with R", CRC Press 2016, 308);
* plot of the cumulative distribution of the covariate at the data points against the cumulative distribution of the covariate at all the spatial location within the study area (rationale: Baddeley et al., "Spatial Point Patterns. Methodology and Applications with R", CRC Press 2016, 184-185);
* plot of the ROC curve, which help assessing the strenght of the dependence on the covariate (Baddeley et al., "Spatial Point Patterns. Methodology and Applications with R", CRC Press 2016, 187-188).

A list is also returned, containing what follows:
* $H0-model: info and relevant statistics regarding the Null Model;
* $H1-model: info and relevant statistics regarding the Alternative Model;
* $Model comparison (LRT): results of the Likelihood Ratio test;
* $AIC-H0: AIC of the Null Model;
* $AIC-H1: AIC of the Atlernative Model;
* $KS test: information regarding the cumulative distribution comparison via Kolmogorov-Smirnov test;
* AUC: the AUC statistics.

<br>

`pointsInPolygons()`: the function allows to test:

* `scenario a`: if there is a significant spatial association between a set of points and a set of polygons, in terms of points falling within the polygons. In other words, it aims at testing whether a set of points falls inside a set of polygons more often than would be expected by chance. The basic assumption is that the polygons are completely contained within the study plot. If the shapefile (of polygon type) representing the study plot is not provided, the calculations use the bounding polygon based on the union the convex hulls of the point and of the polygon feature.
* `scenario b`: if the distribution of points within a set of polygons totally covering the study area can be considered random, or if the observed points count for each polygon is larger or smaller than expected. P values are also reported.

The computations relative to scenario `a` are based on the `dbinom()` and `pbinom()` functions. The probability of observed count within polygons is `dbinom(x, size=n.of.points, prob=p)`, where `x` is the observed number of points within polygons, `n.of.points` is the total number of points, and `p` is the probability that a single point will be found within a polygon, which is equal to the ratio between the area of the polygons and the total area of the study plot. The probability that x or fewer points will be found within the polygons is `pbinom(x, size=n.of.points, prob=p)`.

The calculations relative to the scenario `b` are again based on the binomial distribution: the probability of the observed counts is `dbinom(x, size=n.of.points, prob=p)`, where `x` is the observed number of points within a given polygon, `n.of.points` is the total number of points, and `p` is equal to the size of each polygon relative to sum of the polygons' area. The probability that x or fewer points will be found within a given polygon is `pbinom(x, size=n.of.points, prob=p)`.

For scenario `a` the function produces a plot of the points and polygons (plus the convex hull of the study area), and relevant information are reported at the bottom of the chart itself.
A list is also returned, containing what follows:
* `$Polygons' area`;
* `$Study area's area`;
* `$Total # of points`;
* `$Observed # of points in polygons`;
* `$Expected # of points in polygons`;
* `$Exact probability of observed count within polygons`;
* `$Probability of <= observed count within polygons`;
* `$Probability of >= observed count within polygons`.

For scenario `b` the function returns a plot showing the polygons plus the dots; in each polygon the observed and expected counts are reported, and the p-value of the observed count is indicated. A matrix is also returned, containing what follows:
* `polygons' area`;
* `%area` (size of each polygon relative to sum of the polygons' area; it corresponds to the probability (p) fed into the binomial distribution function);
* `observed number of points`;
* `expected number of points`;
* `probability of observed counts`;
* `probability of observed counts <= than expected`;
* `probability of observed counts >= than expected`.

<br>

`pointsToPointsTess()`: the function can be considered as a special case of the `scenario b` tested by the `pointsInPolygons()` function provided by this same package, with the exception that in this case the polygons are not entered by the use but are internally created by the function aorund the to-feature. The question this function may allow to address is: do the points belonging to a feature dataset tend to occur close to any of the points in another feature dataset than expected if the points would be randomly scattered across the study area? To help addressing this question, the function creates Thiessen polygons around the input `to.feature` and then runs the `pointsInPolygons()` function using its `scenario b`. For further details, see the help documentation of the `pointsInPolygons()` function.

<br>

`ppdPlot()`: function for plotting Posterior Probability Densities for Bayesian modeled 14C dates/parameters. The function allows plot Posterior Probability Densities with a nice outlook thanks to `ggplot2`. It takes as input a dataframe that must be organized as follows (it is rather easy to do that once the data have been exported from OxCal):
* calendar dates (first column to the left);
* posterior probabilities (second column);
* grouping variables (third column), which could contain the names of the events of interest (e.g., phase 1 start, phase 1 end, phase 2 start, phase 2 end, etc).

<br>

`prob.phases.relat()`: function to calculate the Posterior Probability for different chronological relations between two Bayesian radiocarbon phases. The function allows to calculate the posterior probability for different chronological relations between two phases defined via Bayesian radiocarbon modeling. For the results to make sense, the phases have to be defined as independent if one wishes to assess what is the posterior probability for different relative chronological relations between them.
The rationale for this approach is made clear in an article by Buck et al 1992 (https://doi.org/10.1016/0305-4403(92)90025-X), and it runs as follows: *if we do not make any assumption about the relationship between the phases, can we test how likely they are to be in any given order*?

The function takes as input the table provided by the 'OxCal' program as result of the 'Order' query. Once the table as been saved from 'OxCal' in .csv format, you have to feed it in R. A .csv file can be imported into R using (for instance): `mydata <- read.table(file.choose(), header=TRUE, sep=",", dec=".", as.is=T)`
Make sure to insert the phases' parameters (i.e., the starting and ending boundaries of the two phases) in the OxCal's Order query in the following order: StartA, EndA, StartB, EndB. That is, first the start and end of your first phase, then the start and end of the second one. You can give any name to your phases, as long as the order is like the one described.

Given two phases A and B, the function allows to calculate the posterior probability for:
* A being within B
* B being within A
* A starting earlier but overlapping with B
* B starting earlier but overlapping with A
* A being entirely before B
* B being entirely before A
* sA being within B
* eA being within B
* sB being within A
* eB being within A
where 's' and 'e' refer to the starting and ending boundaries of a phase. 

The function will return a table and a dot plot.
Thanks are due to Dr. Andrew Millard (Durham University) for the help provided in working out the operations on probabilities.

<br>

`refNNa()`: allows to perform the refined Nearest Neighbor analysis of point patterns by plotting the cumulative Nearest Neighbour distance, along with an acceptance interval (with significance level equal to 0.05; sensu Baddeley et al., "Spatial Point Patterns. Methodology and Applications with R", CRC Press 2016, 208) based on B (set to 200 by default) realizations of a Complete Spatial Random process. The function also allows to control for a first-order effect (i.e., influence of an underlaying numerical covariate) while performing the analysis. The covariate must be of RasterLayer class.

The function uses a randomized approach to build the mentioned acceptance interval, whereby cumulative distributions of average NN distances of random points are computed across B iterations. In each iteration, a set of random points (with sample size equal to the number of points of the input feature) is drawn.

<br>

`resc.val()`: function that allows to rescale the values of a dataset between a minimum and a maximum that are specified by the user. In doing that, it allows to preserve the shape of the distribution of the original data. The function produces two density charts representing the distribution of the original and of the rescaled dataset. It also returns a dataframe storing the original and rescaled values.

<br>

`robustBAplot()`: function to plot a robust version of the Bland-Altman plot. It returns a chart based on robust (i.e. resistant to outlying values) measures of central tendency and variability: median and Median Absolute Deviation (MAD) (Wilcox R R. 2001. *Fundamentals of modern statistical methods: Substantially improving power and accuracy*. New York: Springer) instead of mean and standard deviation.
The x-axis displays the median of the two variables being compared, while the y-axis represents their difference. A solid horizontal line represents the bias, i.e. the median of the differences reported on the y-axis. Two dashed horizontal lines represent the region in which 95percent of the observations are expected to lie; they are set at the median plus or minus `z*(MAD/0.6745)`.

<br>

`vislim()`: function for computing the limit of visibility of an object given its height. It plots the angular size of an object (in degrees) against the distance from the observer, and computes at which distance from the observer the angular size of the object hits the limit of human visual acuity (0.01667 degrees). 

The function returns:
* a plot displaying the decay in angular size as function of the object's distance from the observer; a black dot represents the distance at which the angular size hits the limit of human visual acuity;
* the value (in km) of the visibility limit.

<br>

## History
`version 0.19`: 
improvements and typos fixes to the help documentation; 
minor improvement to the title of chart returned by the `Aindex()` function; 
minor improvements to the annotations present in some of the plots returned by the `logregr()` function; 
Hosmer-Lemeshow test added to the routines implemented in the `modelvalid()` function; 
minor modification to the title and subtitle of the plot returned by the `distRandCum()` and `refNNa()` functions, where the term "95% confidence envelope" has been changed to "acceptance interval" to comply with the distinction stressed in Baddeley et al., "Spatial Point Patterns. Methodology and Applications with R", CRC Press 2016, page 208;
`moveCost()` function added;
`volc`, `volc.loc`, and `destin_loc` datasets added;
`gdistance` package added to the list of dependencies.

`version 0.18`: 
improvements and typos fixes to the help documentation; under-the-hood performance improvements;improvements to the `distRandCum()` function's chart; fixes aimed at restoring the default graphical device's settings after multiple plots are produced.

`version 0.17`: 
improvements and typos fixes to the help documentation; `distRandCum()` function added.

`version 0.16`: 
`rgdal` package added to the list of dependencies; substantial improvements to the `modelvalid()` function.

`version 0.15`: 
improvements and typos fixes to the help documentation. 
`featClust()` function added. 
Improvements to the `NNa()` function, which now provides the observed and randomized calculation of the Clark and Evans's R statistic, with significance computed via permutation procedure; the observed value and the mean value of the randomized Nearest Neighbor distance are provided within the list returned by the function; 2-tailed p value is not also reported, along with the lower (significance for clustered pattern) and upper (significance for dispersed pattern) tailed p values; histogram of the frequency distribution of the randomized R values is returned instead of a density chart. 
Improvements to the significance calculation in the `distRandSign()` function, which now returns a frequency distribution histogram instead of a density chart. Also, if the analysis takes into account two point patterns, the function now performs an additional test based on a permuted routine along the one provided in earlier versions, which was (and is still) based on a randomized procedure.
Improvements to the `refNNa()` function, which now can control for the effect of a numerical covariate (controlling for a 1st order effect). 
Minor improvements to point pattern plot optionally returned by the `Aindex()` function.

`version 0.14`: 
improvements and typos fixes to the help documentation; `Aindex()` function added.

`version 0.13`: 
improvements and typos fixes to the help documentation; improvements to the permuted p-value calculation in the `chiperm()` function; `landfClass()` and `resc.val()` functions added.

`version 0.12`: 
improvements to the output charts of the `chiperm()`, `distRandSign()`, `NNa()`, and `perm.t.test()` functions; 
`distCovarModel()`, `pointsCovarModel()`, `refNNa()`, and `vislim()` functions added; 
improvement to the randomized p-values calculation in the `chiperm()`, `distRandSign()`, `NNa()`, and `perm.t.test()` functions; 
`NNa()` modified to provide the facility to controll for the effect of a spatial covariate; 
`malta_dtm_40`, `malta_polyg`, `Massachusetts`, `Starbucks`, and `popdensity` datasets added;  
improvements and typos fixes to the help documentation.

`version 0.11`: 
improvements and typos fixes to the help documentation; improvements to the output charts of the `chiperm()` function.

`version 0.10`: 
`distRandSign()` and `pointsInPolygons()` now accept a shapefile (of polygon type) representing the study plot; dot (representing the mean value) added to the density chart returned by the `distRandSign()` function. The `NNa()` function has been added. Progression bar added to a number of functions. Improvements to the functions' help documentation.

`version 0.9`: 
improvements to the help documentation of a numbe of functions; minor inconsistencies fixed.

`version 0.8`: 
improvements to the `pointsToPointsTess()` function; `pumps` and `deaths` datasets added.

`version 0.7`: 
the `pointsInPolygons()` functions has been modified; the function `pointsToPointsTess()` and the `events`, `locations`, and `thiessenpolyg` datasets have been added.

`version 0.6`: 
the `distRandSign()` and `pointsInPolygons()` functions have been added, along with the `springs`, `faults`, `points`, and `polygons` datasets.

`version 0.5`: 
the facility to download wind data on the basis of country code has been added to the `monthlyWind()` function.

`version 0.4.1`: 
fix to `description` file.

`version 0.4`: 
fix to the help documentation of the `aucadj()` function; `monthlyWind()` function added.

`version 0.3`: 
typos correction in the README file, one conflicting keyword removed from the `modelvalid()` function.

`version 0.2`: 
bug fixes, improvement to the help documentation, improvements to the README file.

`version 0.1`: 
first release.

<br>

## Installation
To install the package in R, just follow the few steps listed below:

1) install the `devtools` package:  
```r
install.packages("devtools", dependencies=TRUE)
```
2) load that package: 
```r
library(devtools)
```
3) download the `GmAMisc` package from GitHub via the `devtools`'s command: 
```r
install_github("gianmarcoalberti/GmAMisc@v0.19")
```
4) load the package: 
```r
library(GmAMisc)
```
5) enjoy!
