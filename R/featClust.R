#' R function for points clustering on the basis of planar distance
#'
#' The function provides the facility to cluster the features of the input dataset on the basis of either their (projected) coordinates (for points; SpatialPointsDataFrame class)
#' or of their area (for polygons; SpatialPolygonsDataFrame class).
#'
#' The function internally calculates a distance matrix (based on the Euclidean Distance) on the basis of the points' coordinates or polygons' area.
#' A dendrogram is produced which depicts the hierarchical clustering based (by default) on the Ward's agglomeration method; rectangles identify the selected cluster partition.
#' Besides the dendrogram, a silhouette plot is produced, which allows to measure how 'good' is the selected cluster solution.\cr
#'
#' As for the latter, if the parameter 'part' is left empty (default), an optimal cluster solution is obtained.
#' The optimal partition is selected via an iterative procedure which locates at which cluster solution the highest average silhouette width is achieved.
#' If a user-defined partition is needed, the user can input the desired number of clusters using the parameter 'part'.
#' In either case, an additional plot is returned besides the cluster dendrogram and the silhouette plot; it displays a scatterplot in which the cluster solution (x-axis)
#' is plotted against the average silhouette width (y-axis). A black dot represent the partition selected either by the iterative procedure or by the user.\cr
#'
#' Notice that in the silhouette plot, the labels on the left-hand side of the chart show the point ID number and the cluster to which each point is closer.\cr
#'
#' Also, the function returns a plot showing the input point dataset, with points colored by cluster membership. Two new variables are added to the
#' point shapefile's dataframe, storing a point ID number and the corresponding cluster membership.\cr
#'
#' The silhouette plot is obtained from the silhouette() function out from the 'cluster' package (https://cran.r-project.org/web/packages/cluster/index.html).\cr
#' For a detailed description of the silhouette plot, its rationale, and its interpretation, see:\cr
#' Rousseeuw P J. 1987. "Silhouettes: A graphical aid to the interpretation and validation of cluster analysis", Journal of Computational and Applied Mathematics 20, 53-65 (http://www.sciencedirect.com/science/article/pii/0377042787901257)
#'
#' The function also returns a list storing the following:\cr
#'
#' -$dist.matrix: distance matrix;\cr
#' -$avr.silh.width.by.n.of.clusters: average silhouette width by number of clusters;\cr
#' -$partition.silh.data: silhouette data for the selected partition;\cr
#' -$dataset: the input dataset with two variables added ($feat_ID and $clust, the latter storing the cluster membership).
#'
#' @param x: dataset whose feature are to be clustered; either points (SpatialPointsDataFrame class) or polygons (SpatialPolygonsDataFrame class).
#' @param aggl.meth: agglomeration method ("ward.D2" by default).
#' @param part: desired number of clusters; if NULL (default), an optimal partition is calculated.
#' @param oneplot: TRUE (default) or FALSE if the user wants or does not want the plots to be visualized in a single window.
#' @param cex.dndr.lab: set the size of the labels used in the dendrogram.
#' @param cex.sil.lab: set the size of the labels used in the silhouette plot.
#' @keywords clustering
#' @export
#' @examples
#' data(springs)
#' res <- featClust(springs) #perform the analysis and automatically select an optimal partition
#' res <- featClust(springs, part=3) #as above, but selecting a 3-cluster partition
#'
featClust <- function(x, aggl.meth = "ward.D2", part=NULL, oneplot=TRUE, cex.dndr.lab = 0.85, cex.sil.lab = 0.75) {

  #if the input dataset is a SpatialPolygonsDataframe, create a dataframe storing polygons area;
  #otherwise store points' coordinates
  ifelse(class(x)[1]=="SpatialPolygonsDataFrame",
         df <- as.data.frame(gArea(x, byid=TRUE)),
         df <- as.data.frame(coordinates(x)))

  #store the number of features
  nx <- length(x)

  #calculate the Euclidean distance
  d <- dist(df)

  #calculate the hier. clustering
  fit <- hclust(d, method = aggl.meth)

  #calculate the max number of clusters to be later used in the loop which iterates the calculation of the average silhouette width
  max.ncl <- nx - 1

  #create a slot to store the values of the average silhouette value at increasing numbers of cluster solutions
  #this will be used inside the loop
  sil.width.val <- numeric(max.ncl - 1)

  #create a slot to store the increasing number of clusters whose silhouette is iteratively computed
  sil.width.step <- c(2:max.ncl)

  #calculate the average silhouette width at each increasing cluster solution
  for (i in 2:max.ncl) {
    counter <- i - 1
    clust <- silhouette(cutree(fit, k = i), d)
    sil.width.val[counter] <- mean(clust[, 3])
  }

  #create a dataframe that stores the number of cluster solutions and the corresponding average silhouette values
  sil.res <- as.data.frame(cbind(sil.width.step, sil.width.val))

  #if the user does not enter the desired partition, the latter is equal to the optimal one, ie the one with the largest average silhouette;
  #otherwise use the user-defined partition
  ifelse(is.null(part)==TRUE,
         select.clst.num <- sil.res$sil.width.step[sil.res$sil.width.val == max(sil.res$sil.width.val)],
         select.clst.num <- part)

  #create the final silhouette data using the partition established at the preceding step
  final.sil.data <- silhouette(cutree(fit, k = select.clst.num), d)

  #add the cluster membership to the points coordinates dataframe
  df$clust <- assignCluster(df, df, cutree(fit, k = select.clst.num))

  #add a variable to the input shapefile to store an ID for each feature
  x$feat_ID <- 1:nx

  #add a variable to the input shapefile to store the cluster membership taken from the df dataframe
  x$clust <- df$clust

  #set the layout of the plot visualization according to whether or not the parameter 'oneplot' is set to TRUE
  if(oneplot==TRUE){
    m <- rbind(c(1,2), c(3,4))
    layout(m)
  }

  #plot the points with colour according to the cluster membership
  plot(x,
       pch=20,
       col=x$clust,
       main="Plot of the features colored by cluster membership",
       cex.main=0.90,
       axes=TRUE)

  #plot the dendrogram
  plot(fit, main = "Clusters Dendrogram",
       sub = paste0("\nAgglomeration method: ", aggl.meth),
       xlab = "",
       cex = cex.dndr.lab,
       cex.main = 0.9,
       cex.sub = 0.75)

  #add the rectangles representing the selected partitions
  solution <- rect.hclust(fit, k = select.clst.num, border="black")

  #modify the row names of the final silhouette dataframe to also store the cluster to which each feature is closest
  rownames(final.sil.data) <- paste(1:nx, final.sil.data[, 2], sep = "_")

  #plot the silhouette chart
  plot(final.sil.data,
       cex.names = cex.sil.lab,
       max.strlen = 30,
       nmax.lab = nx + 1,
       main = "Silhouette plot")

  #add a line representing the average silhouette width
  abline(v = mean(final.sil.data[, 3]), lty = 2)

  #plot the average silhouette value at each cluster solution
  plot(sil.res, xlab = "number of clusters",
       ylab = "silhouette width",
       ylim = c(0, 1),
       xaxt = "n",
       type = "b",
       main = "Silhouette width vs. number of clusters",
       sub = paste0("values on the y-axis represent the average of the silhouettes' width distribution at each cluster solution"),
       cex.main=0.90,
       cex.sub = 0.75)

  axis(1, at = 0:max.ncl, cex.axis = 0.7)

  #draw a black dot corresponding to the selected partition
  points(x=select.clst.num,
         y=sil.res[select.clst.num-1,2],
         pch=20)

  #draw the label selecting the row from the sil.res dataframe corresponding to the selected partition
  text(x = select.clst.num,
       y = sil.res[select.clst.num-1,2],
       labels= round(sil.res[select.clst.num-1,2],3),
       cex = 0.65,
       pos = 3,
       offset = 1.2,
       srt = 90)

  return <- list("dist.matrix"=d,
                 "avr.silh.width.by.n.of.clusters"=sil.res,
                 "partition.silh.data"=final.sil.data,
                 "dataset"=x)
}
