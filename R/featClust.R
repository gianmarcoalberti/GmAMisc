#' R function for features clustering on the basis of distances/area
#'
#' The function provides the facility to cluster the features of the input dataset on the basis of either their (projected) coordinates (for points; SpatialPointsDataFrame class)
#' or of their area (for polygons; SpatialPolygonsDataFrame class). If a target feature dataset (to.feat) is provided, the clustering will be based
#' on the distance of the x feature to the nearest to.feature. When a to.feature is specified, the x feature (i.e., the feature that the user wants to cluster)
#' can be either a point (SpatialPointsDataFrame class), or a polyline (SpatialLinesDataFrame class), or a polygon (SpatialPolygonsDataFrame class) feature.
#' Notice that if all the x features overlap with all the to.feature, all the minimum distances will be 0, and the function will trow an error.\cr
#'
#' If the to.feature is not provided, the function internally calculates a distance matrix (based on the Euclidean Distance) on the basis of the points' coordinates or polygons' area.
#' If the to.feature is provided, the distance matrix will be based on the distance of the x feature to the nearest to.feature.
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
#' Also, the function returns a plot showing the input dataset, with features colored by cluster membership. Two new variables are added to the
#' shapefile's dataframe, storing a point ID number and the corresponding cluster membership.\cr
#'
#' The silhouette plot is obtained from the 'silhouette()' function out from the 'cluster' package (https://cran.r-project.org/web/packages/cluster/index.html).\cr
#' For a detailed description of the silhouette plot, its rationale, and its interpretation, see:\cr
#' Rousseeuw P J. 1987. "Silhouettes: A graphical aid to the interpretation and validation of cluster analysis", Journal of Computational and Applied Mathematics 20, 53-65 (http://www.sciencedirect.com/science/article/pii/0377042787901257)
#'
#' For the hierarchical clustering of features, see: Conolly, J., & Lake, M. (2006). Geographic Information Systems in Archaeology. Cambridge: Cambridge University Press, 168-173.\cr.
#'
#' The function also returns a list storing the following:\cr
#'
#' -$dist.matrix: distance matrix;\cr
#' -$avr.silh.width.by.n.of.clusters: average silhouette width by number of clusters;\cr
#' -$partition.silh.data: silhouette data for the selected partition;\cr
#' -$coord.or.area.or.min.dist.by.clust: coordinates, area, or distance to the nearest to.feat coupled with cluster membership;\cr
#' -$dist.stats.by.cluster: by-cluster summary statistics of the x feature distance to the nearest to.feature;\cr
#' -$dataset: the input dataset with two variables added ($feat_ID and $clust, the latter storing the cluster membership).
#'
#' @param x: dataset whose feature are to be clustered; either points (SpatialPointsDataFrame class) or polygons (SpatialPolygonsDataFrame class);
#' if the to.feat is specified, x can also be a polylines feature (SpatialLinesDataFrame class).
#' @param to.feat: dataset (NULL by default) representing the feature the distance toward which is used as basis for clustering x;
#' either points (SpatialPointsDataFrame class), polygons (SpatialPolygonsDataFrame class), or polylines (SpatialLinesDataFrame).
#' @param aggl.meth: agglomeration method ("ward.D2" by default).
#' @param part: desired number of clusters; if NULL (default), an optimal partition is calculated (see Details).
#' @param showID: TRUE (default) or FALSE if the used wants or does not want the ID of the clustered features to be displayed in the plot where
#' the features are colored by cluster membership.
#' @param oneplot: TRUE (default) or FALSE if the user wants or does not want the plots to be visualized in a single window.
#' @param cex.dndr.lab: set the size of the labels used in the dendrogram.
#' @param cex.sil.lab: set the size of the labels used in the silhouette plot.
#' @param cex.feat.lab: set the size of the labels used (if 'showID' is set to TRUE) to show the clustered features' IDs.
#' @param col.feat.lab: set the color of the clustered features' IDs ('black' by default).
#' @param export: TRUE or FALSE (default) if the user wants or does not want the clustered input dataset to be exported;
#' if TRUE, the input dataset with a new variable indicating the cluster membership will be exported as a shapefile.
#' @keywords clustering
#' @export
#' @examples
#' data(springs)
#' res <- featClust(springs) #perform the analysis and automatically select an optimal partition
#'
#' res <- featClust(springs, part=3) #as above, but selecting a 3-cluster partition
#'
#' res <- featClust(springs, faults) #cluster springs on the basis of their distance to the nearest geological fault
#'
#' res <- featClust(polygons, springs) #cluster polygonal areas on the basis of their distance to the nearest spring
#'
#' res<- featClust(destin_loc, volc.loc) #cluster location points on the basis of their distance from a single location point
#'
#' res <- featClust(points, polygons) #cluster points on the basis of their distance to the nearest polygon
#'
featClust <- function(x, to.feat=NULL, aggl.meth = "ward.D2", part=NULL, showID=TRUE, oneplot=TRUE, cex.dndr.lab = 0.85, cex.sil.lab = 0.75, cex.feat.lab=0.65, col.feat.lab="black", export=FALSE) {

  #if there is no to.feature, if the input dataset is a SpatialPolygonsDataframe, create a dataframe storing polygons area;
  #otherwise store points' coordinates; this is skipped if there actually is a to.feature
  if(is.null(to.feat)==TRUE){
    ifelse(class(x)[1]=="SpatialPolygonsDataFrame",
           df <- data.frame(area=gArea(x, byid=TRUE)),
           df <- as.data.frame(coordinates(x)))
  }

  #if there is a to.feature...
  if(is.null(to.feat)==FALSE){
    #calculate all the pairwise distances between x and to.feat
    dist.matrix <- gDistance(x, to.feat, byid=TRUE)
    #then calculate the minimum distances of x toward to.feat by appling the 'min' function column-wise
    df <- as.data.frame(apply(dist.matrix, 2, FUN=min))
    #give name to the df column storing the minimum distances
    colnames(df)[1] <- "minim.dist."
  }

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

  #if there is a to.feat, calculate distance summary statistics by cluster
  if(is.null(to.feat)==FALSE){
    dist.summary <- tapply(df$minim.dist., df$clust, summary)
  } else {
    dist.summary <- "NA"
  }

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
  #if there is not to.feat...
  if(is.null(to.feat)==TRUE){
    #plot the x feature
    plot(x,
         pch=20,
         col=x$clust,
         main="Plot of the features colored by cluster membership",
         cex.main=0.90,
         axes=TRUE)

    #show the ID of the feature if showID is equal to TRUE and the feature is NOT a polyline
    if(showID==TRUE & class(x)[1]!="SpatialLinesDataFrame"){
    text(coordinates(x),
         labels=x$feat_ID,
         pos = 4,
         cex=cex.feat.lab,
         col=col.feat.lab)
    }
  }

  #if there is a to.feat...
  if(is.null(to.feat)==FALSE){
    #calculate the extent of the x and to.feature
    area.x <- (extent(x)[2] - extent(x)[1]) * (extent(x)[4] - extent(x)[3])
    area.to.feat <- (extent(to.feat)[2] - extent(to.feat)[1]) * (extent(to.feat)[4] - extent(to.feat)[3])

    #give precedence to plotting according to which feature has the larger extent
    if(area.x > area.to.feat){
      plot(x,
           pch=20,
           col=x$clust,
           main="Plot of the features colored by cluster membership",
           cex.main=0.90,
           axes=TRUE)

      #show the ID of the feature if showID is equal to TRUE and the feature is NOT a polyline
      if(showID==TRUE & class(x)[1]!="SpatialLinesDataFrame"){
        text(coordinates(x),
             labels=x$feat_ID,
             pos = 4,
             cex=cex.feat.lab,
             col=col.feat.lab)
      }

      plot(to.feat, add=TRUE)

    } else {

      plot(to.feat,
           main="Plot of the features colored by cluster membership",
           cex.main=0.90,
           axes=TRUE)
      plot(x,
           pch=20,
           col=x$clust,
           add=TRUE)

      #show the ID of the feature if showID is equal to TRUE and the feature is NOT a polyline
      if(showID==TRUE & class(x)[1]!="SpatialLinesDataFrame"){
        text(coordinates(x),
             labels=x$feat_ID,
             pos = 4,
             cex=cex.feat.lab,
             col=col.feat.lab)
      }
    }
  }

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
       ylab = "average silhouette width",
       ylim = c(0, 1),
       xaxt = "n",
       type = "b",
       main = "Average silhouette width vs. number of clusters",
       sub = paste0("values on the y-axis represent the average silhouette width distribution at each cluster solution"),
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

  if(export==TRUE){
    writeOGR(x, ".", "clustered_dataset", driver="ESRI Shapefile")
  }

  results <- list("dist.matrix"=d,
                  "avr.silh.width.by.n.of.clusters"=sil.res,
                  "partition.silh.data"=final.sil.data,
                  "coord.or.area.or.min.dist.by.clust"=df,
                  "dist.stats.by.cluster"=dist.summary,
                  "dataset"=x)

  # restore the original graphical device's settings
  par(mfrow = c(1,1))

  return(results)
}
