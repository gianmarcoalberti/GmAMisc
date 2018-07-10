#' R function for Brainerd-Robinson similarity coefficient (and optional clustering)
#'
#' The function allows to calculate the Brainerd-Robinson similarity coefficient, taking as input a cross-tabulation (dataframe), and to optionally
#' perform an agglomerative hierarchical clustering.\cr
#'
#' The function returns produces a correlation matrix in tabular form and a heat-map representing, in a graphical form, the aforementioned correlation matrix.\cr
#'
#' In the heat-map (which is built using the 'corrplot' package), the size and the color of the squares are proportional to the Brainerd-Robinson coefficients, which are also reported by numbers.\cr
#'
#' In order to "penalize" BR similarity coefficient(s) arising from assemblages with unshared categories, the function does what follows:
#' it divides the BR coefficient(s) by the number of unshared categories plus 0.5.
#' The latter addition is simply a means to be still able to penalize coefficient(s) arising from assemblages having just one unshared category.
#' Also note that joint absences will have no weight on the penalization of the coefficient(s).
#' In case of assemblages sharing all their categories, the corrected coefficient(s) turns out to be equal to the uncorrected one.\cr
#'
#' By setting the parameter 'clust' to TRUE, the units for which the BR coefficients have been calculated will be clustered.
#' Notice that the clustering is based on a dissimilarity matrix which is internally calculated as the maximum values of the BR coefficient (i.e., 200 for the normal values, 1 for the rescales values) minus the BR coefficient.
#' This allows a simpler reading of the dendrogram which is produced by the function, where the less dissimilar (i.e., more similar) units will be placed at lower
#' levels, while more dissimilar (i.e., less similar) units will be placed at higher levels within the dendrogram.\cr
#'
#' The latter depicts the hierarchical clustering based (by default) on the Ward's agglomeration method; rectangles identify the selected cluster partition.
#' Besides the dendrogram, a silhouette plot is produced, which allows to measure how 'good' is the selected cluster solution.\cr
#'
#' As for the latter, if the parameter 'part' is left empty (default), an optimal cluster solution is obtained.
#' The optimal partition is selected via an iterative procedure which locates at which cluster solution the highest average silhouette width is achieved.
#' If a user-defined partition is needed, the user can input the desired number of clusters using the parameter 'part'.
#' In either case, an additional plot is returned besides the cluster dendrogram and the silhouette plot; it displays a scatterplot in which the cluster solution (x-axis)
#' is plotted against the average silhouette width (y-axis). A black dot represent the partition selected either by the iterative procedure or by the user.\cr
#'
#' Notice that in the silhouette plot, the labels on the left-hand side of the chart show the units' names and the cluster number to which each unit is closer.\cr
#'
#' The silhouette plot is obtained from the 'silhouette()' function out from the 'cluster' package (https://cran.r-project.org/web/packages/cluster/index.html).\cr
#'
#' For a detailed description of the silhouette plot, its rationale, and its interpretation, see:\cr
#' Rousseeuw P J. 1987. "Silhouettes: A graphical aid to the interpretation and validation of cluster analysis", Journal of Computational and Applied Mathematics 20, 53-65 (http://www.sciencedirect.com/science/article/pii/0377042787901257).\cr
#'
#' The function also provides a Cleveland's dot plots that represent by-cluster proportions. The clustered units are grouped according to their cluster membership,
#' the frequencies are summed, and then expressed as percentages. The latter are represented by the dot plots, along with the average percentage. The latter
#' provides a frame of reference to understand which percentage is below, above, or close to the average. The raw data on which the plots are based are
#' stored within the list returned by the function (see below).\cr
#'
#' The function also returns a list storing the following:\cr
#'
#' -$BR_similarity_matrix: similarity matrix showing the BR coefficients;\cr
#' -$BR_distance_matrix: dissimilarity matrix on which the hierarchical clustering is performed (if selected);\cr
#' -$avr.silh.width.by.n.of.clusters: average silhouette width by number of clusters  (if clustering is selected);\cr
#' -$partition.silh.data: silhouette data for the selected partition (if clustering is selected);\cr
#' -$data.w.cluster.membership: copy of the input data table with an additional column storing the cluster membership for each row (if clustering is selected);\cr
#' -$by.cluster.proportion: data table showing the proportion of column categories across each cluster; rows sum to 100 percent (if clustering is selected).\cr
#'
#' @param data: dataframe containing the dataset (note: assemblages in rows, variables in columns).
#' @param which: takes "rows" (default) if the user wants the coefficients be calculated for the row categories, "cols" if the users wants the coefficients be calculated for the column categories.
#' @param correction: takes FALSE (default) if the user does not want the coefficients to be corrected, while TRUE will provide corrected coefficients.
#' @param rescale: takes FALSE if the user does NOT want the coefficients to be rescaled between 0.0 and 1.0 (i.e., the user will get the original version of the Brainerd-Robinson coefficient (spanning from 0 [maximum dissimilarity] to 200 [maximum similarity]), while TRUE (default) will return rescaled coefficient.
#' @param clust: TRUE (default) or FALSE if the user does or does not want a agglomerative hierarchical clustering to be performed.
#' @param part: desired number of clusters; if NULL (default), an optimal partition is calculated (see Details).
#' @param aggl.meth: agglomeration method ("ward.D2" by default).
#' @param oneplot: TRUE (default) or FALSE if the user wants or does not want the plots to be visualized in a single window.
#' @param cex.dndr.lab: set the size of the labels used in the dendrogram.
#' @param cex.sil.lab: set the size of the labels used in the silhouette plot.
#' @param cex.dot.plt.lab: set the size of the labels used in the Cleveland's dot charts representing the by-cluster proportions.
#' @keywords similarity
#' @export
#' @examples
#' data(assemblage)
#' coeff <- BRsim(data=assemblage, correction=FALSE, rescale=TRUE, clust=TRUE, oneplot=FALSE)
#'
#' library(archdata) #load the 'archdata' package
#' data(Nelson) #load the 'Nelson' dataset out of the 'archdata' package
#' table <- as.data.frame(as.matrix(Nelson[,3:7])) #build a table to examine
#' res <- BRsim(table, which="rows", clust=TRUE, oneplot=FALSE) # perform the analysis and store the results in the 'res' object
#'
BRsim <- function(data, which="rows", correction=FALSE, rescale=TRUE, clust=TRUE, part=NULL, aggl.meth="ward.D2", oneplot=TRUE, cex.dndr.lab = 0.85, cex.sil.lab = 0.75, cex.dot.plt.lab = 0.80) {
  x <- data

  ifelse(which=="rows", x <- x, x <- as.data.frame(t(x)))

  rd <- dim(x)[1]

  results <- matrix(0,rd,rd)

  if (correction == T){
    for (s1 in 1:rd) {
      for (s2 in 1:rd) {
        zero.categ.a <-length(which(x[s1,]==0))
        zero.categ.b <-length(which(x[s2,]==0))
        joint.absence <-sum(colSums(rbind(x[s1,], x[s2,])) == 0)
        if(zero.categ.a==zero.categ.b) {
          divisor.final <- 1
        } else {
          divisor.final <- max(zero.categ.a, zero.categ.b)-joint.absence+0.5
        }
        results[s1,s2] <- round((1 - (sum(abs(x[s1, ] / sum(x[s1,]) - x[s2, ] / sum(x[s2,]))))/2)/divisor.final, digits=3)
      }
    }
  } else {
    for (s1 in 1:rd) {
      for (s2 in 1:rd) {
        results[s1,s2] <- round(1 - (sum(abs(x[s1, ] / sum(x[s1,]) - x[s2, ] / sum(x[s2,]))))/2, digits=3)
      }
    }
  }

  rownames(results) <- rownames(x)

  colnames(results) <- rownames(x)

  col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white", "cyan", "#007FFF", "blue", "#00007F"))

  if (rescale == FALSE) {
    upper <- 200
    results <- results * 200
  } else {
    upper <- 1.0
  }

  if(clust==TRUE & oneplot==TRUE){
    m <- rbind(c(1,2), c(3,4), c(5,6))
    layout(m)
  }

  corrplot(results, method="square", addCoef.col="red", is.corr=FALSE, cl.lim = c(0, upper), col = col1(100), tl.col="black", tl.cex=0.8)

  if(clust==TRUE){
    dist.matrix <- as.dist(max(results) - results, diag=TRUE, upper=TRUE)

    fit <- hclust(dist.matrix, method = aggl.meth)

    #calculate the max number of clusters to be later used in the loop which iterates the calculation of the average silhouette width
    max.ncl <- rd-1

    #create a slot to store the values of the average silhouette value at increasing numbers of cluster solutions
    #this will be used inside the loop
    sil.width.val <- numeric(max.ncl - 1)

    #create a slot to store the increasing number of clusters whose silhouette is iteratively computed
    sil.width.step <- c(2:max.ncl)

    #calculate the average silhouette width at each increasing cluster solution
    for (i in min(sil.width.step):max(sil.width.step)) {
      counter <- i - 1
      clust.sol <- silhouette(cutree(fit, k = i), dist.matrix)
      sil.width.val[counter] <- mean(clust.sol[, 3])
    }

    #create a dataframe that stores the number of cluster solutions and the corresponding average silhouette values
    sil.res <- as.data.frame(cbind(sil.width.step, sil.width.val))

    #if the user does not enter the desired partition, the latter is equal to the optimal one, ie the one with the largest average silhouette;
    #otherwise use the user-defined partition
    ifelse(is.null(part)==TRUE,
           select.clst.num <- sil.res$sil.width.step[sil.res$sil.width.val == max(sil.res$sil.width.val)],
           select.clst.num <- part)

    #create the final silhouette data using the partition established at the preceding step
    final.sil.data <- silhouette(cutree(fit, k = select.clst.num), dist.matrix)

    #add the cluster membership to a copy of the input dataframe
    data.w.cluster.membership <- x
    data.w.cluster.membership$clust<- assignCluster(x, x, cutree(fit, k = select.clst.num))

    #aggregate the rows of the input dataframe by summing by cluster membership
    #x[-ncol(x)] tells R to work on all the dataframe exclusing the last columnn since, as per the preceding step, it contains the
    #cluster membership
    sum.by.clust <- aggregate(data.w.cluster.membership[-ncol(data.w.cluster.membership)], list(data.w.cluster.membership$clust), sum)[,-1]

    #number of rows of the preceding dataframe
    nrows <- nrow(sum.by.clust)

    #add the columns total as a new row of the preceding dataframe
    sum.by.clust[nrows+1,] <- colSums(sum.by.clust, dims=1)

    #turn the counts to percentages
    prop.by.clust <- as.data.frame(round(prop.table(as.matrix(sum.by.clust), 1)*100,3))

    #add row names to the above dataframe:
    #add reference to cluster membership to all the rows but the last one
    rownames(prop.by.clust)[seq(1,nrows,1)] <- paste0("cluster ", seq(1,select.clst.num, 1))

    #add reference to the average percentage as name of the last row of the dataframe
    rownames(prop.by.clust)[nrows+1] <- "average"

    #plot the dendrogram
    plot(fit, main = "Clusters Dendrogram based on dissimilarity (1-BR sim. coeff.)",
         sub = paste0("\nAgglomeration method: ", aggl.meth),
         xlab = "",
         cex = cex.dndr.lab,
         cex.main = 0.9,
         cex.sub = 0.75)

    #add the rectangles representing the selected partitions
    solution <- rect.hclust(fit, k = select.clst.num, border="black")

    #modify the row names of the final silhouette dataframe to also store the cluster to which each feature is closest
    rownames(final.sil.data) <- paste(rownames(x), final.sil.data[, 2], sep = "_")

    #plot the silhouette chart
    plot(final.sil.data,
         cex.names = cex.sil.lab,
         max.strlen = 30,
         nmax.lab = rd + 1,
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

    #plot the Cleveland's dotcharts of the proportions by cluster
    dotchart(as.matrix(prop.by.clust),
             main="By-cluster proportions",
             xlab="%",
             xlim=c(0,100),
             pt.cex=0.9,
             cex=cex.dot.plt.lab,
             pch=20)

    dotchart(as.matrix(t(prop.by.clust)),
             main="By-cluster proportions",
             xlab="%",
             xlim=c(0,100),
             pt.cex=0.9,
             cex=cex.dot.plt.lab,
             pch=20)

  }

  if(clust==FALSE) {
    dist.matrix <- NULL
    sil.res <- NULL
    final.sil.data <- NULL
    prop.by.clust <- NULL
    data.w.cluster.membership <- NULL
  }

  to.be.returned <- list("BR_similarity_matrix"=results,
                         "BR_distance_matrix"=dist.matrix,
                         "avr.silh.width.by.n.of.clusters"=sil.res,
                         "partition.silh.data"=final.sil.data,
                         "data.w.cluster.membership"=data.w.cluster.membership,
                         "by.cluster.proportion"=prop.by.clust)

  # restore the original graphical device's settings
  par(mfrow = c(1,1))

  return(to.be.returned)
}
