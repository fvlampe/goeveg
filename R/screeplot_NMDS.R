#' Scree plot/Stress plot for NMDS
#' @description This function provides a plot of stress values against a given number of tested dimensions (default \code{k = 6}) in NMDS.
#' This scree plot (or stress plot) shows the decrease in ordination stress with an increase in the number of ordination dimensions.
#' It is based on function \code{\link[vegan]{metaMDS}} (\code{vegan} package) and uses the \code{monoMDS} engine.
#' @param matrix Community data, a matrix-like object with samples in rows and species in columns.
#' @param distance Dissimilarity index used in vegdist.
#' @param k Number of dimensions (default \code{k = 6}).
#' @param trymax Maximum number of random configuration for iterative search search of stable solution.
#' @param autotransform Whether to use transformation (see \code{\link[vegan]{metaMDS}}) or not. Default is \code{autotransform = TRUE}.
#' @section Details:
#' The simplest indicator for the goodness of non-metric multidimensional scaling (NMDS) is the stress value. Stress is a value between 0 and 1 and expresses a proportion between the distance in the original dissimilarity matrix and the fitted distance in ordination space.
#' The lower the stress value, the better is the fit. Details and exact formula are found under function \code{\link[vegan]{monoMDS}}.
#' Stress value depends on dimensionality: it is decreasing with increasing dimensionality. On the other hand, stress-reduction does not mean to maximize interpretation capability.
#' Low-dimensional projections are often better to interpret and are so preferable for interpretation issues.
#'
#' A scree plot (or sometimes also called stress plot) is a diagnostic plot to explore both, dimensionality and interpretative value.
#' Often the 'elbow' of the plot is an indicator to determine the optimal number of dimensions to capture most information contained in the data.
#' However, for ecological data this is rarely seen, so that the following rule of thumb under consideration of the interpretability can be used.
#'
#' \cite{Clarke 1993} suggests the following guidelines for stress values:
#' <0.05 = excellent, <0.10 = good, <0.20 = usable without putting too much reliance on details, >0.20 = could be dangerous to interpret, > 0.35 = samples effectively randomly placed.
#' The plot shows the border of the 0.20 stress value limit. Solutions with higher stress values should be interpreted with caution and those with stress above 0.35 are highly suspect.
#'
#' It should be taken into account that the stress value naturally increases with the number of samples and/or variables: “The greater the number of samples, the harder it will be to reflect the complexity of their inter-relationship in a two-dimensional plot” (Clarke 1993, p. 125).
#' Furthermore high stress values can be attributed to all points or only to few or even single points, that represent samples different to the others.
#'
#' The scree plot is not an exclusive method to determine the optimal number of dimensions.
#' The effect of individual points on the stress value can be explored with the \code{\link[vegan:goodness.metaMDS]{goodness}}-function.
#' Large values indicate poor fit and can be easily visualized in an ordination diagram.
#'
#' Another diagnostic plot is the the Shepard diagram (\code{\link[vegan:goodness.metaMDS]{stressplot}}), which is a scatterplot of the (Euclidean) distances in ordination space against the original dissimilarities.
#' However, as the number of pairwise elements increases rapidly with the number of samples, the Shepard diagram will mostly lead to a high ‘correlation-like' fit and lacks therefore a reliable and precise interpretability.
#'
#' @return
#' A numeric vector of length \emph{k} containing stress values for \emph{k} dimensions.
#' @examples
#' ## Use of function with default values
#' screeplot_NMDS(schedenveg)
#'
#' ## Use of function for testing 10 dimensions
#' screeplot_NMDS(schedenveg, k = 10)
#'
#' ## Alternative diagnostic plots
#' library(vegan)
#' nmds <- metaMDS(schedenveg, k = 2)
#'
#' # Draw Shepard plot
#' stressplot(nmds)
#'
#' # Calculate goodness of fit
#' gof <- goodness(object = nmds)
#'
#' # Draw NMDS ordination diagram with sites
#' plot(nmds, display = "sites", type = "n", cex = 0.7)
#' # Add the points with size reflecting goodness of fit (bigger = worse fit)
#' points(nmds, display = "sites", cex = 2*gof/mean(gof))
#'
#' @seealso \code{\link[vegan]{metaMDS}} \code{\link[vegan:goodness.metaMDS]{stressplot}}
#' @references Clarke, K. R. (1993). Non-parametric multivariate analysis of changes in community structure. \emph{Austral J Ecol} \strong{18:} 117-143. \doi{10.1111/j.1442-9993.1993.tb00438.x}
#' @author Friedemann von Lampe, Jenny Schellenberg
#' @export dimcheckMDS screeplot_NMDS
#' @aliases dimcheckMDS
#' @import graphics utils
#' @importFrom vegan metaMDS

screeplot_NMDS <- function(matrix, distance = "bray", k = 6,  trymax = 20, autotransform = TRUE) {
  if(!is.data.frame(matrix)) {
    matrix <- data.frame(matrix)
  }
  stress <- 0

  # create progress bar
  pb <- txtProgressBar(min = 0, max = k, style = 3)

  for (i in 1:k) {
  capture.output(nmds_i <- invisible(metaMDS(matrix, distance = distance, k = i, trymax = trymax, 
                                             engine = "monoMDS", autotransform = autotransform)),
                 file = NULL)

  stress[i] <- nmds_i$stress

  setTxtProgressBar(pb, i)
  }
  
  # Close bar
  close(pb)

  plot(seq(1,k,1), stress, main="Stress value in tested dimensions", xlab="Dimension", ylab="Stress",
       ylim=c(0,0.3), type="b")
  points(seq(1,k,1), stress, pch=21, col="black", bg = "red", cex = 1.2)
  abline(0.2, 0, col="red", lty = 2)

  names(stress) <- 1:k
  stress <- round(stress, 3)
  print(stress)
}

dimcheckMDS = screeplot_NMDS


