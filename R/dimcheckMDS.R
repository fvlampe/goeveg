#' Stress plot/Scree plot for NMDS
#' @description This function provides a simple plot of stress values for a given number of tested dimensions (default \code{k = 6}) in NMDS.
#' This stress plot (or scree plot) shows the decrease in ordination stress with an increase in the number of ordination dimensions.
#' It is based on function \code{\link[vegan]{metaMDS}} (\code{vegan} package) and uses the \code{monoMDS} engine.
#' @param matrix Community data, a matrix-like object with samples in rows and species in columns.
#' @param distance Dissimilarity index used in vegdist.
#' @param k Number of dimensions.
#' @param trymax Maximum number of random configuration for iterative search search of stable solution.
#' @param autotransform Whether to use transformation (see \code{\link[vegan]{metaMDS}}) or not. Default is \code{autotransform = TRUE}.
#' @section Details:
#' Goodness of Non-metric multidimensional scaling (NMDS) is measured by stress value.
#' The lower the stress value, the better fit of original distances/dissimilarities and projected distances in ordination diagram is reached.
#' Stress value depends on dimensionality; it is decreasing with increasing dimensionality. On the other hand, stress-reduction does not mean to maximise interpretation capability.
#' Low-dimensional projections are often better to interprete. and are so preferable for interpretation issues.
#' The stress plot (or sometimes also called scree plot) is a diagnostic plots to explore both, dimensionality and interpretative value.
#' It provides dimension-dependant stress reduction and curve estimate gives indices for meaningful stress reduction with increasing dimensionality.
#' Furthermore, another diagnostic plot for detecting best dimension for projection of NMDS, the Shepard diagram (\code{\link[vegan:goodness.metaMDS]{stressplot}}) is recommended for detecting best dimensionality in NMDS.
#'
#' \cite{Clarke 1993} suggests the following guidelines for acceptable stress values:
#' <0.05 = excellent, <0.10 = good, <0.20 = usable, >0.20 = not acceptable.
#' The plot shows the border of the 0.20 stress value limit. Solutions with higher stress values should be interpreted with caution and those with stress above 0.30 are highly suspect.
#' @examples
#' ## Use of function with default values
#' dimcheckMDS(schedenveg)
#'
#' ## Use of function for testing 10 dimensions
#' dimcheckMDS(schedenveg, k = 10)
#' @seealso \code{\link[vegan]{metaMDS}} \code{\link[vegan:goodness.metaMDS]{stressplot}}
#' @references Clarke, K. R. (1993). Non-parametric multivariate analysis of changes in community structure. \emph{Austral J Ecol} \strong{18:} 117-143.
#' @author Jenny Schellenberg (\email{jschell@gwdg.de}) and Friedemann Goral (\email{fgoral@gwdg.de})
#' @export
#' @import graphics
#' @importFrom vegan metaMDS

dimcheckMDS <- function(matrix, distance = "bray", k = 6,  trymax = 20, autotransform = TRUE) {
  if(!is.data.frame(matrix)) {
    matrix <- data.frame(matrix)
  }
  stress <- 0

  for (i in 1:k) {
  nmds_i <- metaMDS(matrix, distance = distance, k = i, trymax = trymax, engine = "monoMDS", autotransform = autotransform)
  stress[i] <- nmds_i$stress
  }

  plot(seq(1,k,1), stress, main="Stress value in tested dimensions", xlab="Dimension", ylab="Stress",
       ylim=c(0,0.3), type="n")
  lines(seq(1,k,1), stress)
  points(seq(1,k,1), stress, pch=21, col="black", bg = "red", cex = 1.2)
  abline(0.2, 0, col="red", lty = 2)
  print(stress)
}
