#' Plot the stress values of NMDS for tested dimensions
#' @description This function provides a simple plot of stress values for a given number of tested dimensions (default \code{k = 6}) in NMDS.
#' It is based on function \code{\link[vegan]{metaMDS}} (\code{vegan} package) and uses the \code{monoMDS} engine.
#' @param matrix Community data, a matrix-like object with samples in rows and species in columns.
#' @param distance Dissimilarity index used in vegdist.
#' @param k Number of dimensions.
#' @param trymax Maximum number of random configuration for iterative search search of stable solution.
#' @section Details:
#' The plot shows the border of the 0.2 stress value limit.
#' \cite{Clarke 1993} suggests the following guidelines for acceptable stress values:
#' <0.05 = excellent, <0.10 = good, <0.20 = usable, >0.20 = not acceptable.
#' @examples
#' ## Use of function with default values
#' dimcheckMDS(meadows)
#'
#' ## Use of function for testing 10 dimensions
#' dimcheckMDS(meadows, k = 10)
#' @seealso \code{\link[vegan]{metaMDS}} \code{\link[vegan]{stressplot}}
#' @references Clarke, K. R. (1993). Non-parametric multivariate analysis of changes in community structure. \emph{Austral J Ecol} \strong{18:} 117-143.
#' @author Jenny Schellenberg \email{jschell@gwdg.de} and Friedemann Goral \email{fgoral@gwdg.de}
#' @export

dimcheckMDS <- function(matrix, distance = "bray", k = 6,  trymax = 20, , autotransform = FALSE) {
  if(!is.data.frame(matrix)) {
    matrix <- data.frame(matrix)
  }
  stress <- 0

  for (i in 1:k) {
  nmds_i<-metaMDS(matrix, distance = distance, k = i, trymax = trymax, engine = "monoMDS", autotransform = FALSE)
  stress[i]<-nmds_i$stress
  }
  plot(seq(1,k,1), stress, main="Stress value in tested dimensions", xlab="Dimension", ylab="Stress", ylim=c(0,0.3), pch=19, col="black")
  abline(0.2, 0, col="red", lty = 2)
  print(stress)
}
