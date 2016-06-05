#' Plot the stress values of NMDS for tested dimensions
#' @description This function provides a simple plot of stress values for a given number of tested dimensions (default \code{k = 6}) in NMDS.
#' It is based on function \code{\link[vegan]{metaMDS}} (\code{vegan} package) and uses the \code{monoMDS} engine.
#' @param matrix Community data, a matrix-like object with samples in rows and species in columns.
#' @param distance Dissimilarity index used in vegdist.
#' @param k Number of dimensions.
#' @param trymax Maximum number of random starts in search of stable solution.
#' @examples
#' ## Use of function with default values
#' dimcheckMDS(meadows)
#'
#' ## Use of function for testing 10 dimensions
#' dimcheckMDS(meadows, k = 10)
#' @seealso \code{\link[vegan]{metaMDS}} \code{\link[vegan]{stressplot}}
#' @author Friedemann Goral (\email{fgoral@gwdg.de}) and Jenny Schellenberg
#' @export

dimcheckMDS <- function(matrix, distance = "bray", k = 6,  trymax = 20) {
  if(!is.data.frame(matrix)) {
    matrix <- data.frame(matrix)
  }
  stress <- 0

  for (i in 1:k) {
  nmds_i<-metaMDS(matrix, distance = distance, k = i, trymax = trymax, engine = "monoMDS")
  stress[i]<-nmds_i$stress
  }
  plot(seq(1,k,1), stress, main="Stress value in tested dimensions", xlab="Dimension", ylab="Stress", ylim=c(0,0.3), pch=19, col="black")
  abline(0.2, 0, col="red")
  print(stress)
}
