#' Draw rank-abundance curves
#'
#' @description This function draws a rank-abundance curve for community data.
#' If you wish to draw multiple rank-abundance curves for selected samples use \link{racurves}.
#'
#' @param matrix Community data, a matrix-like object with samples in rows.
#' @param main The main title.
#' @param ylog If set on \code{TRUE} the y-axis is displayed on a log-scale.
#' @examples
#' ## Draw simple rank-abundance curve
#' racurve(meadows)
#'
#' ## Draw simple rank-abundance curve with log-scaled axis
#' racurve(meadows, ylog=T)
#'
#' ## Draw multiple rank-abundance curves for selected samples
#' racurve(meadows, samples=c(22,31,45))
#' @seealso \code{\link{racurves}}
#' @author Friedemann Goral (\email{fgoral@gwdg.de})
#' @export

racurve <-  function(matrix, main = "Rank-abundance diagram", ylog = FALSE) {
  if(!is.data.frame(matrix)) {
    matrix <- data.frame(matrix)
  }

  abund <- apply(matrix, 2, sum)
  sum(abund)
  rel.abund <- abund / sum(abund)

  if(ylog == TRUE) {
    plot(sort(rel.abund, decreasing = T), xlab="Abundance Rank", ylab="Relative abundance",
         main=main, log="y")

  } else {
    plot(sort(rel.abund, decreasing = T), xlab="Abundance Rank", ylab="Relative abundance",
     main=main)
  }
}









