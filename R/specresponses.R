#' Draw multiple species response curves
#' @description This function draws multiple species response curves for rough interpretation of species responses to environmental variables or ordination axes.
#' It is based on \code{\link[stats]{smooth.spline}} which fits a cubic smoothing spline to the supplied data.
#' In contrast to \code{\link{specresponse}} this function can draw multiple curves into one plot, but will not draw points.
#' @param matrix Community data, a matrix-like object with samples in rows and species in columns. Response curves are drawn for all selected columns (species).
#' @param var Vector containing environmental variable (per plot) \strong{OR} \code{vegan} ordination result object if \code{method = "ord"}.
#' @param main Optional: Main title.
#' @param xlab Optional: Label of x-axis.
#' @param method Method defining the type of variable. Default \code{method = "env"} fits a response curve to environmental variables. Alternatively \code{method = "ord"} fits a response along ordination axes.
#' @param axis Ordination axis (only if \code{method = "ord"}).
#' @param df Desired equivalent number of degrees of freedom (trace of the smoother matrix).
#' @param bw If set on \code{TRUE} the lines will be drawn in black/white with different line types instead of colours.
#' @section Details:
#' For response curves based on environmental variables the argument \code{var} takes a single vector containing the variable corresponding to the species abundances.
#'
#' For a response to ordination axes (\code{method = "ord"}) the argument \code{var} requires a \code{vegan} ordination result object (e.g. from \code{\link[vegan]{decorana}}, \code{\link[vegan]{cca}}, \code{\link[vegan]{rda}} or \code{\link[vegan]{metaMDS}}).
#' First axis is used as default.
#'
#' A minimum of 10 occurences is recommenced to use response curves. Curves for species with less than 5 occurences are not drawn.
#' It is recommended to filter the vegetation matrix for species with a minimum frequency of 10 before using this function.
#'
#' If you plot a response curve for only one species the use of \code{\link{specresponse}} is recommended for correct labels and the display of points.
#' @seealso \code{\link{specresponse}}
#' @examples
#' ## Species (columns) need to be selected; call names() to get column numbers
#' names(schedenveg)
#' ## Draw multiple species response curves on environmental variable in black/white
#' specresponses(schedenveg[ ,c(9,18,14,19)], schedenenv$height_herb, bw = TRUE)
#'
#' ## Draw multiple species response curves on environmental variable with
#' ## custom x-axis label and lower df
#' specresponses(schedenveg[ ,c(9,18,14,19)], schedenenv$height_herb, df = 3,
#'     xlab = "Height of herb layer (cm)")
#'
#' ## Draw multiple species response curves on ordination axes
#' ## First calculate DCA
#' library(vegan)
#' scheden.dca <- decorana(schedenveg)
#' specresponses(schedenveg[ ,c(9,18,14,19)], scheden.dca, method = "ord")
#' specresponses(schedenveg[ ,c(9,18,14,19)], scheden.dca, method = "ord", axis = 2)
#'
#' ## Plot with manually log-transformed abundances
#' specresponses(log(schedenveg[ ,c(9,18,14,19)]+1), schedenenv$height_herb)
#'
#' @author Friedemann Goral \email{fgoral@gwdg.de}
#' @export

specresponses <- function(matrix, var, main, xlab, method="env", axis=1, df=5, bw = FALSE) {

  if(!is.data.frame(matrix)) {
    matrix <- data.frame(matrix)
  }

  if(missing(main)) {
    main <- "Species response curves"
  }

  if(length(matrix) >= 1) {

    ls <- length(matrix)

    if(method == "env") {

      if(missing(xlab)) {
        xlab <- deparse(substitute(var))
      }

      varrange <- 0
      varrange <- var[matrix[,1]>0]

      if(ls > 1) {
          for(i in 2:ls) {
          varrange <- c(varrange, var[matrix[,i]>0])
          }
      }

      plot(var[matrix[,1]>0], matrix[matrix[,1]>0,1], main = main, type="n",
               xlab = xlab, ylab="Abundance", ylim = c(min(matrix), max(matrix)),
               xlim = c(min(varrange), max(varrange)))

        for(i in 1:ls) {
          if(length(matrix[matrix[,i]>0,i]) < 5) {
            print(paste("FATAL: Only", length(matrix[matrix[,i]>0,i]), "occurences of", names(matrix)[i], ". Response curve cannot be drawn with less than 5 occurences!"))
          } else {
              if(length(matrix[matrix[,i]>0,i]) < 10) {
              print(paste("WARNING: Only", length(matrix[matrix[,i]>0,i]), "occurences of", names(matrix)[i], ". At least 10 occurences recommended to draw response curve."))
            }

            if(bw == T) {
            lines(smooth.spline(var[matrix[,i]>0], matrix[matrix[,i]>0,i], df=df), lty=i)
            } else {
            lines(smooth.spline(var[matrix[,i]>0], matrix[matrix[,i]>0,i], df=df), col=i+1)
            }
          }

        if(bw == T){
          legend("topright", inset=0.02, legend=names(matrix), lty=1:ls,
                 bty = "n", cex = 0.85)
        } else {
          legend("topright", inset=0.02, legend=names(matrix), col=2:(ls+1), lty=1,
                 bty = "n", cex = 0.85)
        }
      }

    } else if(method=="ord") {

      var <- scores(var, display="sites", choices=axis)

      varrange <- var[matrix[,1]>0]

      if(ls > 1) {
        for(i in 2:ls) {
          varrange <- c(varrange, var[matrix[,i]>0])
        }
      }

      if(missing(xlab)) {
        xlab <- paste("Axis", axis, "sample scores")
      }

      plot(var[matrix[,1]>0], matrix[matrix[,1]>0,1], main = main, type="n",
           xlab = xlab, ylab="Abundance", ylim = c(min(matrix), max(matrix)),
           xlim = c(min(varrange), max(varrange)))

      for(i in 1:ls) {
        if(length(matrix[matrix[,i]>0,i]) < 5) {
          print(paste("FATAL: Only", length(matrix[matrix[,i]>0,i]), "occurences of", names(matrix)[i], ". Response curve cannot be drawn with less than 5 occurences!"))
        } else {
          if(length(matrix[matrix[,i]>0,i]) < 10) {
            print(paste("WARNING: Only", length(matrix[matrix[,i]>0,i]), "occurences of", names(matrix)[i], ". At least 10 occurences recommended to draw response curve."))
          }

          if(bw == T) {
            lines(smooth.spline(var[matrix[,i]>0], matrix[matrix[,i]>0,i], df=df), lty=i)
          } else {
            lines(smooth.spline(var[matrix[,i]>0], matrix[matrix[,i]>0,i], df=df), col=i+1)
          }
        }
      }

      if(bw == T){
          legend("topright", inset=0.02, legend=names(matrix), lty=1:ls,
                 bty = "n", cex = 0.85)
        } else {
          legend("topright", inset=0.02, legend=names(matrix), col=2:(ls+1), lty=1,
                 bty = "n", cex = 0.85)
        }

    } else {
      print("FATAL: Method unknown.")
    }
  }  else {
    print("FATAL: No species in matrix.")
  }
}






