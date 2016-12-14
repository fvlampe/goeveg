#' Multiple species response curves
#' @description This function fits  multiple species response curves to visualize species responses to environmental gradients or ordination axes.
#' It is based on Logistic Regression using Generalised Linear Models (GLMs) or Generalized Additive Models (GAMs) with integrated smoothness estimation.
#' In contrast to \code{\link{specresponse}} this function can draw multiple curves into one plot.
#' @param matrix Community data, a matrix-like object with samples in rows and species in columns. Response curves are drawn for all (selected) columns (species).
#' @param var Vector containing environmental variable (per plot) \strong{OR} \code{vegan} ordination result object if \code{method = "ord"}.
#' @param main Optional: Main title.
#' @param xlab Optional: Label of x-axis.
#' @param model Defining the assumed species responses: Default \code{method = "unimodal"} fits unimodal responses. Other methods are \code{method = "linear"} (linear responses), \code{method = "bimodal"} (bimodal responses), \code{method = "auto"} (automatic selection based on AIC) and \code{method = "gam"} (using GAM with smoother).
#' @param method Method defining the type of variable. Default \code{method = "env"} fits a response curve to environmental variables. Alternatively \code{method = "ord"} fits a response along ordination axes.
#' @param axis Ordination axis (only if \code{method = "ord"}).
#' @param points If set on \code{TRUE} the species occurences are shown as points. To avoid overlapping they are shown with vertical offset.
#' @param bw If set on \code{TRUE} the lines will be drawn in black/white with different line types instead of colours.
#' @section Details:
#' For response curves based on environmental gradients the argument \code{var} takes a single vector containing the variable corresponding to the species abundances.
#'
#' For a response to ordination axis (\code{method = "ord"}) the argument \code{var} requires a \code{vegan} ordination result object (e.g. from \code{\link[vegan]{decorana}}, \code{\link[vegan]{cca}}, \code{\link[vegan]{rda}} or \code{\link[vegan]{metaMDS}}).
#' First axis is used as default.
#'
#' If you plot a response curve for only one species the use of \code{\link{specresponse}} is recommended.
#' @seealso \code{\link{specresponse}}
#' @examples
#' ## Species (columns) need to be selected; call names() to get column numbers
#' names(schedenveg)
#' ## Draw multiple species response curves on variable in black/white
#' specresponses(schedenveg[ ,c(9,18,14,19)], schedenenv$height_herb, bw = TRUE)
#'
#'## Draw the same curves based on GAM
#' specresponses(schedenveg[ ,c(9,18,14,19)], schedenenv$height_herb, bw = TRUE, model = "gam")
#'
#' ## Draw multiple species response curves on variable with
#' ## custom x-axis label and points of occurences
#' specresponses(schedenveg[ ,c(9,18,14,19)], schedenenv$height_herb,
#'     xlab = "Height of herb layer (cm)", points = TRUE)
#'
#' ## Draw multiple species response curves on ordination axes
#' ## First calculate DCA
#' library(vegan)
#' scheden.dca <- decorana(schedenveg)
#' specresponses(schedenveg[ ,c(9,18,14,19)], scheden.dca, method = "ord")
#' specresponses(schedenveg[ ,c(9,18,14,19)], scheden.dca, method = "ord", axis = 2)
#'
#'
#' @author Friedemann Goral (\email{fgoral@gwdg.de})
#' @export

specresponses <- function(matrix, var, main, xlab, model = "unimodal", method="env", axis=1, points = FALSE, bw = FALSE) {

  if(!is.data.frame(matrix)) {

    if(missing(main)) {
      main <- deparse(substitute(matrix))
    }

    matrix <- data.frame(matrix)
    ismat <- F

  } else {
    if(missing(main)) {
      main <- "Species response curves"
    }

    ismat <- T
  }

  matrix <- decostand(matrix, method = "pa")

  if(length(matrix) >= 1) {

    ls <- length(matrix)

    if(method == "env") {

      if(missing(xlab)) {
        xlab <- deparse(substitute(var))
      }
    } else if(method == "ord") {

      # Extract site scores from ordination
      var <- as.numeric(scores(var, display="sites", choices=axis))

      # X-axis labeling
      if(missing(xlab)) {
        xlab <- paste("Axis", axis, "sample scores")
      }
    } else {
      stop("Method unknown.")
    }

    plot(var, matrix[,1], main = main, type="n",
         xlab = xlab, ylab="Probability of occurence", ylim = c(0,1))

    for(i in 1:ls) {

      if(length(matrix[matrix[,i]>0,i]) <= 3) {
        # tried warning instead of print
        warning(paste("Only", length(matrix[matrix[,i]>0,i]), "occurences of", names(matrix)[i], "."))
      }

      if(model == "unimodal") {

        specresponse <- suppressWarnings(glm(matrix[,i] ~ var
                                             + I(var^2), family="binomial"))

        xneu <- seq(min(var), max(var), len = 101)
        preds <- predict(specresponse, newdata = data.frame(var = xneu),
                         type="response")

      } else if (model == "linear") {

        specresponse <- suppressWarnings(glm(matrix[,i] ~ var,
                                             family="binomial"))

        xneu <- seq(min(var), max(var), len = 101)
        preds <- predict(specresponse, newdata = data.frame(var = xneu),
                         type="response")

      } else if (model == "bimodal") {

        specresponse <- suppressWarnings(glm(matrix[,i] ~ var + I(var^2) +
                                               I(var^3) + I(var^4), family="binomial"))

        xneu <- seq(min(var), max(var), len = 101)
        preds <- predict(specresponse, newdata = data.frame(var = xneu),
                         type="response")

      }
      else if (model == "auto") {

        full <- suppressWarnings(glm(matrix[,i] ~ var + I(var^2) +
                                       I(var^3), family="binomial"))
        specresponse <- suppressWarnings(step(full, trace = 0))

        xneu <- seq(min(var), max(var), len = 101)
        preds <- predict(specresponse, newdata = data.frame(var = xneu),
                         type="response")

      } else if (model == "gam") {

        specresponse <- suppressWarnings(gam(matrix[,i] ~ s(var, k = 6),
                                             family="binomial"))

        xneu <- seq(min(var), max(var), len = 101)
        preds <- predict(specresponse, newdata = data.frame(var = xneu),
                         type="response")

      } else {
        stop("Model unknown.")
      }

      if(bw == T) {
        if(points == TRUE) {
          if(i > 1) {
            matrix[,i][matrix[,i]==1] <- matrix[,i][matrix[,i]==1] - 0.02*(i-1)
            matrix[,i][matrix[,i]==0] <- matrix[,i][matrix[,i]==0] + 0.02*(i-1)
          }

          col <- col2rgb("black")
          points(var, matrix[,i], pch = i,
                 col = rgb(col[1], col[2], col[3], maxColorValue = 255, alpha = 50))
        }

        lines(preds ~ xneu, lty=i)
      } else {

        if(points == TRUE) {
          if(i > 1) {
            matrix[,i][matrix[,i]==1] <- matrix[,i][matrix[,i]==1] - 0.02*(i-1)
            matrix[,i][matrix[,i]==0] <- matrix[,i][matrix[,i]==0] + 0.02*(i-1)
          }

          col <- col2rgb(i+1)
          points(var, matrix[,i],
                 col = rgb(col[1], col[2], col[3], maxColorValue = 255, alpha = 50),
                 pch = 16)
        }

        lines(preds ~ xneu, col = i+1)
      }
    }

    if(ismat == T) {

      if(bw == T){
        if(points == TRUE) {
          legend("topright", inset=0.02, legend=names(matrix), lty=1:ls, pch=1:ls,
                 bty = "n", cex = 0.85)
        } else {
          legend("topright", inset=0.02, legend=names(matrix), lty=1:ls,
                 bty = "n", cex = 0.85)
        }
      } else {
        legend("topright", inset=0.02, legend=names(matrix), col=2:(ls+1), lty=1,
               bty = "n", cex = 0.85)
      }
    }

  }  else {
    stop("No species in matrix.")
  }
}






