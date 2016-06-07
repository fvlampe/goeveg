#' Draw species response curves
#' @description This function draws species response curves for rough interpretation of a species response to environmental variables or ordination axes.
#' It is based on \code{\link[stats]{smooth.spline}} which fits a cubic smoothing spline to the supplied data.
#' For drawing multiple curves into one plot use \code{\link{specresponses}}.
#' @param species Vector containing species abundances (per plot).
#' @param var Vector containing environmental variable (per plot) \strong{OR} \code{vegan} ordination result object if \code{method = "ord"}.
#' @param main Optional: Main title.
#' @param xlab Optional: Label of x-axis.
#' @param method Method defining the type of variable. Default \code{method = "env"} fits a response curve to environmental variables. Alternatively \code{method = "ord"} fits a response along ordination axes.
#' @param axis Ordination axis (only if \code{method = "ord"}).
#' @param df Desired equivalent number of degrees of freedom (trace of the smoother matrix). Default = 5.
#' @param ylog If set on \code{TRUE} the y-axis is displayed on a log-scale.
#' @section Details:
#' For response curves based on environmental variables the argument \code{var} takes a single vector containing the variable corresponding to the species abundances.
#'
#' For a response to ordination axes (\code{method = "ord"}) the argument \code{var} requires a \code{vegan} ordination result object (e.g. from \code{\link[vegan]{decorana}}, \code{\link[vegan]{cca}}, \code{\link[vegan]{rda}} or \code{\link[vegan]{metaMDS}}).
#' First axis is used as default.
#' @seealso \code{\link{specresponses}}
#' @examples
#' ## Draw species response curve on environmental variable
#' specresponse(meadows$CarxDemi, env$pH)
#'
#' ## Draw species response curve on environmental variable with
#' ## custom labels, lower df and log-scaled y-axis
#' specresponse(meadows$CarxEchi, env$Ca, main = "Carex echinata",
#'        xlab = "Ca", df = 3, ylog = T)
#'
#' ## Draw species response curve on ordination axes
#' specresponse(meadows$CarxDemi, meadows.dca, method = "ord")
#' specresponse(meadows$CarxDemi, meadows.dca, method = "ord", axis = 2)
#' @author Friedemann Goral \email{fgoral@gwdg.de} and Jenny Schellenberg
#' @export

specresponse <- function(species, var, main, xlab, method="env", axis=1, df=5, ylog = FALSE) {

  if(missing(main)) {
    main <- deparse(substitute(species))
  }

  if(method == "env") {

    if(missing(xlab)) {
      xlab <- deparse(substitute(var))
    }

    if(ylog == TRUE) {
        plot(var[species>0], species[species>0], pch=19, col="dimgrey", main = main,
             xlab = xlab, ylab="Abundance", log="y")

        lines(smooth.spline(var[species>0], species[species>0], df=df), col=2)


    } else {

        plot(var[species>0], species[species>0], pch=19, col="dimgrey", main = main,
             xlab = xlab, ylab="Abundance")

        lines(smooth.spline(var[species>0], species[species>0], df=df), col=2)
    }

  } else if(method=="ord") {

    frame <- data.frame(cbind(scores = scores(var, display="sites", choices=axis), species=species))
    frame <- frame[frame$species>0, ]
    frame$scores <- frame$scores + abs(min(frame$scores))

    if(missing(xlab)) {
      xlab <- paste("Axis", axis, "sample scores")
    }

    if(ylog == TRUE) {

        plot(frame$scores, frame$species, pch=19, col="dimgrey",
             ylab="Abundance", xlab=xlab, main=main, log="y")

        lines(smooth.spline(frame$scores, frame$species, df=df), col=2)

    } else {

        plot(frame$scores, frame$species, pch=19, col="dimgrey",
             ylab="Abundance", xlab=xlab, main=main)

        lines(smooth.spline(frame$scores, frame$species, df=df), col=2)

    }

  } else {
    print("FATAL: Method unknown.")
  }
}




