#' Species response curves
#' @description This function fits species response curves to visalize a species response to environmental gradients or ordination axes.
#' It is based on Logistic Regression using Generalised Linear Models (GLMs) or Generalized Additive Models (GAMs) with integrated smoothness estimation.
#' For drawing multiple curves into one plot use \code{\link{specresponses}}.
#' @param species Vector containing species abundances (per plot).
#' @param var Vector containing environmental variable (per plot) \strong{OR} \code{vegan} ordination result object if \code{method = "ord"}.
#' @param main Optional: Main title.
#' @param xlab Optional: Label of x-axis.
#' @param model Defining the assumed species response: Default \code{method = "unimodal"} fits a unimodal response. Other methods are \code{method = "linear"} (linear response), \code{method = "bimodal"} (bimodal response), \code{method = "auto"} (automatic selection based on AIC) and \code{method = "gam"} (using GAM with smoother)
#' @param method Method defining the type of variable. Default \code{method = "env"} fits a response curve to environmental variables. Alternatively \code{method = "ord"} fits a response along an ordination axis.
#' @param axis Ordination axis (only if \code{method = "ord"}).
#' @section Details:
#' For response curves based on environmental variable gradients the argument \code{var} takes a single vector containing the variable corresponding to the species abundances. The direct response of a species to the environmental variable is shown.
#'
#' For a response to ordination axis (\code{method = "ord"}) the argument \code{var} requires a \code{vegan} ordination result object (e.g. from \code{\link[vegan]{decorana}}, \code{\link[vegan]{cca}}, \code{\link[vegan]{rda}} or \code{\link[vegan]{metaMDS}}).
#' Response of species to axis is shown; note that meaning of axis differ in unconstrained and constrained methods.
#' First axis is used as default.
#' @return
#' Returns an (invisible) object containing the model results
#' @seealso \code{\link{specresponses}}
#' @examples
#' ## Draw species response curve on environmental variable
#' specresponse(schedenveg$ArrElat, schedenenv$soil_depth)
#'
#' ## Draw species response curve on environmental variable with
#' ## custom labels
#' specresponse(schedenveg$ArrElat, schedenenv$soil_depth, main = "Arrhenatherum elatius",
#'        xlab = "Soil depth")
#'
#' ## Draw species response curve on ordination axes
#' ## First calculate DCA
#' library(vegan)
#' scheden.dca <- decorana(schedenveg)
#'
#' # Using a linear model on first axis
#' specresponse(schedenveg$ArrElat, scheden.dca, method = "ord", model = "linear")
#' # Using an unimodal model on second axis
#' specresponse(schedenveg$ArrElat, scheden.dca, method = "ord", axis = 2)
#' @author Friedemann Goral (\email{fgoral@gwdg.de}) and Jenny Schellenberg
#' @export

specresponse <- function(species, var, main, xlab, model="unimodal", method="env", axis=1) {

  # Add species variable name as title
  if(missing(main)) {
    main <- deparse(substitute(species))
  }

  # Convert to presence/absence
  species <- as.numeric(species>0)

  if(method == "env") {

    # X-axis labeling
    if(missing(xlab)) {
      xlab <- deparse(substitute(var))
    }
  } else if(method=="ord") {

    # Extract site scores from ordination
    var <- as.numeric(scores(var, display="sites", choices=axis))

    # X-axis labeling
    if(missing(xlab)) {
      xlab <- paste("Axis", axis, "sample scores")
    }
  } else {
    stop("Method unknown.")
  }

  # Main plot
  col <- col2rgb("black")

  plot(var, species, pch = 19,
       col = rgb(col[1], col[2], col[3], maxColorValue = 255, alpha = 50),
       main = main, xlab = xlab, ylab="Probability of occurence")

  # GLM models
  if(model == "unimodal") {

    specresponse <- suppressWarnings(glm(species ~ var + I(var^2),
                                         family="binomial"))

    xneu <- seq(min(var), max(var), len = 101)
    preds <- predict(specresponse, newdata = data.frame(var = xneu),
                     type="response")

  } else if (model == "linear") {

    specresponse <- suppressWarnings(glm(species ~ var,
                                         family="binomial"))

    xneu <- seq(min(var), max(var), len = 101)
    preds <- predict(specresponse, newdata = data.frame(var = xneu),
                     type="response")

  } else if (model == "bimodal") {

    specresponse <- suppressWarnings(glm(species ~ var + I(var^2) +
                                           I(var^3) + I(var^4), family="binomial"))

    xneu <- seq(min(var), max(var), len = 101)
    preds <- predict(specresponse, newdata = data.frame(var = xneu),
                     type="response")


  } else if (model == "auto") {

    full <- suppressWarnings(glm(species ~ var + I(var^2) +
                                   I(var^3), family="binomial"))
    specresponse <- suppressWarnings(step(full, trace = 0))

    xneu <- seq(min(var), max(var), len = 101)
    preds <- predict(specresponse, newdata = data.frame(var = xneu),
                     type="response")

  } else if (model == "gam") {

    specresponse <- suppressWarnings(mgcv::gam(species ~ s(var, k = 6),
                                         family="binomial"))

    xneu <- seq(min(var), max(var), len = 101)
    preds <- predict(specresponse, newdata = data.frame(var = xneu),
                     type="response")

  } else {
    stop("Model unknown.")
  }

  lines(preds ~ xneu, col = "red")

  invisible(specresponse)
}




