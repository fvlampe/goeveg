#' Species response curves
#' @description This function fits species response curves to visualize species responses to environmental gradients or ordination axes.
#' It is based on Logistic Regression using Generalized Linear Models (GLMs) or Generalized Additive Models (GAMs) with integrated smoothness estimation.
#' The function can draw response curves for single or multiple species.
#' @param species Species data (either a community matrix object with samples in rows and species in columns - response curves are drawn for all (selected) columns; or a single vector containing species abundances per plot).
#' @param var Vector containing environmental variable (per plot) \strong{OR} \code{vegan} ordination result object if \code{method = "ord"}.
#' @param main Optional: Main title.
#' @param xlab Optional: Label of x-axis.
#' @param model Defining the assumed species response: Default \code{model = "auto"} selects the model automatically based on AIC. Other methods are \code{model = "linear"} (linear response), \code{model = "unimodal"} (unimodal response), \code{model = "bimodal"} (bimodal response) and \code{model = "gam"} (using GAM with regression smoother).
#' @param method Method defining the type of variable. Default \code{method = "env"} fits a response curve to environmental variables. Alternatively \code{method = "ord"} fits a response along ordination axes.
#' @param axis Ordination axis (only if \code{method = "ord"}).
#' @param points If set on \code{TRUE} the species occurrences are shown as transparent points (the darker the point the more samples at this x-value). To avoid overlapping they are shown with vertical offset.
#' @param bw If set on \code{TRUE} the lines will be drawn in black/white with different line types instead of colors.
#' @param lwd Optional: Graphical parameter defining the line width.
#' @section Details:
#' For response curves based on environmental gradients the argument \code{var} takes a single vector containing the variable corresponding to the species abundances.
#'
#' For a response to ordination axis (\code{method = "ord"}) the argument \code{var} requires a \code{vegan} ordination result object (e.g. from \code{\link[vegan]{decorana}}, \code{\link[vegan]{cca}} or \code{\link[vegan]{metaMDS}}).
#' First axis is used as default.
#'
#' By default the response curves are drawn with automatic GLM model selection based on AIC out of GLMs with 1 - 3 polynomial degrees (thus excluding bimodal responses which must be manually defined). The GAM model is more flexible and choses automatically between an upper limit of 3 - 6 degrees of freedom for the regression smoother.
#'
#' Available information about species is reduced to presence-absence as species abundances can contain much noise (being affected by complex factors) and the results of Logistic Regression are easier to interpret showing the "probabilities of occurrence".
#' Be aware that response curves are only a simplification of reality (model) and their shape is strongly dependent on the available dataset.
#' @return
#' Returns character string with information on model type and parameters per species.
#'
#' No return if model type predefined \emph{linear} or \emph{unimodal}.
#' @examples
#' ## Draw species response curve for one species on environmental variable
#' ## with points of occurrences
#' specresponse(schedenveg$ArrElat, schedenenv$soil_depth, points = TRUE)
#'
#' ## Draw species response curve on environmental variable with custom labels
#' specresponse(schedenveg$ArrElat, schedenenv$soil_depth, points = TRUE,
#'        main = "Arrhenatherum elatius", xlab = "Soil depth")
#'
#' ## Draw species response curve on ordination axes
#' ## First calculate DCA
#' library(vegan)
#' scheden.dca <- decorana(schedenveg)
#'
#' # Using a linear model on first axis
#' specresponse(schedenveg$ArrElat, scheden.dca, method = "ord", model = "linear")
#' # Using an unimodal model on second axis
#' specresponse(schedenveg$ArrElat, scheden.dca, method = "ord", axis = 2, model = "unimodal")
#'
#' ## Community data: species (columns) need to be selected; call names() to get column numbers
#' names(schedenveg)
#' ## Draw multiple species response curves on variable in black/white
#' specresponse(schedenveg[ ,c(9,18,14,19)], schedenenv$height_herb, bw = TRUE)
#'
#' ## Draw the same curves based on GAM
#' specresponse(schedenveg[ ,c(9,18,14,19)], schedenenv$height_herb, bw = TRUE, model = "gam")
#'
#' ## Draw multiple species response curves on variable with
#' ## custom x-axis label and points of occurrences
#' specresponse(schedenveg[ ,c(9,18,14,19)], schedenenv$height_herb,
#'     xlab = "Height of herb layer (cm)", points = TRUE)
#'
#' ## Draw multiple species response curves on ordination axes
#' specresponse(schedenveg[ ,c(9,18,14,19)], scheden.dca, method = "ord")
#' specresponse(schedenveg[ ,c(9,18,14,19)], scheden.dca, method = "ord", axis = 2)
#'
#' @author Friedemann Goral (\email{fgoral@gwdg.de})
#' @export
#' @import graphics
#' @import stats
#' @import grDevices
#' @importFrom mgcv gam
#' @importFrom vegan scores decostand

specresponse <- function(species, var, main, xlab, model = "auto", method = "env", axis = 1, points = FALSE, bw = FALSE, lwd = NULL) {

  if(!is.data.frame(species)) {

    if(missing(main)) {
      main <- deparse(substitute(species))
    }

    specnames <- deparse(substitute(species))
    species <- data.frame(species)
    ismat <- F

  } else {
    if(missing(main)) {
      main <- "Species response curves"
    }
    specnames <- names(species)
    ismat <- T
  }

  species <- decostand(species, method = "pa")

  if(length(species) >= 1) {

    ls <- length(species)

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

    plot(var, species[,1], main = main, type="n",
         xlab = xlab, ylab="Probability of occurrence", ylim = c(0,1))

    for(i in 1:ls) {

      if(length(species[species[,i]>0,i]) <= 5) {
        # tried warning instead of print
        warning(paste("Only", length(species[species[,i]>0,i]), "occurrences of", names(species)[i], "."))
      }

      if(model == "unimodal") {

        specresponse <- suppressWarnings(glm(species[,i] ~ poly(var, 2),
                                             family="binomial"))

      } else if (model == "linear") {

        specresponse <- suppressWarnings(glm(species[,i] ~ var,
                                             family="binomial"))

      } else if (model == "bimodal") {

        specresponse <- suppressWarnings(glm(species[,i] ~ poly(var, 4),
                                             family="binomial"))

      }
      else if (model == "auto") {

        glm.1 <- suppressWarnings(glm(species[,i] ~ poly(var, 1), family="binomial"))
        glm.2 <- suppressWarnings(glm(species[,i] ~ poly(var, 2), family="binomial"))
        glm.3 <- suppressWarnings(glm(species[,i] ~ poly(var, 3), family="binomial"))
        glm.AIC <- c(extractAIC (glm.1)[2], extractAIC (glm.2)[2],
                     extractAIC (glm.3)[2])

        switch(c(1,2,3)[glm.AIC==min(glm.AIC)],
               {specresponse <- glm.1; deg<-1},
               {specresponse <- glm.2; deg<-2},
               {specresponse <- glm.3; deg<-3})

        print(paste("GLM with", deg, "degrees fitted for", specnames[i]))

      } else if (model == "gam") {

        gam.list <- list()
        gam.AIC <- 0

        for(n in 1:4) {
          gam.list[[paste("gam", n, sep=".")]] <- mgcv::gam(species[,i] ~ s(var, k = n+2), family='binomial')
          gam.AIC[n] <- extractAIC(gam.list[[n]])[2]
        }

        switch(c(1,2,3,4)[gam.AIC==min(gam.AIC)],
               {specresponse <- gam.list[[1]]; deg<-3},
               {specresponse <- gam.list[[2]]; deg<-4},
               {specresponse <- gam.list[[3]]; deg<-5},
               {specresponse <- gam.list[[4]]; deg<-6})

        print(paste("GAM with", deg, "knots fitted for", specnames[i]))

      } else {
        stop("Model unknown.")
      }

      xneu <- seq(min(var), max(var), len = 101)
      preds <- predict(specresponse, newdata = data.frame(var = xneu),
                       type="response")

      if(bw == T) {
        if(points == TRUE) {
          if(i > 1) {
            species[,i][species[,i]==1] <- species[,i][species[,i]==1] - 0.02*(i-1)
            species[,i][species[,i]==0] <- species[,i][species[,i]==0] + 0.02*(i-1)
          }

          col <- col2rgb("black")
          points(var, species[,i], pch = i,
                 col = rgb(col[1], col[2], col[3], maxColorValue = 255, alpha = 50))
        }

        lines(preds ~ xneu, lty=i, lwd = lwd)
      } else {

        if(points == TRUE) {
          if(i > 1) {
            species[,i][species[,i]==1] <- species[,i][species[,i]==1] - 0.02*(i-1)
            species[,i][species[,i]==0] <- species[,i][species[,i]==0] + 0.02*(i-1)
          }

          col <- col2rgb(i)
          points(var, species[,i],
                 col = rgb(col[1], col[2], col[3], maxColorValue = 255, alpha = 50),
                 pch = 16)
        }

        lines(preds ~ xneu, col = i, lwd = lwd)
      }
    }

    if(ismat == T) {

      if(bw == T){
        if(points == TRUE) {
          legend("topright", inset=0.02, legend=names(species), lty=1:ls, pch=1:ls,
                 bty = "n", cex = 0.85)
        } else {
          legend("topright", inset=0.02, legend=names(species), lty=1:ls,
                 bty = "n", cex = 0.85)
        }
      } else {
        legend("topright", inset=0.02, legend=names(species), col=1:ls, lty=1,
               bty = "n", cex = 0.85)
      }
    }

  }  else {
    stop("No species in matrix.")
  }
}






