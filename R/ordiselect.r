#' Select species for ordination plots
#'
#' @description This function simplifies the selection of relevant species in ordination diagrams. It works with results from the \code{vegan} package. The selection can be based upon cover abundance values and/or species fit.
#' The resulting object contains the names of the selected species and can be used for the \code{select} argument in ordination plots.
#' @param matrix Community data, a matrix-like object with samples in rows and species in columns.
#' @param ord \code{Vegan} ordination result object (e.g. from \code{\link[vegan]{decorana}}, \code{\link[vegan]{cca}}, \code{\link[vegan]{rda}} or \code{\link[vegan]{metaMDS}}).
#' @param ablim Proportion of species with highest abundances displayed.
#' @param fitlim Proportion of species with best fit displayed.
#' @param choices Axes shown.
#' @param method The species fit method: \code{"axes"} or \code{"vars"}.
#' @param env Fitted environmental variabes (must be result object of \code{\link[vegan]{envfit}}). Only used if \code{method = "vars"}.
#' @param freq It this is set on \code{TRUE}, the frequency of species is used instead of cover-abundance values.
#' @section Details:
#' Two methods for species fit are implemented.
#' \itemize{\item In eigenvalue-based methods the species fit is usually based on the fit to the axes, i.e. species scores are used (\code{method = "axes"}). This is the default method.
#' \item In distance-based methods axes have no direct meaning so it seams more reasonable to fit the species to the environmental variables (\code{method = "vars"}).
#' Here the Euclidean distance between species and environment variable centroids is calculated.}
#' If \code{method = "vars"} is used, the environmental variables need to be fitted with \code{\link[vegan]{envfit}} and the result of this function must be provided to the \code{env} argument.
#' The function uses only significant variables (p < 0.05).
#' If there are no significant variables, only cover abundance values are used and a warning message is displayed.
#'
#' Despite these recommendations, the two described methods work well both in eigenvalue-based and in distance-based ordinations.
#' If axes fit should be applied on distance-based ordination, species scores need to be calculated during the analysis, e.g. by selecting \code{wascores = TRUE} in \code{\link[vegan]{metaMDS}}.
#'
#' If no limit is defined for one of the arguments all species are displayed.
#' @examples
#' ## Select the 30% most abundant species and call the result
#' limited <- ordiselect(meadows, meadows.dca, ablim = 0.3)
#' limited
#'
#' ## Use the result in plotting with orditorp
#' orditorp(meadows.dca, display = "species", select = limited,
#'    pch = 3, col = "red", cex = 0.7, air = 0.5)
#'
#' ## Select the 30% most frequent species with 50% best axis fit
#' limited <- ordiselect(meadows, meadows.dca, ablim = 0.3,
#'    fitlim = 0.5, freq = T)
#'
#' ## Select the 30% most abundant species with 60% best environmental fit
#' ## in NDMS for axes 1 & 3
#' limited13 <- ordiselect(meadows, nmds3, ablim = 0.3, fitlim = 0.6,
#'    choices = c(1,3), method = "vars", env = env13)
#' @author Friedemann Goral (\email{fgoral@gwdg.de}) and Jenny Schellenberg
#' @export


ordiselect <-  function(matrix, ord, ablim = 1, fitlim = 1, choices = c(1,2), method = "axes", env, freq = FALSE) {
  if(!is.data.frame(matrix)) {
    matrix <- data.frame(matrix)
  }

  if(freq == F) {
    abund <- apply(matrix, 2, sum)
    } else {
    abund <- apply(matrix>0, 2, sum)
    }

  if(method == "axes") {

    scores <- data.frame(scores(ord, display="species", choices=choices))
    scores1 <- scores[,1]
    scores2 <- scores[,2]

    rownames(scores[abund >= quantile(abund, 1 - ablim) &
                    (scores1 >= quantile(scores1, 1-0.5*fitlim) |
                       scores1 <= quantile(scores1,0.5*fitlim) |
                       scores2 >= quantile(scores2, 1-0.5*fitlim) |
                       scores2 <= quantile(scores2,0.5*fitlim))
                  ,])

    }  else if(method=="vars") {

      if(class(env) != "envfit") {
        print("FATAL: Fitted environmental variables are no result of envfit()")
      } else {

      scores_spec <- data.frame(scores(ord, display="species", choices=choices))

      sig <- env$vectors$pvals
      scores_env <- data.frame(scores(env, display="vectors"))
      scores_env <- scores_env[sig<0.05,]

      if(nrow(scores_env) == 0) {

        print("WARNING: No significant environmental variables. Only abundance limit used for selection.")

        names(abund[abund >= quantile(abund, 1 - ablim)])

      } else {

        euclid<-data.frame(rdist(scores_spec, scores_env))
        names(euclid)<-rownames(scores_env)
        rownames(euclid)<-names(matrix)

        best_fit<-euclid[abund > quantile(abund, 1 - ablim) &
                       (euclid <= quantile(unlist(euclid), 0.5*fitlim) |
                       euclid >= quantile(unlist(euclid), 1-0.5*fitlim)),]
        rownames(best_fit[complete.cases(best_fit),])
      }
      }
    }  else {
      print("FATAL: Selected Method unknown")
  }
}









