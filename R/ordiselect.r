#' Select species for ordination plots
#'
#' @description This function simplifies the selection of relevant species in ordination diagrams. It works with result objects from \code{vegan} package. The selection can be based upon cover abundances, frequency values and/or species fit to multivariate analysis.
#' The resulting object contains a vector of names of the selected species and can be used for the \code{select} argument in ordination plots.
#' @param matrix Community data, a matrix-like object with samples in rows and species in columns.
#' @param ord \code{vegan} ordination result object (e.g. from \code{\link[vegan]{decorana}}, \code{\link[vegan]{cca}}, \code{\link[vegan]{rda}} or \code{\link[vegan]{metaMDS}}).
#' @param ablim Proportion (\code{c(0,1)} of species with highest abundances to be displayed.
#' @param fitlim Proportion (\code{c(0,1)} of species with best fit to be displayed.
#' @param choices Axes shown.
#' @param method The species fit method: \code{"axes"} or \code{"vars"}. See details for methods.
#' @param env Fitted environmental variabes (result object of \code{\link[vegan]{envfit}}). Only used if \code{method = "vars"}.
#' @param p.max Significance limit for variables used in \code{method = "vars"}.
#' @param freq Whether to use cover abundances (= default) or frequencies of \code{matrix}. If \code{TRUE}, frequencies of species are used.
#' @section Details:
#' Two methods for species fit are implemented.
#' \itemize{\item In \code{method = "axes"} species scores are used for selecting best fitting species. This is the default method. The basing assumption is that species that show high correlations to ordination axes have good fit. High scores along ordination axes mean high correlation. In this method, all species with high correlations to ordination axes will be filtered.
#' \item In \code{method = "vars"} environmental variables are used for selecting best fitting species. This is a distance-based approach for showing the species with best species-environment-correlation in ordination diagram. Therefore Euclidean distances between species and environment variable centroids are calculated. Only high-responding species with very close or very far distances were considered.
#' If \code{method = "vars"} is used, the environmental variables need to be fitted with \code{\link[vegan]{envfit}} and the result of this function must be provided to the \code{env} argument.
#'
#' The two described methods work well both in eigenvalue-based and in distance-based ordinations.
#' But note, that the distance-based approach for species fit is recommended for also distance-based methods (e.g. NMDS), as axes are meaningless.
#' If axes fit should be applied on distance-based ordination, species scores need to be calculated during the analysis, e.g. by selecting \code{wascores = TRUE} in \code{\link[vegan]{metaMDS}}.
#' On the other hand, distance calculation may be meaningless in Eigenvalue-based approaches.
#' However, both methods provide good option of objective reduction of visible species in ordination plot for better interpretation issues.
#'
#' The \code{p.max} argument allows selection of only significant variables, default is \code{p.max = 0.05}.
#'
#' The default for \code{matrix} input is a cover-abundance-matrix. This matrix should also be used for ordination.
#'
#' If no limit is defined for one of the arguments \code{ablim, fitlim}, all species are displayed.

#' @examples
#' ## Calculate DCA
#' library(vegan)
#' scheden.dca <- decorana(schedenveg)
#'
#' ## Select the 30% most abundant species and call the result
#' limited <- ordiselect(schedenveg, scheden.dca, ablim = 0.3)
#' limited
#'
#' ## Use the result in plotting
#' plot(scheden.dca, display="n")
#' points(scheden.dca, display="sites")
#' points(scheden.dca, display="species",
#'    select = limited, pch=3, col="red", cex=0.7)
#' ordipointlabel(scheden.dca, display="species",
#'    select = limited, col="red", cex=0.7, add = TRUE)
#'
#' ## Select the 30% most frequent species with 50% best axis fit
#' limited <- ordiselect(schedenveg, scheden.dca, ablim = 0.3,
#'    fitlim = 0.5, freq = TRUE)
#'
#' ## Select the 30% most abundant species with 60% best environmental fit
#' ## in NDMS for axes 1 & 3
#' nmds <- metaMDS(schedenveg, k = 3)   # run NMDS
#' env13 <- envfit(nmds, schedenenv[,2:10], choices=c(1,3))
#' limited13 <- ordiselect(schedenveg, nmds, ablim = 0.3, fitlim = 0.6,
#'    choices = c(1,3), method = "vars", env = env13)
#' @author Friedemann Goral \email{fgoral@gwdg.de} and Jenny Schellenberg
#' @export


ordiselect <-  function(matrix, ord, ablim = 1, fitlim = 1, choices = c(1,2), method = "axes", env, p.max = 0.05, freq = FALSE) {
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
      scores_env <- scores_env[sig<p.max,]

      if(nrow(scores_env) == 0) {

        print("WARNING: No significant environmental variables. Only abundance limit used for selection.")

        names(abund[abund >= quantile(abund, 1 - ablim)])

      } else if(nrow(scores_env) == 1) {

        euclid<-data.frame(rdist(scores_spec, scores_env))
        names(euclid)<-rownames(scores_env)
        rownames(euclid)<-names(matrix)
        rownames(data.frame(euclid[abund > quantile(abund, 1 - ablim) &
                           (euclid <= quantile(unlist(euclid), 0.5*fitlim) |
                            euclid >= quantile(unlist(euclid), 1-0.5*fitlim)), ,drop=FALSE]))

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









