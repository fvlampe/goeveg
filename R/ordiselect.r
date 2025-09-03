#' Species selection for ordination plots
#'
#' @description This function simplifies the selection of relevant species in ordination diagrams. It works with result objects from the \code{vegan} package. The selection can be based upon cover abundances, frequency values and/or species fit to multivariate analysis (see Details).
#' The result is a vector of names of the selected species and can be used for the \code{select} argument in ordination plots.
#' @param matrix Community data, a matrix-like object with samples in rows and species in columns.
#' @param ord \code{vegan} ordination result object (e.g. from \code{\link[vegan]{decorana}}, \code{\link[vegan]{cca}} or \code{\link[vegan]{metaMDS}}).
#' @param ablim Proportion of species with highest abundances to be displayed. Value between 0 and 1. Use negative sign for selection of lowest abundances, i.e. rarest species.
#' @param fitlim Proportion of species with best fit to be displayed. Value between 0 and 1.
#' @param choices Axes shown.
#' @param freq Whether to use cover abundances (= default) or frequencies of \code{matrix}. If \code{TRUE}, frequencies of species are used.
#' @param na.rm Set to \code{TRUE} if your ordination object contains NA (e.g. due to selection)
#' @param method The species fit method: \code{"axes"} (= default) or \code{"factors"}. See details for methods.
#' @param env Fitted environmental variables (result object of \code{\link[vegan]{envfit}}) for \code{method = "factors"}. Only factor variables are used.
#' @param p.max Significance limit for variables used in \code{method = "factors"}.
#' @section Details:
#' Two methods for species fit are implemented.
#' \itemize{\item In \code{method = "axes"} (default) species scores are used for selecting best fitting species. The basic assumption is that species that show high correlations to ordination axes provide a good fit to the assumed gradients, Hence high scores along ordination axes mean high correlation. All species with highest axis scores, defined by the threshold given in argument \code{fitlim}, will be filtered from the total ordination result.
#' \item In \code{method = "factors"}, Euclidean distances between species and environmental variable centroids are calculated. Only factor variables are used from \code{\link[vegan]{envfit}} output. The species with smallest distances, defined by \code{fitlim} argument as a threshold, will be filtered from the ordination result.
#' The \code{p.max} argument allows selection of only significant variables, default is \code{p.max = 0.05}.}
#'
#' The species fit methods work well both in eigenvalue-based and in distance-based ordinations and provide good option of objective reduction of visible species in ordination plot for better interpretation issues.
#' If axes fit should be applied on distance-based ordination, species scores need to be calculated during the analysis, e.g. by selecting \code{wascores = TRUE} in \code{\link[vegan]{metaMDS}}. It is mostly recommendable to combine the species fit limit with an abundance limit to avoid overinterpretation of rare species.
#'
#' For the abundance limit, note that the final proportion of the selected species may be higher than the indicated proportion if there are identical values in the abundances.
#' For selection of least abundant (rarest) species you can use a negative sign, e.g. \code{ablim = -0.3} for the 30 percent least abundant species.
#'
#' If both limits are defined only species meeting both conditions are selected.
#' If no limit is defined for one of the arguments \code{ablim, fitlim}, all species are displayed.
#'
#' The default for \code{matrix} input is a cover-abundance-matrix. This matrix should also be used for ordination.
#'
#' @return
#' A vector of variable length containing the names of selected species from matrix.
#' @examples
#' ## Calculate DCA
#' library(vegan)
#' scheden.dca <- decorana(schedenveg)
#'
#' ## Select the 30% most abundant species and call the result
#' limited <- ordiselect(schedenveg, scheden.dca, ablim = 0.3)
#' limited
#'
#' # Use the result in plotting
#' plot(scheden.dca, display="n")
#' points(scheden.dca, display="sites")
#' points(scheden.dca, display="species",
#'    select = limited, pch = 3, col = "red", cex = 0.7)
#' ordipointlabel(scheden.dca, display="species",
#'     select = limited, col="red", cex = 0.7, add = TRUE)
#'
#'
#' ## Select the 70% of the species with the best fit to the axes (highest species scores)
#' ## AND belonging to the 30% most frequent species
#' limited <- ordiselect(schedenveg, scheden.dca, ablim = 0.3,
#'    fitlim = 0.7, freq = TRUE)
#'
#' ## Select the 30% least frequent species and call the result
#' limited <- ordiselect(schedenveg, scheden.dca, ablim = -0.3, freq = TRUE)
#' limited
#'
#' ## Select the 20% of species with the best fit to community assignment
#' ## AND belonging to the 50% most abundant
#' ## in NDMS for axes 1 & 3
#' nmds <- metaMDS(schedenveg, k = 3)   # run NMDS
#' env13 <- envfit(nmds, schedenenv, choices = c(1, 3))
#' limited13 <- ordiselect(schedenveg, nmds, method = "factors",
#'                        fitlim = 0.1, ablim = 1,
#'                        choices = c(1,3), env = env13)
#'
#' # Use the result in plotting
#' plot(nmds, display="sites", choices = c(1, 3))
#' plot(env13, p.max = 0.05)
#' points(nmds, display="species", choices = c(1,3),
#'     select = limited13, pch = 3, col="red", cex=0.7)
#' ordipointlabel(nmds, display="species", choices = c(1,3),
#'     select = limited13, col="red", cex=0.7, add = TRUE)
#'
#' @author Friedemann von Lampe (\email{fvonlampe@uni-goettingen.de}), Jenny Schellenberg
#' @export
#' @import stats
#' @importFrom fields rdist
#' @importFrom vegan scores


ordiselect <-  function(matrix, ord, ablim = 1, fitlim = 1, choices = c(1, 2), freq = FALSE, na.rm = FALSE, method = "axes", env, p.max = 0.05) {
  if(!is.data.frame(matrix)) {
    matrix <- data.frame(matrix)
  }

  if(freq == F) {
    abund <- apply(matrix, 2, sum)
  } else {
    abund <- apply(matrix>0, 2, sum)
  }
  
  abneg = FALSE
  if(ablim < 0) abneg = TRUE
  
  # ablim = 40/ncol(matrix)  # possibility to calculate automatic abundance limit for a maximum of 40 species

  scores <- data.frame(scores(ord, display = "species", choices = choices))

  if(method == "axes") {
    # Calculate absolute values of scores
    scores_abs <- abs(scores)
    # Extract maximum score value per species
    scores_max <- apply(scores_abs, 1, max)

    # Select according to abundance and fit limits
    if(ablim >= 0 & ablim <= 1) {
      selected <- rownames(scores[abund >= quantile(abund, 1 - ablim) &
                      (scores_max >= quantile(scores_max, 1 - fitlim, na.rm = na.rm))
                    , ])
    } else if (ablim < 0 & ablim >= -1) {
      ablim <- abs(ablim)
      selected <- rownames(scores[abund <= quantile(abund, ablim) &
                                    (scores_max >= quantile(scores_max, 1 - fitlim, na.rm = na.rm))
                                  , ])
    } else {
      stop("Abundance limit must be between 0 and 1 (for most abundant species) OR -1 and 0 (for least abundant species)")
    }

    per <- round(length(selected) / length(abund) * 100, digits = 1)
    
    # Messages
    print(paste0(length(selected), " species selected (", per, "% of total number of species)."))
    
    if(fitlim < 1 & ablim < 1) {
      print(paste0("All species selected which belong to the ", ablim*100, "% ", ifelse(abneg == TRUE, "least ", "most "), 
                   ifelse(freq == FALSE, "abundant ", "frequent "),
                   "species and to the ", fitlim*100, "% of species with the highest absolute axis scores."))
    }
    
    if(fitlim == 1 & ablim < 1) {
      print(paste0("All species selected which belong to the ", ablim*100, "% ", ifelse(abneg == TRUE, "least ", "most "), 
                   ifelse(freq == FALSE, "abundant ", "frequent "),
                   "species."))
    }
    
    if(fitlim < 1 & ablim == 1) {
      print(paste0("All species selected which belong to the ", fitlim*100, "% of species with the highest absolute axis scores."))
    }
    
    if(fitlim == 1 & ablim == 1) {
      print(paste0("All species selected."))
    }
    
    selected


  }  else if(method=="factors") {
    if(inherits(env, "envfit")) {
        sig <- env$factors$pvals
        if(length(sig) == 0) {
          stop("No factors detected. Check envfit-generated input table
               or use method = 'axes' instead")
          } else if(any(env$factors$pvals < p.max)) {
            # Create dataframe with all centroids
            vars <- data.frame(factor = env$factors$var.id)
            vars
            # Add significance values to all centroids
            sig <- data.frame(factor = names(env$factors$pvals), pval = env$factors$pvals)
            sigtab <- merge(vars, sig, by = "factor")

            # Create dataframe with all centroid scores
            scores_env <- data.frame(scores(env, display = "factors"))

            # Select only scores for signifikant factors
            scores_env <- scores_env[sigtab$pval < p.max, ]

            # Calculate distance between species scores and factor centroids for each species
            euclid <- data.frame(fields::rdist(scores, scores_env))
            names(euclid) <- rownames(scores_env)
            rownames(euclid) <- names(matrix)

            # Calculate minimum distance to any centroid for each species
            euclid_min <- apply(euclid, 1, min)
            euclid_min


            # Select species with lowest distances based on given quantile
            # For positive and negative abundance limits
            if(ablim >= 0 & ablim <= 1) {

              selected <- rownames(euclid[abund >= quantile(abund, 1 - ablim) &
                                          (euclid_min <= quantile(euclid_min, fitlim, na.rm = na.rm))
                                        , ])
            } else if (ablim < 0 & ablim >= -1) {

              ablim <- abs(ablim)
              selected <- rownames(euclid[abund <= quantile(abund, ablim) &
                                            (euclid_min <= quantile(euclid_min, fitlim, na.rm = na.rm))
                                          , ])
            } else {
              stop("Abundance limit must be between 0 and 1 (for most abundant species) OR -1 and 0 (for least abundant species)")
            }

            per <- round(length(selected) / length(abund) * 100, digits = 1)

            # Messages
            print(paste0(length(selected), " species selected (", per, "% of total number of species)."))
            
            if(fitlim < 1 & ablim < 1) {
              print(paste0("All species selected which belong to the ", ablim*100, "% ", ifelse(abneg == TRUE, "least ", "most "), ifelse(freq == FALSE, "abundant ", "frequent "),
                           "species and to the ", fitlim*100, "% of species with the smallest distance to variable centroids."))
            }
            
            if(fitlim == 1 & ablim < 1) {
              print(paste0("All species selected which belong to the ", ablim*100, "% ", ifelse(abneg == TRUE, "least ", "most "), ifelse(freq == FALSE, "abundant ", "frequent "),
                           "species."))
            }
            
            if(fitlim < 1 & ablim == 1) {
              print(paste0("All species selected which belong to the ", fitlim*100, "% of species with the smallest distance to variable centroids."))
            }
            
            if(fitlim == 1 & ablim == 1) {
              print(paste0("All species selected."))
            }
            
            selected
          }
        else {
          stop("No significant factors on the given significance level. Check envfit-generated input table
               or use method = 'axes' instead")
        }
    } else {
      stop("Fitted environmental variables are no result of envfit()")
    }
  }
  else {
      stop("Selected method unknown. Use 'axes' or 'factors'")
  }
}



