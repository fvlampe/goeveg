#' Synoptic tables and calculation of group-wise frequencies, fidelity and
#' differential species character
#'
#' @description
#' Synoptic tables are a tool for the visualization and interpretation of previously
#' defined plant species groups, e.g. from cluster analysis, classification methods or
#' pre-defined categories, e.g. spatial distribution units.
#' They help to determine characteristic patterning of species occurrences in plant communities
#' by calculating group-wise percentage or absolute frequencies, mean/median cover values, fidelity
#' (phi) or differential species character.
#'
#' \code{syntable} calculates synoptic tables, using vegetation sample data and a vector of groups identity.
#' The vegetation data can either be provided as species-sample matrix (\emph{default}) or as long-format vegetation 
#' data (one row per species occurrence) (\code{long = TRUE}). 
#' 
#' The unordered output table can be sorted automatically with \code{\link[goeveg]{synsort}} function.
#' 
##' @param vegdata A data-frame-like object. Can be:
#' \itemize{
#'   \item default: a species-sample matrix with species in columns and samples in rows.
#'         Species names must be column names, sample (row) names are optional.
#'   \item with \code{long = TRUE}: a long-format table with at least three columns.
#'         The first three columns must contain sample names, taxon names, and abundances, respectively.
#' }
#' Missing values (NA) will be transformed to 0.
#' If non-numeric abundance values are present, the matrix will be transformed to presence/absence with all non-zero values defined as 1.
#' @param long Logical. If \code{TRUE}, \code{vegdata} is treated as long-format data otherwise as species-sample matrix (\emph{default})
#' @param groups Integer or character vector/factor with groups identity. Ensure matching order of
#' groups identity and samples in vegdata for correct allocation of groups numbers to samples.
#' @param abund Type of abundances. Define whether vegetation data are percentage cover (\code{abund = "percentage"}, default)
#' or presence/absence data (\code{abund = "pa"}, with values 0/1). You may use function \code{\link{cov2per}} to transform
#' cover-abundance values from different scales into percentage cover.
#' @param type Type of synoptic table output \code{type = c("percfreq", "totalfreq", "mean",
#' "median", "diffspec", "phi")}. See Details.
#' @param digits Integer indicating the number of decimal places to be displayed in result tables (default 0)
#' @param phi_method Fidelity measure when \code{type = "phi"}. Choose \code{"default"} for the
#'  classical binary phi coefficient, \code{"cover"} for a cover-weighted phi based on the
#'  correlation between species cover and group membership, \code{"ochiai"} for the Ochiai
#'  coefficient or \code{"uvalue"} for the hypergeometric \eqn{u}-value.
#' @param phi_transform Transformation of cover values when \code{phi_method = "cover"}.
#'  Options are \code{"none"}, \code{"sqrt"} or \code{"log"} (using \code{log(x + 1)}).
#' @param phi_standard Standardisation of group sizes for fidelity calculation. Use \code{"none"}
#'  for the original group sizes, \code{"rarefy"} for random rarefaction to the smallest group
#'  (repeated \code{phi_reps} times) or \code{"adjust"} for analytical adjustment of phi.
#' @param phi_reps Number of repetitions for random rarefaction (default 100).
#'
#' @section Details:
#' For synoptic table calculation, six types are available.
#'   \itemize{
#'   \item \code{type = "percfreq" } Creates a percentage frequency table \emph{(default)}
#'   \item \code{type = "totalfreq" } Creates an absolute frequency table
#'   \item \code{type = "mean" }  Calculates mean of species values given in \code{matrix} per group
#'   \item \code{type = "median" }  Calculates median of species values given in \code{matrix} per
#'    group
#'   \item \code{type = "phi" } Calculates species fidelity. The default corresponds to the
#'    classical phi coefficient (Sokal & Rohlf 1995, Bruelheide 2000) with values between -1 and 1.
#'    Alternatively, a cover-weighted phi based on correlation, the Ochiai coefficient or the
#'    hypergeometric \eqn{u}-value can be selected via \code{phi_method} (see Chytry et al., 2002). Group size effects can
#'    be handled by \code{phi_standard}.
#'    }
#'
#' For sorting the output synoptic table, use \code{\link{synsort}} function, providing several
#' options.
#'
#' @return
#' The function returns an (invisible) list of result components.
#'   \item{\code{$syntable }}{unordered synoptic table for given species and groups}
#'   \item{\code{$samplesize }}{total number of samples per group}
#'
#'Additionally for differential species character calculation:
#'   \item{\code{$onlydiff }}{Synoptic table only with differential species}
#'   \item{\code{$others }}{List of non-differential species}
#'   \item{\code{$differentials }}{Lists differential species for each group}
#'
#'
#' @references
#' Bruelheide, H. (2000): A new measure of fidelity and its application to defining species groups.
#'  \emph{Journal of Vegetation Science} \strong{11}: 167-178. \doi{https://doi.org/10.2307/3236796}
#'
#' Chytry, M., Tichy, L., Holt, J., Botta-Dukat, Z. (2002): Determination of diagnostic species with
#'  statistical fidelity measures. \emph{Journal of Vegetation Science} \strong{13}: 79-90. \doi{https://doi.org/10.1111/j.1654-1103.2002.tb02025.x}
#'
#' Sokal, R.R. & Rohlf, F.J. (1995): Biometry. 3rd edition Freemann, New York.
#'
#' Tsiripidis, I., Bergmeier, E., Fotiadis, G. & Dimopoulos, P. (2009): A new algorithm for the
#' determination of differential taxa. \emph{Journal of Vegetation Science} \strong{20}: 233-240. \doi{https://doi.org/10.1111/j.1654-1103.2009.05273.x}
#'
#' @author Jenny Schellenberg (\email{jschell@gwdg.de}) and Friedemann von Lampe
#' @seealso \code{\link{synsort}}
#' @examples
#' ## Synoptic table of Scheden vegetation data
#' library(cluster)
#' pam1 <- pam(schedenveg, 4)  # PAM clustering with 4 clusters output
#'
#'
#' ## 1) Unordered synoptic percentage frequency table
#' percfreq <- syntable(schedenveg, pam1$clustering, abund = "percentage",
#'                          type = "percfreq")
#'                          percfreq                   # view results
#'
#'
#' ## 2) Differential species analysis
#' differential <- syntable(schedenveg, pam1$clustering, abund = "percentage",
#'                              type = "diffspec")
#' # show complete table with differential character of species
#' differential$syntable
#' # list differential species for second cluster
#' differential$differentials[2]
#'
#'
#' ## 3) Synoptic table with phi fidelity
#' phitable <- syntable(schedenveg, pam1$clustering, abund = "percentage",
#'                          type = "phi")
#' phitable
#' 
#' ## 3b) Cover-weighted phi with square-root transformed cover
#' phicover <- syntable(schedenveg, pam1$clustering, abund = "percentage",
#'                       type = "phi", phi_method = "cover", phi_transform = "sqrt")
#' phicover
#'
#' ## 3c) Ochiai coefficient
#' ochiai <- syntable(schedenveg, pam1$clustering, abund = "percentage",
#'                     type = "phi", phi_method = "ochiai")
#' ochiai
#'
#' ## 3d) Hypergeometric u-value
#' phiu <- syntable(schedenveg, pam1$clustering, abund = "percentage",
#'                  type = "phi", phi_method = "uvalue")
#' phiu
#'
#'
#' ## 4) Synoptic percentage frequency table based on historical classification from 1997
#' percfreq <- syntable(schedenveg, schedenenv$comm, abund = "percentage",
#'                          type = "percfreq")
#' percfreq
#'
#' @export


syntable <- function(vegdata, groups, abund = "percentage",
                           type = "percfreq", digits = 0, long = FALSE,
                           phi_method = "default", phi_transform = "none",
                           phi_standard = "none", phi_reps = 100) {
  
  
  if (isTRUE(long)) {
    if (type == "diffspec") {
      stop("Differential species analysis ('diffspec') is not available for long-format data.\n",
           "Pivot your data to a wide species-sample matrix first (species in columns).")
    }
    
    res <- syntable_long(
      vegdata, groups,
      abund = abund,
      type = type,
      digits = digits,
      phi_method = phi_method,
      phi_transform = phi_transform,
      phi_standard = phi_standard,
      phi_reps = phi_reps
    )
  } else {
    res <- syntable_wide(
      matrix = vegdata, groups,
      abund = abund,
      type = type,
      digits = digits,
      phi_method = phi_method,
      phi_transform = phi_transform,
      phi_standard = phi_standard,
      phi_reps = phi_reps
    )
  }
  
  return(invisible(res))
}