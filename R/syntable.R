#' Synoptic tables and calculation of cluster-wise frequencies, fidelity and
#' differential species character
#'
#' @description
#' Synoptic tables are a tool for the visualization and interpretation of previously
#' defined plant species groups (clusters), e.g. from cluster analysis, classification methods or
#' pre-defined categories, e.g. spatial distribution units.
#' They help to determine characteristic patterning of species occurrences in plant communities
#' by calculating cluster-wise percentage or absolute frequencies, mean/median cover values, fidelity
#' (phi) or differential species character.
#'
#' \code{syntable} calculates synoptic tables, using vegetation sample data and a vector of cluster identity.
#' The vegetation data can either be provided as species-sample matrix (\emph{default}) or as long-format vegetation 
#' data (one row per species occurrence) (\code{long = TRUE}). 
#' 
#' The unordered output table can be sorted automatically with \code{\link[goeveg]{synsort}} function.
#' 
##' @param vegdata A data-frame–like object. Can be:
#' \itemize{
#'   \item default: a species–sample matrix with species in columns and samples in rows.
#'         Species and sample names must be defined as column and row names.
#'   \item with \code{long = TRUE}: a long-format table with at least three columns.
#'         The first three columns must contain sample names, taxon names, and abundances, respectively.
#' }
#' Missing values (NA) will be transformed to 0.
#' If non-numeric abundance values are present, the matrix will be transformed to presence/absence with all non-zero values defined as 1.
#' @param long Logical. If \code{TRUE}, \code{vegdata} is treated as long-format data otherwise as species-sample matrix (\emph{default})
#' @param cluster Integer or character vector/factor with classification cluster identity. Ensure matching order of
#' cluster identity and samples in vegdata for correct allocation of cluster numbers to samples.
#' @param abund Type of abundances. Define whether vegetation data are percentage cover (\code{abund = "percentage"}, default)
#' or presence/absence data (\code{abund = "pa"}, with values 0/1). You may use function \code{\link{cov2per}} to transform
#' cover-abundance values from different scales into percentage cover.
#' @param type Type of synoptic table output \code{type = c("percfreq", "totalfreq", "mean",
#' "median", "diffspec", "phi")}. See Details.
#' @param digits Integer indicating the number of decimal places to be displayed in result tables (default 0)
#'
#' @section Details:
#' For synoptic table calculation, six types are available.
#'   \itemize{
#'   \item \code{type = "percfreq" } Creates a percentage frequency table \emph{(default)}
#'   \item \code{type = "totalfreq" } Creates an absolute frequency table
#'   \item \code{type = "mean" }  Calculates mean of species values given in \code{matrix} per cluster
#'   \item \code{type = "median" }  Calculates median of species values given in \code{matrix} per
#'    cluster
#'   \item \code{type = "diffspec" } Calculates differential character of species according to
#'    Tsiripidis et al. 2009, with resulting character p = positive, n = negative, pn = positive-
#'    negative or no differential character (-). Consider that differential character is always
#'    restricted to some and not necessarily all of the other units, thus considering percentage
#'    frequency is essential for correct interpretation of the diagnostic species character.
#'    This calculation needs at least three groups. Only available for (\code{long = FALSE}).
#'   \item \code{type = "phi" } Calculates fidelity measure phi (algorithm basing on Sokal & Rohlf
#'    1995, Bruelheide 2000). Values are ranging between -1 and 1 with high values near 1 indicating
#'    high fidelity.
#'    }
#'
#' For sorting the output synoptic table, use \code{\link{synsort}} function, providing several
#' options.
#'
#' @return
#' The function returns an (invisible) list of result components.
#'   \item{\code{$syntable }}{unordered synoptic table for given species and clusters}
#'   \item{\code{$samplesize }}{total number of samples per cluster}
#'
#'Additionally for differential species character calculation:
#'   \item{\code{$onlydiff }}{Synoptic table only with differential species}
#'   \item{\code{$others }}{List of non-differential species}
#'   \item{\code{$differentials }}{Lists differential species for each cluster}
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
#' ## 1) Unordered synoptic percentage frequency table
#' percfreq <- syntable(schedenveg, pam1$clustering, abund = "percentage",
#'                          type = "percfreq")
#'                          percfreq                   # view results
#'
#' ## 2) Differential species analysis
#' differential <- syntable(schedenveg, pam1$clustering, abund = "percentage",
#'                              type = "diffspec")
#' # show complete table with differential character of species
#' differential$syntable
#' # list differential species for second cluster
#' differential$differentials[2]
#'
#' ## 3) Synoptic table with phi fidelity
#' phitable <- syntable(schedenveg, pam1$clustering, abund = "percentage",
#'                          type = "phi")
#' phitable
#'
#' ## 4) Synoptic percentage frequency table based on historical classification from 1997
#' percfreq <- syntable(schedenveg, schedenenv$comm, abund = "percentage",
#'                          type = "percfreq")
#' percfreq
#'
#' @export


syntable <- function(vegdata, cluster, long = FALSE, abund = "percentage",
                           type = "percfreq", digits = 0) {
  if (long) {
    
    if (type == "diffspec") {
      stop("Differential species analysis not available for long-format data.")
    }
    
    syntable_long(vegdata, cluster, abund = abund, type = type, digits = digits)
    
  } else {
    
    matrix <- vegdata
    syntable_wide(matrix, cluster, abund = abund, type = type, digits = digits)
  }
}
