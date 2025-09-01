#' Synoptic tables and calculation of group-wise frequencies, fidelity and
#' differential species character
#'
#' @description
#' Synoptic tables summarize previously defined plant community groups, e.g., from
#  cluster analysis, classification methods or pre-defined strata, such as spatial distribution units.
#' 
#' They help to determine characteristic patterning of species occurrences in plant communities
#' by calculating group-wise percentage or absolute frequencies, mean/median cover values, fidelity
#' or differential species character.
#'
#' \code{syntable} calculates synoptic tables, using vegetation sample data and a vector of groups identity.
#' The vegetation data can either be provided as species-sample matrix (\emph{default}) or as long-format vegetation 
#' data (one row per species occurrence) (\code{long = TRUE}). 
#' 
#' The unordered output table can be sorted automatically with \code{\link[goeveg]{synsort}} function.
#' 
##' @param vegdata A data-frame-like object. Either:
#' \itemize{
#'   \item default: a species-sample matrix with species in columns and samples in rows.
#'         Species names must be column names, sample (row) names are optional.
#'   \item with \code{long = TRUE}: a long-format table with at least three columns: sample ID, taxon name, and abundance/cover
#'         (these must be the first three columns).
#' }
#' Missing values (NA) are converted to 0.
#' If non-numeric abundance values are present, the matrix will be transformed to presence/absence with all non-zero values defined as 1.
#' @param long Logical. If \code{TRUE}, \code{vegdata} is treated as long-format data; otherwise as species-sample matrix (\emph{default})
#' @param groups Group identities for samples. For wide data (\emph{default}): a vector/factor of length
#'   \code{nrow(vegdata)} (one group per sample/row). For long data: either a
#'   per-\emph{row} vector (length \code{nrow(vegdata)}), a per-\emph{sample} vector
#'   (named by sample IDs, or in the order of first appearance of samples), or
#'   \code{NULL} when \code{group_col} is used. Within each sample, exactly one group
#'   must be defined. 
#' @param abund Type of abundances: percentage cover (\code{"percentage"}, default)
#'   or presence/absence (\code{"pa"} with values 0/1). Use \code{\link{cov2per}}
#'   to transform cover-abundance scales to percentage cover if needed.
#' @param type Output type. One of \code{c("percfreq","totalfreq","mean","median","diffspec","phi")}.
#'   See \strong{Details}.
#' @param digits Integer indicating the number of decimal places to be displayed in result tables (default 0; for phi 3)
#' @param phi_method Fidelity measure when \code{type = "phi"}. One of:
#'   \itemize{
#'     \item \code{"default"}: binary phi coefficient (classical)
#'     \item \code{"ochiai"}: Ochiai coefficient
#'     \item \code{"uvalue"}: hypergeometric \eqn{u}-value
#'   }
#' @param group_col (Long data only) Optional name of a column
#'   in \code{vegdata} that contains the group labels. When supplied, \code{groups}
#'   may be \code{NULL}.
#'
#' @section Details:
#' For synoptic table calculation, six types are available.
#'   \itemize{
#'   \item \code{type = "percfreq" }: percentage frequency of occurrence per group \emph{(default)}
#'   \item \code{type = "totalfreq" }: absolute frequency (number of plots with presence) per group
#'   \item \code{type = "mean" }  mean cover per group (\code{abund = "percentage"} only)
#'   \item \code{type = "median" }  median cover per group (\code{abund = "percentage"} only)
#'   \item \code{type = "phi" } Calculates species fidelity. The default corresponds to the
#'    binary phi coefficient (Sokal & Rohlf 1995, Bruelheide 2000) with values between -1 and 1.
#'    Alternatively, the Ochiai coefficient or the hypergeometric \eqn{u}-value can be selected via \code{phi_method} 
#'    (see Chytry et al., 2002). Handling of group size effects are currently not implemented. 
#'   \item \code{type = "diffspec" } Calculates differential character of species according to
#'    Tsiripidis et al. 2009, with resulting character p = positive, n = negative, pn = positive-
#'    negative or no differential character (-). Consider that differential character is always
#'    restricted to some and not necessarily all of the other units, thus considering percentage
#'    frequency is essential for correct interpretation of the diagnostic species character.
#'    This calculation needs at least 3 groups.
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
#' percfreq                   # view results
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
#'
#' ## 3b) Hypergeometric u-value
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


syntable <- function(vegdata, 
                     groups = NULL, 
                     abund = "percentage",
                     type = "percfreq", 
                     digits = NULL, 
                     long = FALSE,
                     group_col     = NULL,
                     phi_method = "default"
                     ) {
  
  # Validate args (clear errors)
  type         <- match.arg(type,         c("percfreq","totalfreq","mean","median","diffspec","phi"))
  abund        <- match.arg(abund,        c("percentage","pa"))
  phi_method   <- match.arg(phi_method,   c("default","ochiai","uvalue"))

  
  if (isTRUE(long)) {
    if (type == "diffspec") {
      stop("Differential species analysis ('diffspec') is not available for long-format data.\n",
           "Pivot your data to a wide species-sample matrix first (species in columns).")
    }
    
    # For long data: require either a groups vector (per-row or per-sample) or a group_col
    if (is.null(groups) && is.null(group_col)) {
      stop("For long-format data, provide either 'groups' (per-row or per-sample) or 'group_col' (the name of a column in 'vegdata').")
    }
    
    res <- syntable_long(
      vegdata   = vegdata,
      groups    = groups,
      abund     = abund,
      type      = type,
      digits    = digits,
      phi_method    = phi_method,
      group_col     = group_col
    )
  } else {
    
    # Wide matrix requires a groups vector (per-sample)
    if (is.null(groups)) {
      stop("For wide (matrix) data, 'groups' must be provided (one group per sample/row).")
    }
    
    
    res <- syntable_wide(
      matrix    = vegdata,
      groups    = groups,
      abund     = abund,
      type      = type,
      digits    = digits,
      phi_method    = phi_method
    )
  }
  
  return(invisible(res))
}