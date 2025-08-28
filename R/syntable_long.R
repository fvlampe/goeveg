#' Synoptic tables and calculation of cluster-wise frequencies and fidelity for
#' long-format vegetation databases
#'
#' @description
#' Synoptic tables are a tool for the visualization and interpretation of previously
#' defined plant species groups (clusters), e.g. from cluster analysis, classification methods or
#' pre-defined categories, e.g. spatial distribution units.
#' They help to determine characteristic patterning of species occurrences in plant communities
#' by calculating cluster-wise percentage or absolute frequencies, mean/median cover values or fidelity
#' (phi) values. 
#' 
#' `syntable_long` calculates synoptic tables from long-format vegetation data (one row per species occurrence)
#' where the first column contains the sample identifier. A parallel vector provides the cluster identity for each sample.
#' It is therefore particularly suited for comprehensive vegetation databases, where usual species-sample matrices are very large.
#' The unordered output table can be sorted automatically with \code{\link[goeveg]{synsort}} function.
#'  
#' For data in wide species-by-sample matrix form, use the companion function \code{\link[goeveg]{syntable}}.
#'
#' @param vegdata A data-frame–like object in long format, with at least the columns. The first three columns must contain
#'   sample names, taxon names and abundances, respectively.
#' @param cluster Integer or character vector/factor with classification
#'  cluster identity for each sample. Ensure matching order of cluster identity and samples in vegdata for correct allocation of cluster numbers to samples.
#' @param abund Type of abundances. Define whether the third column contains
#'   percentage cover (`abund = "percentage"`, default) or presence/absence
#'   data (`abund = "pa"`, with values 0/1).
#' @param type Type of synoptic table output
#'   `type = c("percfreq", "totalfreq", "mean", "median", "phi")`. See Details.
#' @param digits Integer indicating the number of decimal places to be
#'   displayed in result tables (default 0).
#'   
#' @section Details:
#' For synoptic table calculation, five types are available.
#'   \itemize{
#'   \item \code{type = "percfreq" } Creates a percentage frequency table \emph{(default)}
#'   \item \code{type = "totalfreq" } Creates an absolute frequency table
#'   \item \code{type = "mean" }  Calculates mean of species values given in \code{matrix} per cluster
#'   \item \code{type = "median" }  Calculates median of species values given in \code{matrix} per
#'    cluster
#'   \item \code{type = "phi" } Calculates fidelity measure phi (algorithm basing on Sokal & Rohlf
#'    1995, Bruelheide 2000). Values are ranging between -1 and 1 with high values near 1 indicating
#'    high fidelity.
#'    }
#'
#' For sorting the output synoptic table, use \code{\link{synsort}} function, providing several
#' options.
#'
#' @return An (invisible) list with components:
#'   \item{\code{$syntable }}{unordered synoptic table for given species and clusters}
#'   \item{\code{$samplesize }}{total number of samples per cluster}
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
#' @author Friedemann von Lampe and Jenny Schellenberg
#' @seealso [syntable()], [synsort()]
#' @import data.table
#' @export
syntable_long <- function(vegdata, cluster, abund = "percentage", type = "percfreq", digits = 0) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required for syntable_long.")
  }
  dt <- data.table::as.data.table(vegdata)
  if (ncol(dt) < 3) {
    stop("Input table must contain at least three columns (sample, taxon, abundance).")
  }
  data.table::setnames(dt, 1:3, c("Sample", "TaxonName", "Abund"))
  
  narep <- FALSE
  if (any(dt$Abund == "", na.rm = TRUE)) {
    dt$Abund[dt$Abund == ""] <- NA
    narep <- TRUE
  }
  if (any(is.na(dt$Abund))) {
    dt$Abund[is.na(dt$Abund)] <- 0
    narep <- TRUE
  }
  if (narep) {
    print("NA and/or empty character values replaced by 0.")
  }
  
  if (any(is.na(cluster))) {
    stop("NA values in cluster not allowed.")
  }
  
  if (any(is.na(suppressWarnings(as.numeric(dt$Abund))))) {
    warning("Non-numeric abundance values transformed into 1. Using presence/absence scale.")
    dt$Abund[dt$Abund != 0] <- 1
    abund <- "pa"
  }
  dt$Abund <- as.numeric(dt$Abund)
  
  if (is.null(names(cluster))) {
    if (length(cluster) != length(unique(dt$Sample))) {
      stop("Cluster vector length must equal number of unique samples.")
    }
    names(cluster) <- unique(dt$Sample)
  }
  cluster_vec <- factor(cluster)
  group <- levels(cluster_vec)
  dt[, cluster := factor(cluster_vec[as.character(Sample)], levels = group)]
  if (any(is.na(dt$cluster))) {
    stop("Cluster vector does not cover all samples.")
  }
  samplesize <- table(cluster_vec)
  
  if (!abund %in% c("percentage", "pa")) {
    stop("Argument 'abund' must be either percentages ('percentage') or presence/absence data ('pa')")
  }
  
  if (type == "totalfreq") {
    freq_dt <- dt[Abund > 0, .(freq = data.table::uniqueN(Sample)), by = .(TaxonName, cluster)]
    syntab <- data.table::dcast(freq_dt, TaxonName ~ cluster, value.var = "freq", fill = 0)
    syntab_df <- as.data.frame(syntab)
    rownames(syntab_df) <- syntab_df$TaxonName
    syntab_df$TaxonName <- NULL
    syntab_df <- syntab_df[, group, drop = FALSE]
    results <- list("syntable" = syntab_df, "samplesize" = samplesize)
    
  } else if (type == "percfreq") {
    freq_dt <- dt[Abund > 0, .(freq = data.table::uniqueN(Sample)), by = .(TaxonName, cluster)]
    syntab <- data.table::dcast(freq_dt, TaxonName ~ cluster, value.var = "freq", fill = 0)
    syntab_df <- as.data.frame(syntab)
    rownames(syntab_df) <- syntab_df$TaxonName
    syntab_df$TaxonName <- NULL
    for (g in group) {
      syntab_df[[g]] <- round(syntab_df[[g]] * 100 / samplesize[g], digits = digits)
    }
    results <- list("syntable" = syntab_df, "samplesize" = samplesize)
    
  } else if (type == "mean" && abund == "percentage") {
    sum_dt <- dt[, .(sum = sum(Abund)), by = .(TaxonName, cluster)]
    syntab <- data.table::dcast(sum_dt, TaxonName ~ cluster, value.var = "sum", fill = 0)
    syntab_df <- as.data.frame(syntab)
    rownames(syntab_df) <- syntab_df$TaxonName
    syntab_df$TaxonName <- NULL
    for (g in group) {
      syntab_df[[g]] <- round(syntab_df[[g]] / samplesize[g], digits = digits)
    }
    results <- list("syntable" = syntab_df, "samplesize" = samplesize)
    
  } else if (type == "mean" && abund == "pa") {
    stop("Cannot calculate mean abundance in clusters with presence/absence values.")
    
  } else if (type == "median" && abund == "percentage") {
    species <- sort(unique(dt$TaxonName))
    ids <- unique(dt$Sample)
    complete <- data.table::CJ(Sample = ids, TaxonName = species)
    complete <- merge(complete, dt[, .(Sample, TaxonName, Abund)],
                      by = c("Sample", "TaxonName"), all.x = TRUE)
    complete[, cluster := factor(cluster_vec[as.character(Sample)], levels = group)]
    complete[is.na(Abund), Abund := 0]
    med_dt <- complete[, .(value = stats::median(Abund)), by = .(TaxonName, cluster)]
    syntab <- data.table::dcast(med_dt, TaxonName ~ cluster, value.var = "value", fill = 0)
    syntab_df <- as.data.frame(syntab)
    rownames(syntab_df) <- syntab_df$TaxonName
    syntab_df$TaxonName <- NULL
    syntab_df <- round(syntab_df, digits = digits)
    results <- list("syntable" = syntab_df, "samplesize" = samplesize)
    
  } else if (type == "median" && abund == "pa") {
    stop("Cannot calculate median abundance in clusters with presence/absence values.")
    
  } else if (type == "phi") {
    freq_dt <- dt[Abund > 0, .(freq = data.table::uniqueN(Sample)), by = .(TaxonName, cluster)]
    freq <- data.table::dcast(freq_dt, TaxonName ~ cluster, value.var = "freq", fill = 0)
    freq_df <- as.data.frame(freq)
    rownames(freq_df) <- freq_df$TaxonName
    freq_df$TaxonName <- NULL
    freq_df <- freq_df[, group, drop = FALSE]
    N <- sum(samplesize)
    n <- rowSums(freq_df)
    phitab <- freq_df
    for (i in seq_along(group)) {
      phitab[, i] <- (N * freq_df[, i] - n * samplesize[i]) /
        sqrt(n * samplesize[i] * (N - n) * (N - samplesize[i]))
    }
    results <- list("syntable" = phitab, "samplesize" = samplesize)
    
  } else {
    stop("Cannot calculate synoptic table. Define correct type of species matrix values to use (abund = c('percentage', 'pa')).\nCheck correct type of synoptic table output type (type = c('totalfreq', 'percfreq', 'mean', 'median', 'phi')).")
  }
  
  return(invisible(results))
}