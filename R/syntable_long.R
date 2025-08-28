#' Synoptic tables for long-format vegetation databases
#'
#' @description
#' `syntable_long` calculates synoptic tables for long-format vegetation data
#' where species occur in rows and samples are given by `RELEVE_NR`. 
#'
#' @param db Data frame in long format containing at least the columns
#'   `RELEVE_NR`, `TaxonName` and `Cover_Perc`.
#' @param cluster Integer or character vector/factor with classification
#'   cluster identity for each `RELEVE_NR`.
#' @param abund Type of abundances. Define whether `Cover_Perc` contains
#'   percentage cover (`abund = "percentage"`, default) or presence/absence
#'   data (`abund = "pa"`, with values 0/1).
#' @param type Type of synoptic table output
#'   `type = c("percfreq", "totalfreq", "mean", "median", "phi")`.
#' @param digits Integer indicating the number of decimal places to be
#'   displayed in result tables (default 0).
#'
#' @return An (invisible) list with components:
#'   \item{\code{$syntable }}{unordered synoptic table for given species and clusters}
#'   \item{\code{$samplesize }}{total number of samples per cluster}
#'
#' @seealso [syntable()], [synsort()]
#' @import data.table
#' @export
syntable_long <- function(db, cluster, abund = "percentage", type = "percfreq", digits = 0) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required for syntable_long.")
  }
  dt <- data.table::as.data.table(db)
  if (!all(c("RELEVE_NR", "TaxonName", "Cover_Perc") %in% names(dt))) {
    stop("Input table must contain columns 'RELEVE_NR', 'TaxonName' and 'Cover_Perc'.")
  }
  
  narep <- FALSE
  if (any(dt$Cover_Perc == "", na.rm = TRUE)) {
    dt$Cover_Perc[dt$Cover_Perc == ""] <- NA
    narep <- TRUE
  }
  if (any(is.na(dt$Cover_Perc))) {
    dt$Cover_Perc[is.na(dt$Cover_Perc)] <- 0
    narep <- TRUE
  }
  if (narep) {
    print("NA and/or empty character values replaced by 0.")
  }
  
  if (any(is.na(cluster))) {
    stop("NA values in cluster not allowed.")
  }
  
  if (any(is.na(suppressWarnings(as.numeric(dt$Cover_Perc))))) {
    warning("Non-numeric cover values transformed into 1. Using presence/absence scale.")
    dt$Cover_Perc[dt$Cover_Perc != 0] <- 1
    abund <- "pa"
  }
  dt$Cover_Perc <- as.numeric(dt$Cover_Perc)
  
  if (is.null(names(cluster))) {
    if (length(cluster) != length(unique(dt$RELEVE_NR))) {
      stop("Cluster vector length must equal number of unique RELEVE_NR.")
    }
    names(cluster) <- unique(dt$RELEVE_NR)
  }
  cluster_vec <- factor(cluster)
  group <- levels(cluster_vec)
  dt[, cluster := factor(cluster_vec[as.character(RELEVE_NR)], levels = group)]
  if (any(is.na(dt$cluster))) {
    stop("Cluster vector does not cover all RELEVE_NR.")
  }
  samplesize <- table(cluster_vec)
  
  if (!abund %in% c("percentage", "pa")) {
    stop("Argument 'abund' must be either percentages ('percentage') or presence/absence data ('pa')")
  }
  
  if (type == "totalfreq") {
    freq_dt <- dt[Cover_Perc > 0, .(freq = data.table::uniqueN(RELEVE_NR)), by = .(TaxonName, cluster)]
    syntab <- data.table::dcast(freq_dt, TaxonName ~ cluster, value.var = "freq", fill = 0)
    syntab_df <- as.data.frame(syntab)
    rownames(syntab_df) <- syntab_df$TaxonName
    syntab_df$TaxonName <- NULL
    syntab_df <- syntab_df[, group, drop = FALSE]
    results <- list("syntable" = syntab_df, "samplesize" = samplesize)
    
  } else if (type == "percfreq") {
    freq_dt <- dt[Cover_Perc > 0, .(freq = data.table::uniqueN(RELEVE_NR)), by = .(TaxonName, cluster)]
    syntab <- data.table::dcast(freq_dt, TaxonName ~ cluster, value.var = "freq", fill = 0)
    syntab_df <- as.data.frame(syntab)
    rownames(syntab_df) <- syntab_df$TaxonName
    syntab_df$TaxonName <- NULL
    for (g in group) {
      syntab_df[[g]] <- round(syntab_df[[g]] * 100 / samplesize[g], digits = digits)
    }
    results <- list("syntable" = syntab_df, "samplesize" = samplesize)
    
  } else if (type == "mean" && abund == "percentage") {
    sum_dt <- dt[, .(sum = sum(Cover_Perc)), by = .(TaxonName, cluster)]
    syntab <- data.table::dcast(sum_dt, TaxonName ~ cluster, value.var = "sum", fill = 0)
    syntab_df <- as.data.frame(syntab)
    rownames(syntab_df) <- syntab_df$TaxonName
    syntab_df$TaxonName <- NULL
    for (g in group) {
      syntab_df[[g]] <- round(syntab_df[[g]] / samplesize[g], digits = digits)
    }
    results <- list("syntable" = syntab_df, "samplesize" = samplesize)
    
  } else if (type == "mean" && abund == "pa") {
    stop("Cannot calculate mean cover in clusters with presence/absence values.")
    
  } else if (type == "median" && abund == "percentage") {
    species <- sort(unique(dt$TaxonName))
    ids <- unique(dt$RELEVE_NR)
    complete <- data.table::CJ(RELEVE_NR = ids, TaxonName = species)
    complete <- merge(complete, dt[, .(RELEVE_NR, TaxonName, Cover_Perc)],
                      by = c("RELEVE_NR", "TaxonName"), all.x = TRUE)
    complete[, cluster := factor(cluster_vec[as.character(RELEVE_NR)], levels = group)]
    complete[is.na(Cover_Perc), Cover_Perc := 0]
    med_dt <- complete[, .(value = stats::median(Cover_Perc)), by = .(TaxonName, cluster)]
    syntab <- data.table::dcast(med_dt, TaxonName ~ cluster, value.var = "value", fill = 0)
    syntab_df <- as.data.frame(syntab)
    rownames(syntab_df) <- syntab_df$TaxonName
    syntab_df$TaxonName <- NULL
    syntab_df <- round(syntab_df, digits = digits)
    results <- list("syntable" = syntab_df, "samplesize" = samplesize)
    
  } else if (type == "median" && abund == "pa") {
    stop("Cannot calculate median cover in clusters with presence/absence values.")
    
  } else if (type == "phi") {
    freq_dt <- dt[Cover_Perc > 0, .(freq = data.table::uniqueN(RELEVE_NR)), by = .(TaxonName, cluster)]
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