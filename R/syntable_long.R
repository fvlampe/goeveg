#' Synoptic tables and calculation of cluster-wise frequencies and fidelity for
#' long-format vegetation databases
#' 
#' @import data.table
#' @keywords internal
#' @noRd


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