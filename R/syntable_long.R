#' Synoptic tables and calculation of cluster-wise frequencies and fidelity for
#' long-format vegetation databases
#' 
#' @import data.table
#' @keywords internal
#' @noRd


syntable_long <- function(vegdata, cluster, abund = "percentage", type = "percfreq", digits = 0,
                          phi_method = "default", phi_transform = "none",
                          phi_standard = "none", phi_reps = 100) {
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
    species <- sort(unique(dt$TaxonName))
    N_full <- length(cluster_vec)
    samplesize_full <- table(cluster_vec)
    
    calc_phi <- function(dt_sub) {
      N <- length(unique(dt_sub$Sample))
      samplesize <- table(dt_sub$cluster)
      
      if (phi_method %in% c("default", "ochiai", "uvalue")) {
        freq_dt <- dt_sub[Abund > 0,
                          .(freq = data.table::uniqueN(Sample)),
                          by = .(TaxonName, cluster)]
        freq_mat <- data.table::dcast(freq_dt, TaxonName ~ cluster,
                                      value.var = "freq", fill = 0)
        freq_mat <- merge(data.table::data.table(TaxonName = species),
                          freq_mat, by = "TaxonName", all.x = TRUE)
        freq_mat[is.na(freq_mat)] <- 0
        a_mat <- as.matrix(freq_mat[, -1, with = FALSE])
        rownames(a_mat) <- freq_mat$TaxonName
        n_vec <- rowSums(a_mat)
        phitab <- matrix(0, nrow = length(species), ncol = length(group),
                         dimnames = list(species, group))
        for (i in seq_along(group)) {
          b <- samplesize[group[i]]
          for (k in seq_along(species)) {
            a <- a_mat[k, i]
            n <- n_vec[k]
            if (phi_method == "default") {
              den <- sqrt(n * b * (N - n) * (N - b))
              phitab[k, i] <- if (den == 0) 0 else (N * a - n * b) / den
            } else if (phi_method == "ochiai") {
              den <- sqrt(b * n)
              phitab[k, i] <- if (den == 0) 0 else a / den
            } else if (phi_method == "uvalue") {
              den <- sqrt(n * b * (N - b) * (N - n) / (N * (N - 1)))
              phitab[k, i] <- if (den == 0) 0 else (a - n * b / N) / den
            }
          }
        }
        return(list(phitab = phitab, samplesize = samplesize, N = N))
        
      } else if (phi_method == "cover") {
        if (abund != "percentage") {
          stop("Cover-based fidelity requires percentage abundances.")
        }
        dt_sub[, cov := switch(phi_transform,
                               none = Abund,
                               sqrt = sqrt(Abund),
                               log = log(Abund + 1),
                               stop("Unknown phi_transform"))]
        total_stats <- dt_sub[, .(sum_x = sum(cov), sum_x2 = sum(cov^2)),
                              by = TaxonName]
        group_stats <- dt_sub[, .(sum_x = sum(cov)),
                              by = .(TaxonName, cluster)]
        total_stats <- merge(data.table::data.table(TaxonName = species),
                             total_stats, by = "TaxonName", all.x = TRUE)
        group_stats <- data.table::dcast(group_stats, TaxonName ~ cluster,
                                         value.var = "sum_x", fill = 0)
        group_stats <- merge(data.table::data.table(TaxonName = species),
                             group_stats, by = "TaxonName", all.x = TRUE)
        total_stats[is.na(total_stats)] <- 0
        group_stats[is.na(group_stats)] <- 0
        phitab <- matrix(0, nrow = length(species), ncol = length(group),
                         dimnames = list(species, group))
        for (i in seq_along(group)) {
          b <- samplesize[group[i]]
          sum_x_group <- group_stats[[group[i]]]
          sum_x_total <- total_stats$sum_x
          sum_x2_total <- total_stats$sum_x2
          num <- N * sum_x_group - sum_x_total * b
          den <- sqrt((N * sum_x2_total - sum_x_total^2) * b * (N - b))
          phitab[, i] <- ifelse(den == 0, 0, num / den)
        }
        return(list(phitab = phitab, samplesize = samplesize, N = N))
        
      } else {
        stop("Unknown phi_method")
      }
    }
    
    if (phi_standard == "rarefy") {
      m <- min(samplesize_full)
      phisum <- matrix(0, nrow = length(species), ncol = length(group),
                       dimnames = list(species, group))
      for (r in seq_len(phi_reps)) {
        ids <- unlist(lapply(group, function(g) {
          sample(unique(dt[cluster == g, Sample]), m)
        }))
        sub_dt <- dt[Sample %in% ids]
        res <- calc_phi(sub_dt)
        phisum <- phisum + res$phitab
      }
      phitab <- phisum / phi_reps
      samplesize_out <- rep(m, length(group))
      names(samplesize_out) <- group
    } else {
      res <- calc_phi(dt)
      phitab <- res$phitab
      samplesize_out <- samplesize_full
      if (phi_standard == "adjust" && phi_method %in% c("default", "uvalue")) {
        m <- min(samplesize_out)
        for (i in seq_along(group)) {
          b <- samplesize_out[i]
          phitab[, i] <- phitab[, i] *
            sqrt(m * (N_full - m) / (b * (N_full - b)))
        }
      }
    }
    
    results <- list("syntable" = as.data.frame(phitab),
                    "samplesize" = samplesize_out)
    
  } else {
    stop("Cannot calculate synoptic table. Define correct type of species matrix values to use (abund = c('percentage', 'pa')).\nCheck correct type of synoptic table output type (type = c('totalfreq', 'percfreq', 'mean', 'median', 'phi')).")
  }
  
  return(invisible(results))
}