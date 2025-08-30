#' Synoptic tables and calculation of group-wise frequencies and fidelity for
#' long-format vegetation databases
#'
#' @import data.table
#' @keywords internal
#' @noRd
syntable_long <- function(vegdata, groups, abund = "percentage", type = "percfreq", digits = 0,
                          phi_method = "default", phi_transform = "none",
                          phi_standard = "none", phi_reps = 100) {

  # avoid R CMD check notes for data.table
  Sample <- TaxonName <- Abund <- groups_dt <- cov <- sum_x <- sum_x2 <- value <- NULL
  
  # --- Input & cleaning -------------------------------------------------------
  
  dt <- data.table::as.data.table(vegdata)
  if (ncol(dt) < 3) stop("Input table must contain at least three columns (sample, taxon, abundance).")
  data.table::setnames(dt, 1:3, c("Sample", "TaxonName", "Abund"))
  
  # empties -> NA
  if (any(dt$Abund == "", na.rm = TRUE)) dt$Abund[dt$Abund == ""] <- NA
  
  # numeric/PA handling
  Abund_char <- as.character(dt$Abund)
  Abund_num  <- suppressWarnings(as.numeric(Abund_char))
  nonnum     <- is.na(Abund_num) & !is.na(Abund_char)
  
  if (any(nonnum)) {
    # Any non-zero numeric OR any non-numeric string != "0" => presence
    pres <- rep(FALSE, nrow(dt))
    pres[!is.na(Abund_num)] <- Abund_num[!is.na(Abund_num)] != 0
    pres[nonnum] <- Abund_char[nonnum] != "0"
    pres[is.na(Abund_char)] <- FALSE
    dt[, Abund := as.numeric(pres)]
    abund <- "pa"
    message("Non-numeric abundance values detected -> using presence/absence scale.")
  } else {
    if (any(is.na(Abund_num))) {
      Abund_num[is.na(Abund_num)] <- 0
      message("NA values in abundances replaced by 0.")
    }
    dt[, Abund := Abund_num]
  }
  
  # Validate groups vector
  if (any(is.na(groups))) stop("NA values in 'groups' are not allowed.")
  
  # Ensure name-based mapping Sample -> groups
  if (is.null(names(groups))) {
    samp_order <- unique(as.character(dt$Sample))
    if (length(groups) != length(samp_order)) {
      stop("groups vector length must equal number of unique samples if unnamed.")
    }
    names(groups) <- samp_order
  }
  grp_factor <- factor(groups)
  grp_levels <- levels(grp_factor)
  
  # Attach groups to rows
  dt[, groups_dt := factor(grp_factor[as.character(Sample)], levels = grp_levels)]
  if (any(is.na(dt$groups_dt))) stop("groups vector does not cover all samples present in 'vegdata'.")
  
  # Samplesize per group (counts SAMPLES, not records)
  samplesize_vec <- table(grp_factor)  # ordered like 'grp_levels'
  
  if (!abund %in% c("percentage", "pa")) {
    stop("Argument 'abund' must be either 'percentage' (cover) or 'pa' (presence/absence).")
  }
  
  # Helper to coerce output with consistent column order
  .finalize_syntab <- function(syntab_dt) {
    syntab_df <- as.data.frame(syntab_dt)
    rownames(syntab_df) <- syntab_df$TaxonName
    syntab_df$TaxonName <- NULL
    syntab_df <- syntab_df[, grp_levels, drop = FALSE]
    syntab_df
  }
  
  # --- Calculations by 'type' -------------------------------------------------
  
  if (type == "totalfreq") {
    
    # % frequency counts of samples with presence
    dt_pa <- unique(dt[Abund > 0, .(TaxonName, groups = groups_dt, Sample)])
    freq_dt <- dt_pa[, .N, by = .(TaxonName, groups)]
    data.table::setnames(freq_dt, "N", "freq")
    syntab <- data.table::dcast(freq_dt, TaxonName ~ groups, value.var = "freq", fill = 0)
    results <- list("syntable" = .finalize_syntab(syntab),
                    "samplesize" = samplesize_vec)
    
  } else if (type == "percfreq") {
    
    dt_pa <- unique(dt[Abund > 0, .(TaxonName, groups = groups_dt, Sample)])
    freq_dt <- dt_pa[, .N, by = .(TaxonName, groups)]
    data.table::setnames(freq_dt, "N", "freq")
    syntab <- data.table::dcast(freq_dt, TaxonName ~ groups, value.var = "freq", fill = 0)
    syntab_df <- .finalize_syntab(syntab)
    for (g in grp_levels) {
      syntab_df[[g]] <- round(syntab_df[[g]] * 100 / as.numeric(samplesize_vec[g]), digits = digits)
    }
    results <- list("syntable" = syntab_df,
                    "samplesize" = samplesize_vec)
    
  } else if (type == "mean" && abund == "percentage") {
    
    sum_dt <- dt[, .(sum = sum(Abund)), by = .(TaxonName, groups = groups_dt)]
    syntab <- data.table::dcast(sum_dt, TaxonName ~ groups, value.var = "sum", fill = 0)
    syntab_df <- .finalize_syntab(syntab)
    for (g in grp_levels) {
      syntab_df[[g]] <- round(syntab_df[[g]] / as.numeric(samplesize_vec[g]), digits = digits)
    }
    results <- list("syntable" = syntab_df,
                    "samplesize" = samplesize_vec)
    
  } else if (type == "mean" && abund == "pa") {
    
    stop("Cannot calculate mean abundance in groups with presence/absence values.")
    
  } else if (type == "median" && abund == "percentage") {
    
    species <- sort(unique(dt$TaxonName))
    ids <- unique(dt$Sample)
    complete <- data.table::CJ(Sample = ids, TaxonName = species)
    complete <- merge(complete, dt[, .(Sample, TaxonName, Abund)],
                      by = c("Sample", "TaxonName"), all.x = TRUE)
    complete[, groups_dt := factor(grp_factor[as.character(Sample)], levels = grp_levels)]
    complete[is.na(Abund), Abund := 0]
    med_dt <- complete[, .(value = stats::median(Abund)), by = .(TaxonName, groups = groups_dt)]
    syntab <- data.table::dcast(med_dt, TaxonName ~ groups, value.var = "value", fill = 0)
    syntab_df <- .finalize_syntab(syntab)
    syntab_df <- round(syntab_df, digits = digits)
    results <- list("syntable" = syntab_df,
                    "samplesize" = samplesize_vec)
    
  } else if (type == "median" && abund == "pa") {
    
    stop("Cannot calculate median abundance in groups with presence/absence values.")
    
  } else if (type == "phi") {
    
    # ----- Fidelity -----------------------------------------------------------
    
    species <- sort(unique(dt$TaxonName))
    
    calc_phi <- function(dt_sub, samp_ids = NULL) {
      # Determine N and samplesize from SAMPLES
      if (is.null(samp_ids)) {
        samp_unique <- unique(dt_sub[, .(Sample, groups = groups_dt)])
        N <- data.table::uniqueN(samp_unique$Sample)
        samplesize <- table(factor(samp_unique$groups, levels = grp_levels))
      } else {
        N <- length(unique(samp_ids))
        cl_map <- grp_factor[as.character(samp_ids)]
        samplesize <- table(factor(cl_map, levels = grp_levels))
      }
      if (N <= 1) {
        phitab0 <- matrix(0, nrow = length(species), ncol = length(grp_levels),
                          dimnames = list(species, grp_levels))
        return(list(phitab = phitab0, samplesize = samplesize, N = N))
      }
      
      if (phi_method %in% c("default", "ochiai", "uvalue")) {
        # Build presence counts per (Taxon, group)
        dt_pa <- unique(dt_sub[Abund > 0, .(TaxonName, groups = groups_dt, Sample)])
        a_dt  <- dt_pa[, .N, by = .(TaxonName, groups)]
        freq_mat <- data.table::dcast(a_dt, TaxonName ~ groups, value.var = "N", fill = 0)
        freq_mat <- merge(data.table::data.table(TaxonName = species),
                          freq_mat, by = "TaxonName", all.x = TRUE)
        freq_mat[is.na(freq_mat)] <- 0
        
        a_mat <- as.matrix(freq_mat[, -1, with = FALSE])       # species x grp_levels
        rownames(a_mat) <- freq_mat$TaxonName
        n_vec <- rowSums(a_mat)                                # species totals
        b_vec <- as.numeric(samplesize)                         # group sizes
        Nv    <- as.numeric(N)
        
        if (phi_method == "default") {
          num <- Nv * a_mat - n_vec %o% b_vec
          den <- sqrt((n_vec * (Nv - n_vec)) %o% (b_vec * (Nv - b_vec)))
        } else if (phi_method == "ochiai") {
          num <- a_mat
          den <- sqrt(n_vec %o% b_vec)
        } else { # uvalue
          num <- a_mat - (n_vec %o% (b_vec / Nv))
          den <- sqrt((n_vec * (Nv - n_vec)) %o% (b_vec * (Nv - b_vec) / (Nv * (Nv - 1))))
        }
        
        phitab <- num / den
        phitab[!is.finite(phitab)] <- 0
        dimnames(phitab) <- list(rownames(a_mat), colnames(a_mat))
        return(list(phitab = phitab, samplesize = samplesize, N = N))
        
      } else if (phi_method == "cover") {
        if (abund != "percentage") stop("Cover-based fidelity requires percentage abundances.")
        dt_sub[, cov := switch(phi_transform,
                               none = Abund,
                               sqrt = sqrt(Abund),
                               log  = log(Abund + 1),
                               stop("Unknown phi_transform"))]
        
        total_stats <- dt_sub[, .(sum_x = sum(cov), sum_x2 = sum(cov^2)), by = TaxonName]
        total_stats <- merge(data.table::data.table(TaxonName = species),
                             total_stats, by = "TaxonName", all.x = TRUE)
        total_stats[is.na(total_stats)] <- 0
        
        group_stats <- dt_sub[, .(sum_x = sum(cov)), by = .(TaxonName, groups = groups_dt)]
        group_stats <- data.table::dcast(group_stats, TaxonName ~ groups, value.var = "sum_x", fill = 0)
        group_stats <- merge(data.table::data.table(TaxonName = species),
                             group_stats, by = "TaxonName", all.x = TRUE)
        group_stats[is.na(group_stats)] <- 0
        
        sx   <- as.matrix(group_stats[, -1, with = FALSE])     # species x grp_levels
        rownames(sx) <- group_stats$TaxonName
        sxt  <- total_stats$sum_x
        sx2t <- total_stats$sum_x2
        b_vec <- as.numeric(samplesize)
        Nv <- as.numeric(N)
        
        num <- Nv * sx - sxt %o% b_vec
        den <- sqrt((Nv * sx2t - sxt^2) %o% (b_vec * (Nv - b_vec)))
        phitab <- num / den
        phitab[!is.finite(phitab)] <- 0
        dimnames(phitab) <- list(rownames(sx), colnames(sx))
        return(list(phitab = phitab, samplesize = samplesize, N = N))
        
      } else {
        stop("Unknown phi_method. Use 'default', 'cover', 'ochiai', or 'uvalue'.")
      }
    }
    
    # Sample IDs per group from the groups VECTOR
    ids_by_group <- lapply(grp_levels, function(g) names(grp_factor)[grp_factor == g])
    
    if (phi_standard == "rarefy") {
      m <- min(lengths(ids_by_group))
      if (m < 1) stop("At least one group has zero samples; cannot rarefy.")
      phisum <- matrix(0, nrow = length(species), ncol = length(grp_levels),
                       dimnames = list(species, grp_levels))
      for (r in seq_len(phi_reps)) {
        samp_ids <- unlist(lapply(ids_by_group, function(v) sample(v, m)))
        sub_dt <- dt[Sample %in% samp_ids]
        res <- calc_phi(sub_dt, samp_ids = samp_ids)
        phisum <- phisum + res$phitab
      }
      phitab <- phisum / phi_reps
      samplesize_out <- rep(m, length(grp_levels))
      names(samplesize_out) <- grp_levels
      
    } else {
      # No rarefaction: compute directly on full data
      res <- calc_phi(dt)
      phitab <- res$phitab
      samplesize_out <- res$samplesize
      
      if (phi_standard == "adjust") {
        if (phi_method %in% c("default", "uvalue")) {
          # analytical adjustment only for these methods
          N_adj <- res$N
          m <- min(as.numeric(samplesize_out))
          for (i in seq_along(grp_levels)) {
            b <- as.numeric(samplesize_out[i])
            if (is.finite(b) && b > 0 && b < N_adj) {
              phitab[, i] <- phitab[, i] * sqrt(m * (N_adj - m) / (b * (N_adj - b)))
            } else {
              phitab[, i] <- 0
            }
          }
        } else {
          message("phi_standard = 'adjust' is not defined for phi_method = '",
                  phi_method, "'. Returning unadjusted values. ",
                  "Use phi_standard = 'rarefy' if you need size standardisation.")
        }
      } else if (phi_standard != "none" && phi_standard != "rarefy") {
        stop("Unknown phi_standard. Use 'none', 'rarefy', or 'adjust'.")
      }
      
    }
    
    results <- list("syntable" = as.data.frame(phitab),
                    "samplesize" = samplesize_out)
    
  } else {
    stop("Cannot calculate synoptic table. Check 'abund' (percentage|pa) and 'type' (totalfreq|percfreq|mean|median|phi).")
  }
  
  return(invisible(results))
}
