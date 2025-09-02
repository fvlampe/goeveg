#' Synoptic tables and calculation of group-wise frequencies and fidelity for
#' long-format vegetation databases
#'
#' @import data.table
#' @keywords internal
#' @noRd
syntable_long <- function(vegdata,
                          groups,
                          abund = "percentage",
                          type = "percfreq",
                          digits = NULL,
                          phi_method = "default",
                          phi_standard = c("none", "target", "all"),
                          phi_target_size = NULL,  # % of N for the target group; default = 100/G
                          phi_alpha = NULL,
                          group_col = NULL) {
  
  phi_standard <- match.arg(phi_standard)
  
  # avoid R CMD check notes for data.table
  Sample <- TaxonName <- Abund <- groups_dt <- cov <- sum_x <- sum_x2 <- value <- grp <- n <- NULL 
  
  # Rounding: 3 for phi when 'digits' not supplied; 0 otherwise
  if (is.null(digits)) {
    digits_nonphi <- 0L
    digits_phi    <- 3L
  } else {
    d <- as.integer(digits)
    digits_nonphi <- d
    digits_phi    <- d
  }
  
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
    # presence if numeric non-zero OR non-numeric string != "0"
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
  
  if (!abund %in% c("percentage", "pa")) {
    stop("Argument 'abund' must be either 'percentage' (cover) or 'pa' (presence/absence).")
  }
  
  # --- Groups mapping (accept column, per-row vector, or per-sample vector) ---
  # Result of this block:
  #   dt$groups_dt : factor with the group of each row
  #   sample_map   : per-sample table with a single group per Sample (factor)
  #   grp_levels   : group level names
  if (!is.null(group_col)) {
    if (!group_col %in% names(dt))
      stop("group_col '", group_col, "' not found in vegdata.")
    if (any(is.na(dt[[group_col]])))
      stop("NA values in group_col are not allowed.")
    dt[, groups_dt := factor(.SD[[1]]), .SDcols = group_col]
    
  } else if (!missing(groups) && length(groups) == nrow(dt)) {
    if (any(is.na(groups))) stop("NA values in 'groups' are not allowed.")
    dt[, groups_dt := factor(groups)]
    
  } else if (!missing(groups) && !is.null(names(groups))) {
    # named per-sample vector
    if (!all(unique(dt$Sample) %in% names(groups)))
      stop("Named 'groups' must include all sample IDs in column 1 (Sample).")
    dt[, groups_dt := factor(groups[as.character(Sample)])]
    if (any(is.na(dt$groups_dt)))
      stop("'groups' vector does not cover all samples present in 'vegdata'.")
    
  } else if (!missing(groups) && length(groups) == data.table::uniqueN(dt$Sample)) {
    # per-sample vector in the order of first appearance of each Sample
    sample_order <- dt[, .SD[1L], by = Sample][, Sample]
    dt[, groups_dt := factor(groups[match(Sample, sample_order)])]
    if (any(is.na(dt$groups_dt)))
      stop("'groups' vector (per-sample ordered) could not be mapped to all samples.")
    
  } else {
    stop("Provide 'groups' via group_col=, or as a vector: ",
         "per-row (length nrow(vegdata)) or per-sample (named or length = number of unique samples).")
  }
  
  # Validate: each Sample must have exactly one group (no within-sample conflicts)
  consistency <- dt[, .(n = data.table::uniqueN(groups_dt)), by = Sample]
  if (any(consistency$n != 1L)) {
    bad <- consistency[n != 1L, Sample]
    stop("Group assignment is not consistent within samples. ",
         "These samples have >1 group: ", paste(head(bad, 10), collapse = ", "),
         if (length(bad) > 10) " ...")
  }
  
  # Build per-sample map and canonical levels
  sample_map <- dt[, .(grp = groups_dt[1L]), by = Sample]
  sample_map[, grp := factor(grp)]    # drop to observed levels
  grp_levels  <- levels(sample_map$grp)
  dt[, groups_dt := factor(groups_dt, levels = grp_levels)]
  
  # Samplesize per group (counts SAMPLES, not rows)
  samplesize_vec <- table(sample_map$grp)  # ordered like grp_levels
  
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
    
    # absolute frequency = number of plots with presence per group
    dt_pa <- unique(dt[Abund > 0, .(TaxonName, groups = groups_dt, Sample)])
    freq_dt <- dt_pa[, .N, by = .(TaxonName, groups)]
    data.table::setnames(freq_dt, "N", "freq")
    syntab <- data.table::dcast(freq_dt, TaxonName ~ groups, value.var = "freq", fill = 0)
    results <- list("syntable"  = .finalize_syntab(syntab),
                    "samplesize" = samplesize_vec)
    
  } else if (type == "percfreq") {
    
    dt_pa <- unique(dt[Abund > 0, .(TaxonName, groups = groups_dt, Sample)])
    freq_dt <- dt_pa[, .N, by = .(TaxonName, groups)]
    data.table::setnames(freq_dt, "N", "freq")
    syntab <- data.table::dcast(freq_dt, TaxonName ~ groups, value.var = "freq", fill = 0)
    syntab_df <- .finalize_syntab(syntab)
    for (g in grp_levels) {
      denom <- as.numeric(samplesize_vec[g])
      syntab_df[[g]] <- if (denom > 0) round(syntab_df[[g]] * 100 / denom, digits = digits_nonphi) else NA
    }
    results <- list("syntable"  = syntab_df,
                    "samplesize" = samplesize_vec)
    
  } else if (type == "mean" && abund == "percentage") {
    
    sum_dt <- dt[, .(sum = sum(Abund)), by = .(TaxonName, groups = groups_dt)]
    syntab <- data.table::dcast(sum_dt, TaxonName ~ groups, value.var = "sum", fill = 0)
    syntab_df <- .finalize_syntab(syntab)
    for (g in grp_levels) {
      denom <- as.numeric(samplesize_vec[g])
      syntab_df[[g]] <- if (denom > 0) round(syntab_df[[g]] / denom, digits = digits_nonphi) else NA
    }
    results <- list("syntable"  = syntab_df,
                    "samplesize" = samplesize_vec)
    
  } else if (type == "mean" && abund == "pa") {
    
    stop("Cannot calculate mean abundance in groups with presence/absence values.")
    
  } else if (type == "median" && abund == "percentage") {
    
    species <- sort(unique(dt$TaxonName))
    ids <- unique(dt$Sample)
    complete <- data.table::CJ(Sample = ids, TaxonName = species)
    complete <- merge(complete, dt[, .(Sample, TaxonName, Abund)],
                      by = c("Sample", "TaxonName"), all.x = TRUE)
    # attach groups for these sample IDs
    gmap <- sample_map[, .(Sample, grp)]
    complete <- merge(complete, gmap, by = "Sample", all.x = TRUE)
    setnames(complete, "grp", "groups_dt")
    complete[is.na(Abund), Abund := 0]
    med_dt <- complete[, .(value = stats::median(Abund)), by = .(TaxonName, groups = groups_dt)]
    syntab <- data.table::dcast(med_dt, TaxonName ~ groups, value.var = "value", fill = 0)
    syntab_df <- round(.finalize_syntab(syntab), digits = digits_nonphi)
    results <- list("syntable"  = syntab_df,
                    "samplesize" = samplesize_vec)
    
  } else if (type == "median" && abund == "pa") {
    
    stop("Cannot calculate median abundance in groups with presence/absence values.")
    
  } else if (type == "phi") {
    
    # ----- Fidelity -----------------------------------------------------------
    species <- sort(unique(dt$TaxonName))
    
    calc_phi <- function(dt_sub, samp_ids = NULL) {
      # Determine N and samplesize from SAMPLES (plots), not rows
      if (is.null(samp_ids)) {
        samp_unique <- unique(dt_sub[, .(Sample, groups = groups_dt)])
        N <- data.table::uniqueN(samp_unique$Sample)
        samplesize <- table(factor(samp_unique$groups, levels = grp_levels))
      } else {
        N <- length(unique(samp_ids))
        cl_map <- sample_map[Sample %in% samp_ids, grp]
        samplesize <- table(factor(cl_map, levels = grp_levels))
      }
      if (N <= 1) {
        phitab0 <- matrix(0, nrow = length(species), ncol = length(grp_levels),
                          dimnames = list(species, grp_levels))
        return(list(phitab = phitab0, samplesize = samplesize, N = N))
      }
      
      if (phi_method %in% c("default", "ochiai", "uvalue")) {
        # presence counts per (Taxon, group)
        dt_pa <- unique(dt_sub[Abund > 0, .(TaxonName, groups = groups_dt, Sample)])
        a_dt  <- dt_pa[, .N, by = .(TaxonName, groups)]
        freq_mat <- data.table::dcast(a_dt, TaxonName ~ groups, value.var = "N", fill = 0)
        freq_mat <- merge(data.table::data.table(TaxonName = species),
                          freq_mat, by = "TaxonName", all.x = TRUE)
        freq_mat[is.na(freq_mat)] <- 0
        
        a_mat <- as.matrix(freq_mat[, -1, with = FALSE])   # species x groups
        rownames(a_mat) <- freq_mat$TaxonName
        n_vec <- rowSums(a_mat)                            # species totals (length = #species)
        b_vec <- as.numeric(samplesize)                    # group sizes          (length = G)
        Nv    <- as.numeric(N)
        
        # one-tailed Fisher test on original counts
        # If phi_alpha is provided, compute p = P[X >= a] with X ~ Hypergeom(N, b, n)
        # using the original (non-equalized) 2x2 counts, and later zero-out cells with p >= alpha.
        do_fisher <- !is.null(phi_alpha)
        if (do_fisher) {
          alpha <- as.numeric(phi_alpha)
          if (!is.finite(alpha) || alpha <= 0 || alpha >= 1)
            stop("phi_alpha must be a number in (0, 1).")
          
          p_mat <- matrix(1, nrow = nrow(a_mat), ncol = ncol(a_mat),
                          dimnames = dimnames(a_mat))
          for (j in seq_along(b_vec)) {
            # one-tailed (greater): P(X >= a_ij) with X ~ Hypergeom(N, m=b_j, k=n_i)
            p_mat[, j] <- stats::phyper(q = a_mat[, j] - 1L,
                                        m = b_vec[j],
                                        n = Nv - b_vec[j],
                                        k = n_vec,
                                        lower.tail = FALSE)
          }
        }
        
        
        # equalization following Tichy & Chytry 2006 
        # phi_standard: "none", "target", "all"
        # phi_target_size: percentage of N for the target group (default = 100/G)
        G <- length(b_vec)
        
        # r_in = a / b  (species × groups), Σ r_in across groups per species
        Rin <- sweep(a_mat, 2, b_vec, "/")
        Rin[!is.finite(Rin)] <- 0
        sumR <- rowSums(Rin)
        
        if (phi_standard == "all") {
          ## Target group gets s_t percent, others equally share the rest.
          s_pct <- if (is.null(phi_target_size)) 100 / G else as.numeric(phi_target_size)
          if (!is.finite(s_pct) || s_pct <= 0 || s_pct >= 100)
            stop("phi_target_size must be a percentage in (0, 100).")
          s_t <- (s_pct / 100) * Nv                       # conceptual size of target group
          s_o <- ((100 - s_pct) / 100) * Nv / (G - 1)     # conceptual size of each non-target
          N_eq <- s_t + (G - 1) * s_o                     # total conceptual plots
          
          # a'_ij and n'_ij for each species i & column j (target = j)
          Aeq_all <- Rin * s_t
          # n'_ij = s_t * r_in_ij + s_o * (Σ r_in_i• - r_in_ij)
          Neq_all <- matrix(s_o * sumR, nrow = nrow(Rin), ncol = G) + (s_t - s_o) * Rin
          
          if (phi_method == "ochiai") {
            den <- sqrt(Neq_all * s_t)
            den[den == 0] <- Inf
            phitab <- Aeq_all / den
          } else {
            num <- N_eq * Aeq_all - Neq_all * s_t
            den <- sqrt(Neq_all * (N_eq - Neq_all) * s_t * (N_eq - s_t))
            den[den == 0] <- Inf
            phi_mat <- num / den
            phitab  <- if (phi_method == "uvalue") phi_mat * sqrt(max(N_eq - 1, 0)) else phi_mat
          }
          
        } else if (phi_standard == "target") {
          ## Equalize ONLY the target group to s_t; pool all other groups into "outside".
          s_pct <- if (is.null(phi_target_size)) 100 / G else as.numeric(phi_target_size)
          if (!is.finite(s_pct) || s_pct <= 0 || s_pct >= 100)
            stop("phi_target_size must be a percentage in (0, 100).")
          s_t <- (s_pct / 100) * Nv
          
          # r_out per target column: (n - a)/(N - b)
          Rout <- sweep(a_mat, 2, b_vec, function(a, b) {
            num <- (n_vec - a); den <- (Nv - b)
            ifelse(den > 0, num / den, 0)
          })
          
          Aeq <- Rin * s_t                            # a'
          Neq <- Aeq + Rout * (Nv - s_t)              # n'
          if (phi_method == "ochiai") {
            den <- sqrt(Neq * s_t); den[den == 0] <- Inf
            phitab <- Aeq / den
          } else {
            num     <- Nv * Aeq - Neq * s_t
            den     <- sqrt(Neq * (Nv - Neq) * s_t * (Nv - s_t))
            den[den == 0] <- Inf
            phi_mat <- num / den
            phitab  <- if (phi_method == "uvalue") phi_mat * sqrt(max(Nv - 1, 0)) else phi_mat
          }
          
        } else {  # phi_standard == "none"
          if (phi_method == "ochiai") {
            phitab <- a_mat / sqrt(n_vec %o% b_vec)
          } else {
            num <- Nv * a_mat - n_vec %o% b_vec
            den <- sqrt((n_vec * (Nv - n_vec)) %o% (b_vec * (Nv - b_vec)))
            den[!is.finite(den) | den == 0] <- Inf
            phi_mat <- num / den
            phitab  <- if (phi_method == "uvalue") phi_mat * sqrt(max(Nv - 1, 0)) else phi_mat
          }
        }
        
        # zero-out non-significant cells if Fisher test requested
        if (do_fisher) {
          phitab[p_mat >= alpha] <- 0
        }
        
        # clean + dimnames, then return from calc_phi()
        phitab[!is.finite(phitab)] <- 0
        dimnames(phitab) <- list(rownames(a_mat), colnames(a_mat))
        return(list(phitab = phitab, samplesize = samplesize, N = N))
        
        
      } else {
        stop("Unknown phi_method. Use 'default', 'ochiai', or 'uvalue'.")
      }
      
    }
    
    # Sample IDs per group (from per-sample map)
    ids_by_group <- lapply(grp_levels, function(g) sample_map[grp == g, Sample])
    
    res <- calc_phi(dt)
    phitab <- res$phitab
    samplesize_out <- res$samplesize
    
    phitab <- round(phitab, digits = digits_phi)
    
    results <- list("syntable"  = as.data.frame(phitab),
                    "samplesize" = samplesize_out)
    
  } else {
    stop("Cannot calculate synoptic table. Check 'abund' (percentage|pa) and 'type' (totalfreq|percfreq|mean|median|phi).")
  }
  
  return(invisible(results))
}
