#' Synoptic tables and calculation of group-wise frequencies, fidelity and
#' differential species character (wide species-sample matrix)
#'
#' @keywords internal
#' @noRd
syntable_wide <- function(matrix, groups, abund = "percentage", type = "percfreq", digits = 0,
                          phi_method = "default", phi_transform = "none",
                          phi_standard = "none", phi_reps = 100) {
  
  # --- Input & cleaning -------------------------------------------------------
  X <- as.matrix(matrix)
  if (ncol(X) < 1 || nrow(X) < 1)
    stop("Input 'matrix' must have >= 1 row (samples) and >= 1 column (species).")
  if (is.null(colnames(X)))
    stop("Column names (species) must be set.")
  
  # Align 'groups' to sample order (row names optional)
  if (any(is.na(groups))) stop("NA values in 'groups' are not allowed.")
  if (!is.null(names(groups))) {
    if (!is.null(rownames(X))) {
      # align by sample names
      if (!all(rownames(X) %in% names(groups)))
        stop("Named 'groups' must include all sample (row) names of 'matrix'.")
      groups <- groups[rownames(X)]
    } else {
      # no row names: fall back to position
      if (length(groups) != nrow(X))
        stop("Length of 'groups' must equal number of rows (samples).")
      groups <- unname(groups)
      message("Row names are missing; ignoring names(groups) and aligning by row order.")
    }
  } else {
    if (length(groups) != nrow(X))
      stop("Length of 'groups' must equal the number of rows (samples).")
  }
  grp_factor <- factor(groups)
  grp_levels  <- levels(grp_factor)
  if (length(grp_levels) < 1) stop("No grp_levels in 'groups'.")
  
  # Handle empty strings and NAs; robust PA conversion
  Xc <- matrix(as.character(X), nrow = nrow(X), ncol = ncol(X),
               dimnames = dimnames(X))
  if (any(Xc == "", na.rm = TRUE)) Xc[Xc == ""] <- NA
  Xnum <- suppressWarnings(matrix(as.numeric(Xc), nrow = nrow(X), ncol = ncol(X),
                                  dimnames = dimnames(X)))
  nonnum <- is.na(Xnum) & !is.na(Xc)
  
  if (any(nonnum)) {
    pres <- matrix(FALSE, nrow = nrow(X), ncol = ncol(X), dimnames = dimnames(X))
    pres[!is.na(Xnum)] <- Xnum[!is.na(Xnum)] != 0
    pres[nonnum] <- Xc[nonnum] != "0"
    X <- matrix(0, nrow = nrow(X), ncol = ncol(X), dimnames = dimnames(X))
    X[pres] <- 1
    abund <- "pa"
    message("Non-numeric abundance values detected -> using presence/absence scale.")
  } else {
    Xnum[is.na(Xnum)] <- 0
    X <- Xnum
  }
  
  if (!abund %in% c("percentage", "pa"))
    stop("Argument 'abund' must be either 'percentage' (cover) or 'pa' (presence/absence).")
  
  # Basic sizes and group indicator
  N <- nrow(X); p <- ncol(X)
  Z <- stats::model.matrix(~ grp_factor - 1)
  colnames(Z) <- grp_levels
  b_vec <- colSums(Z)                    # group sizes
  
  # Helper: finalize p x G numeric matrix -> data.frame with species rows
  .finalize <- function(M_pg) {
    df <- as.data.frame(M_pg)
    rownames(df) <- colnames(X)          # species
    colnames(df) <- grp_levels
    df
  }
  
  # --- Calculations by 'type' -------------------------------------------------
  if (type == "totalfreq") {
    P <- X > 0
    counts_Gp <- crossprod(Z, P)         # G x p
    results <- list("syntable" = .finalize(t(counts_Gp)),
                    "samplesize" = b_vec)
    
  } else if (type == "percfreq") {
    P <- X > 0
    counts_Gp <- crossprod(Z, P)
    percfreq <- t(counts_Gp)                           # p x G
    percfreq <- sweep(percfreq, 2, b_vec, "/", check.margin = FALSE) * 100
    percfreq[, b_vec == 0] <- NA                 # avoid div/0 artefacts
    percfreq <- round(percfreq, digits = digits)
    results <- list("syntable" = .finalize(percfreq),
                    "samplesize" = b_vec)
    
  } else if (type == "mean" && abund == "percentage") {
    sums_Gp <- crossprod(Z, X)
    means   <- t(sweep(sums_Gp, 1, b_vec, "/", check.margin = FALSE))
    means[, b_vec == 0] <- NA
    means   <- round(means, digits = digits)
    results <- list("syntable" = .finalize(means),
                    "samplesize" = b_vec)
    
  } else if (type == "mean" && abund == "pa") {
    stop("Cannot calculate mean abundance in groups with presence/absence values.")
    
  } else if (type == "median" && abund == "percentage") {
    med <- matrix(NA, nrow = p, ncol = length(grp_levels))
    for (j in seq_along(grp_levels)) {
      idx <- which(grp_factor == grp_levels[j]) 
      med[, j] <- if (length(idx)) apply(X[idx, , drop = FALSE], 2, stats::median) else NA
    }
    med <- round(med, digits = digits)
    results <- list("syntable" = .finalize(med),
                    "samplesize" = b_vec)
    
  } else if (type == "median" && abund == "pa") {
    stop("Cannot calculate median abundance in groups with presence/absence values.")
    
  } else if (type == "diffspec") {
    # needs >= 3 grp_levels and no empty grp_levels
    if (length(grp_levels) < 3)
      stop("'diffspec' requires at least 3 grp_levels.")
    if (any(b_vec == 0))
      stop("At least one group has zero samples; cannot run 'diffspec' (division by zero).")
    
    P <- X > 0
    counts_Gp <- crossprod(Z, P)               # G x p
    syn <- t(counts_Gp)                        # p x G
    syn <- sweep(syn, 2, b_vec, "/", check.margin = FALSE) * 100
    syn <- round(syn, digits = digits)
    
    pb <- utils::txtProgressBar(min = 1, max = 23, style = 3)
    
    small <- matrix(NA, nrow = nrow(syn), ncol = ncol(syn), dimnames = dimnames(syn))
    for (i in seq_len(nrow(syn))) small[i, ] <- sort(as.numeric(syn[i, ]))
    utils::setTxtProgressBar(pb, 3)
    
    Avsmall <- matrix(NA, nrow = nrow(syn), ncol = ncol(syn))
    Avsmall[, 1] <- small[, 1]
    for (i in seq_len(nrow(syn))) for (k in 2:ncol(syn)) {
      Avsmall[i, k] <- mean(as.numeric(small[i, 1:k]))
    }
    utils::setTxtProgressBar(pb, 4)
    
    CondAvsmall <- matrix("", nrow = nrow(syn), ncol = ncol(syn))
    for (k in seq_len(nrow(syn))) {
      CondAvsmall[k, 2] <- if (small[k, 2] >= 2 * Avsmall[k, 1] + 20) "1" else "0"
    }
    for (i in 3:ncol(syn)) for (k in seq_len(nrow(syn))) {
      CondAvsmall[k, i] <- if (CondAvsmall[k, i - 1] == "1") "1"
      else if (small[k, i] >= 2 * Avsmall[k, i - 1] + 20) "1" else "0"
    }
    utils::setTxtProgressBar(pb, 5)
    
    Avsmall2 <- matrix(NA, nrow = nrow(syn), ncol = ncol(syn))
    Avsmall2[, 1] <- Avsmall[, 1]
    for (i in 2:ncol(syn)) for (k in seq_len(nrow(syn))) {
      Avsmall2[k, i] <- if (CondAvsmall[k, i] == "1") Avsmall2[k, i - 1] else Avsmall[k, i]
    }
    utils::setTxtProgressBar(pb, 6)
    
    large <- matrix(NA, nrow = nrow(syn), ncol = ncol(syn), dimnames = dimnames(syn))
    for (i in seq_len(nrow(syn))) large[i, ] <- sort(as.numeric(syn[i, ]), decreasing = TRUE)
    utils::setTxtProgressBar(pb, 7)
    
    Avlarge <- matrix(NA, nrow = nrow(syn), ncol = ncol(syn))
    Avlarge[, 1] <- large[, 1]
    for (i in seq_len(nrow(syn))) for (k in 2:ncol(syn)) {
      Avlarge[i, k] <- mean(as.numeric(large[i, 1:k]))
    }
    utils::setTxtProgressBar(pb, 8)
    
    CondAvlarge <- matrix("", nrow = nrow(syn), ncol = ncol(syn))
    for (k in seq_len(nrow(syn))) {
      CondAvlarge[k, 2] <- if (large[k, 2] <= (Avlarge[k, 1] / 2) - 10) "1" else "0"
    }
    for (i in 3:ncol(syn)) for (k in seq_len(nrow(syn))) {
      CondAvlarge[k, i] <- if (CondAvlarge[k, i - 1] == "1") "1"
      else if (large[k, i] <= (Avlarge[k, i - 1] / 2) - 10) "1" else "0"
    }
    utils::setTxtProgressBar(pb, 9)
    
    Avlarge2 <- matrix(NA, nrow = nrow(syn), ncol = ncol(syn))
    Avlarge2[, 1] <- Avlarge[, 1]
    for (i in 2:ncol(syn)) for (k in seq_len(nrow(syn))) {
      Avlarge2[k, i] <- if (CondAvlarge[k, i] == "1") Avlarge2[k, i - 1] else Avlarge[k, i]
    }
    utils::setTxtProgressBar(pb, 10)
    
    CondAvsmall2 <- matrix("", nrow = nrow(syn), ncol = ncol(syn))
    for (i in 2:ncol(syn)) for (k in seq_len(nrow(syn))) {
      CondAvsmall2[k, i] <- if (small[k, i] > (Avlarge2[k, i - 1] / 2) - 10) "1" else "0"
    }
    utils::setTxtProgressBar(pb, 11)
    
    Avsmall3 <- matrix(NA, nrow = nrow(syn), ncol = ncol(syn))
    Avsmall3[, 1] <- Avsmall2[, 1]
    for (i in 2:ncol(syn)) for (k in seq_len(nrow(syn))) {
      Avsmall3[k, i] <- if (CondAvsmall2[k, i] == "1") Avsmall3[k, i - 1] else Avsmall2[k, i]
    }
    utils::setTxtProgressBar(pb, 12)
    
    CondAvlarge2 <- matrix("", nrow = nrow(syn), ncol = ncol(syn))
    for (i in 2:ncol(syn)) for (k in seq_len(nrow(syn))) {
      CondAvlarge2[k, i] <- if (large[k, i] < 2 * Avsmall2[k, i - 1] + 20) "1" else "0"
    }
    utils::setTxtProgressBar(pb, 13)
    
    Avlarge3 <- matrix(NA, nrow = nrow(syn), ncol = ncol(syn))
    Avlarge3[, 1] <- Avlarge2[, 1]
    for (i in 2:ncol(syn)) for (k in seq_len(nrow(syn))) {
      Avlarge3[k, i] <- if (CondAvlarge2[k, i] == "1") Avlarge3[k, i - 1] else Avlarge2[k, i]
    }
    utils::setTxtProgressBar(pb, 14)
    
    Small_Avsmall3 <- matrix("", nrow = nrow(syn), ncol = ncol(syn))
    for (i in seq_len(ncol(syn))) for (k in seq_len(nrow(syn))) {
      if (is.na(syn[k, 1])) {
        Small_Avsmall3[k, i] <- ""
      } else if (syn[k, i] >= 2 * max(Avsmall3[k, ], na.rm = TRUE) + 20) {
        Small_Avsmall3[k, i] <- "p"
      } else if (syn[k, i] < 2 * max(Avsmall3[k, ], na.rm = TRUE) + 20) {
        Small_Avsmall3[k, i] <- "n"
      } else stop("Calculation of differential species failed (Small_Avsmall3).")
    }
    utils::setTxtProgressBar(pb, 15)
    
    Large_AvLarge3 <- matrix("", nrow = nrow(syn), ncol = ncol(syn))
    for (i in seq_len(ncol(syn))) for (k in seq_len(nrow(syn))) {
      if (is.na(syn[k, 1])) {
        Large_AvLarge3[k, i] <- ""
      } else if (syn[k, i] >  min(Avlarge3[k, ], na.rm = TRUE) / 2 - 10) {
        Large_AvLarge3[k, i] <- "p"
      } else if (syn[k, i] <= min(Avlarge3[k, ], na.rm = TRUE) / 2 - 10) {
        Large_AvLarge3[k, i] <- "n"
      } else stop("Calculation of differential species failed (Large_AvLarge3).")
    }
    utils::setTxtProgressBar(pb, 16)
    
    Sum_cond <- matrix("", nrow = nrow(syn), ncol = ncol(syn))
    for (i in seq_len(ncol(syn))) for (k in seq_len(nrow(syn))) {
      Sum_cond[k, i] <- paste0(Small_Avsmall3[k, i], Large_AvLarge3[k, i])
    }
    utils::setTxtProgressBar(pb, 17)
    
    Results <- matrix("-", nrow = nrow(syn), ncol = ncol(syn), dimnames = dimnames(syn))
    for (i in seq_len(ncol(syn))) for (k in seq_len(nrow(syn))) {
      Results[k, i] <- switch(Sum_cond[k, i],
                              "pp" = "p", "nn" = "n", "pn" = "pn", "-")
    }
    utils::setTxtProgressBar(pb, 18)
    
    subresult <- matrix(0L, nrow = nrow(syn), ncol = ncol(syn))
    subresult[Results %in% c("n", "p", "pn")] <- 1L
    onlydiff <- Results[rowSums(subresult) > 0, , drop = FALSE]
    others   <- rownames(Results)[rowSums(subresult) == 0]
    
    # Per-group lists
    positivespec <- matrix(NA, nrow = length(grp_levels), ncol = 200)
    negativespec <- matrix(NA, nrow = length(grp_levels), ncol = 200)
    positivenegativespec <- matrix(NA, nrow = length(grp_levels), ncol = 200)
    for (i in seq_along(grp_levels)) {
      sp <- rownames(Results)
      pos <- sp[Results[, i] == "p"]; neg <- sp[Results[, i] == "n"]; pn <- sp[Results[, i] == "pn"]
      if (length(pos)) positivespec[i, seq_along(pos)] <- pos
      if (length(neg)) negativespec[i, seq_along(neg)] <- neg
      if (length(pn)) positivenegativespec[i, seq_along(pn)] <- pn
    }
    
    diffspeclist <- setNames(vector("list", length(grp_levels)), grp_levels)
    for (i in seq_along(grp_levels)) {
      pos <- positivespec[i, ][!is.na(positivespec[i, ])]
      neg <- negativespec[i, ][!is.na(negativespec[i, ])]
      pn  <- positivenegativespec[i, ][!is.na(positivenegativespec[i, ])]
      diffspeclist[[i]] <- list("positive diff"          = if (length(pos)) pos else noquote("no positive diff species"),
                                "negative diff"          = if (length(neg)) neg else noquote("no negative diff species"),
                                "positive/negative diff" = if (length(pn))  pn  else noquote("no positive/negative diff species"))
    }
    
    utils::setTxtProgressBar(pb, 23); close(pb)
    
    results <- list("syntable"   = .finalize(Results),
                    "onlydiff"   = as.data.frame(onlydiff),
                    "others"     = others,
                    "samplesize" = b_vec,
                    "differentials" = diffspeclist)
    
  } else if (type == "phi") {
    # ----- Fidelity -----------------------------------------------------------
    if (phi_standard == "rarefy") {
      m <- min(b_vec)
      if (m < 1) stop("At least one group has zero samples; cannot rarefy.")
      phisum <- matrix(0, nrow = p, ncol = length(grp_levels))
      rownames(phisum) <- colnames(X); colnames(phisum) <- grp_levels
      for (r in seq_len(phi_reps)) {
        idx <- unlist(lapply(grp_levels, function(g) sample(which(grp_factor == g), m)))
        sub <- X[idx, , drop = FALSE]
        sub_cl <- grp_factor[idx]   
        tmp <- syntable_wide(sub, sub_cl, abund = abund, type = "phi",
                             digits = digits, phi_method = phi_method,
                             phi_transform = phi_transform,
                             phi_standard = "none", phi_reps = phi_reps)
        phisum <- phisum + as.matrix(tmp$syntable)
      }
      phitab <- phisum / phi_reps
      samplesize <- rep(m, length(grp_levels)); names(samplesize) <- grp_levels
      results <- list("syntable" = as.data.frame(phitab),
                      "samplesize" = samplesize)
      
    } else {
      if (phi_method %in% c("default", "ochiai", "uvalue")) {
        P <- X > 0
        a_Gp <- crossprod(Z, P)                 # G x p
        a_pg <- t(a_Gp)                         # p x G
        n_vec <- colSums(P)                     # species totals (length p)
        b <- b_vec                              # group sizes (length G)
        
        if (phi_method == "default") {
          num <- N * a_pg - n_vec %o% b
          den <- sqrt((n_vec * (N - n_vec)) %o% (b * (N - b)))
        } else if (phi_method == "ochiai") {
          num <- a_pg
          den <- sqrt(n_vec %o% b)
        } else { # uvalue
          num <- a_pg - (n_vec %o% (b / N))
          den <- sqrt((n_vec * (N - n_vec)) %o% (b * (N - b) / (N * (N - 1))))
        }
        
        phitab <- num / den
        phitab[!is.finite(phitab)] <- 0
        rownames(phitab) <- colnames(X); colnames(phitab) <- grp_levels
        
        if (phi_standard == "adjust") {
          if (phi_method %in% c("default", "uvalue")) {
            m <- min(b)
            for (j in seq_along(b)) {
              bj <- b[j]
              if (is.finite(bj) && bj > 0 && bj < N) {
                phitab[, j] <- phitab[, j] * sqrt(m * (N - m) / (bj * (N - bj)))
              } else {
                phitab[, j] <- 0
              }
            }
          } else {
            message("phi_standard = 'adjust' is not defined for phi_method = '",
                    phi_method, "'. Returning unadjusted values. ",
                    "Use phi_standard = 'rarefy' if you need size standardisation.")
          }
        }
        
        results <- list("syntable" = as.data.frame(phitab),
                        "samplesize" = b)
        
      } else if (phi_method == "cover") {
        if (abund != "percentage")
          stop("Cover-based fidelity requires percentage abundances.")
        C <- switch(phi_transform,
                    none = X,
                    sqrt = sqrt(X),
                    log  = log(X + 1),
                    stop("Unknown phi_transform"))
        sx_Gp <- crossprod(Z, C)                 # G x p
        sx_pg <- t(sx_Gp)                        # p x G
        sxt   <- colSums(C)
        sx2t  <- colSums(C^2)
        b <- b_vec
        
        num <- N * sx_pg - sxt %o% b
        den <- sqrt((N * sx2t - sxt^2) %o% (b * (N - b)))
        phitab <- num / den
        phitab[!is.finite(phitab)] <- 0
        rownames(phitab) <- colnames(X); colnames(phitab) <- grp_levels
        
        if (phi_standard == "adjust") {
          message("phi_standard = 'adjust' is not defined for phi_method = 'cover'. ",
                  "Returning unadjusted values. Use 'rarefy' for size standardisation.")
        }
        
        results <- list("syntable" = as.data.frame(phitab),
                        "samplesize" = b)
      } else {
        stop("Unknown phi_method. Use 'default', 'cover', 'ochiai', or 'uvalue'.")
      }
    }
  } else {
    stop("Cannot calculate synoptic table. Check 'abund' (percentage|pa) and 'type' (totalfreq|percfreq|mean|median|diffspec|phi).")
  }
  
  return(invisible(results))
}
