#' Heterogeneity-constrained random resampling (HCR)
#' @description
#' Performs heterogeneity-constrained random (HCR) resampling (Lengyel, Chytrý & Tichý, 2011) of community data.
#' Within each stratum (e.g., grid cell), many random subsets of plots are evaluated and the subset with
#' the highest mean dissimilarity and the lowest variance of dissimilarities is retained. Optionally, the
#' number of plots per stratum is adapted from the stratum’s mean pairwise dissimilarity (\eqn{\beta}-diversity).
#' 
#'
#' @param data_wide a data-frame like object with the following column contents: 
#'   \itemize{
#'     \item column 1: sample ids
#'     \item column 2: strata
#'     \item columns 3...n: species.
#'   }
#' @param transform One of \code{c("none","sqrt","log1p","binary")}. If "binary", values become 0/1 and
#'   \code{vegan::vegdist(binary = TRUE)} is used.
#' @param score_dist Dissimilarity method for trial scoring; any method accepted by
#'   \code{vegan::vegdist} (e.g., "bray","jaccard", "hellinger", "euclidean", "canberra",
#'   "gower", "kulczynski","morisita","horn","mountford","raup","binomial",
#'   "chao","cao", …).
#' @param beta_dist One of \code{c("bray","jaccard")} for per-stratum mean dissimilarity used to calculate the 
#' adaptive number of plots.
#' With \code{transform="binary"}, "bray" equals Sørensen.
#' @param adaptive_n Logical. If TRUE, adapt the number of plots per stratum from
#'   \code{beta_mean * max_plots} bounded to \code{[min_plots, max_plots]}; if FALSE, use fixed \code{n_plots}.
#' @param n_plots Fixed number of plots per stratum when \code{adaptive_n=FALSE}. If \code{NA},
#'   defaults to \code{max_plots} (capped at stratum size).
#' @param min_plots,max_plots Global default min/max number of plots per stratum
#' @param min_stratum_n Minimum stratum size under which the whole stratum is selected (no resampling).
#' @param trials Number of random trials per stratum (default 1000).
#' @param group_vec Optional vector (length \code{nrow(data_wide)}) assigning each sample to a higher-level
#'   group (e.g., country, region). Used only if \code{adaptive_n=TRUE}.
#' @param group_limits Optional \code{data.frame} with group-specific limits. The first column
#'   must contain group names; it must also contain numeric columns named \code{"min_plots"} and \code{"max_plots"}.
#'   Other columns are ignored.
#' @param write_csv Optional file path to write a CSV with columns \code{sample_id, selected}. If \code{NULL}, no file.
#' @param progress Show a text progress bar (default: \code{interactive()}).
#' @param seed Optional integer seed for reproducibility of random subset trials.
#'   
#' @section Details:
#' The algorithm follows Lengyel, Chytrý & Tichý (2011) and was based upon the JUICE implementation (Tichý, 2002). 
#' For speed, it precomputes per-stratum distance matrices (once) and reuses them across trials, which
#' enables large numbers of trials (default \code{trials = 1000}). 
#' 
#' Within each stratum candidate subsets are scored using \code{score_dist} by high mean dissimilarity and low variance of dissimilarities.
#' 
#' If \code{adaptive_n = TRUE} (default), the target number of plots is computed as a linear function of the mean pairwise 
#' dissimilarity (\eqn{\beta}-diversity; \code{beta_dist}) and the maximum number of plots (\code{beta_mean * max_plots}; Wiser & de Cáceres, 2013) and then
#' bounded to \code{[min_plots, max_plots]} and the stratum size. 
#' 
#' Additionally group-specific limits for minimum and maximum numbers of plots per stratum can be supplied via 
#' \code{group_vec} and \code{group_limits}. Each sample is assigned to a higher-level group 
#' (e.g., country or region), and the minimum and maximum number of plots are defined per group. 
#' This allows, for example, larger plot limits to be set for larger countries or regions.
#' 
#' @return A \code{data.frame} with \code{sample_id} and \code{selected} (0/1).
#'   Attributes: \code{selected_rows} (logical) and \code{params}.
#'   
#' @author Friedemann von Lampe
#'
#' @references 
#' Lengyel, A., Chytrý, M., & Tichý, L. (2011). Heterogeneity-constrained random resampling of phytosociological databases.
#' \emph{Journal of Vegetation Science}, \strong{22(1)}, 175–183. \doi{10.1111/j.1654-1103.2010.01225.x}
#' 
#' Tichý, L. (2002). JUICE, software for vegetation classification.  \emph{Journal of Vegetation Science}, \strong{13(3)}, 451. 
#' \doi{10.1658/1100-9233(2002)013[0451:JSFVC]2.0.CO;2}
#' 
#' Wiser, S. K., & de Cáceres, M. (2013). Updating vegetation classifications: an example with New Zealand's woody vegetation. 
#' \emph{Journal of Vegetation Science}, \strong{24(1)}, 80–93.  \doi{10.1111/j.1654-1103.2012.01450.x}
#' 
#' @export
#' @importFrom stats var


hcr_resampling <- function(
    data_wide,
    transform     = c("none","sqrt","log1p","binary"),
    score_dist    = "bray",
    beta_dist     = c("bray","jaccard"),
    adaptive_n    = TRUE,
    n_plots       = NA_integer_,
    min_plots     = 10L,
    max_plots     = 100L,
    min_stratum_n = 10L,
    trials        = 1000L,
    group_vec     = NULL,
    group_limits  = NULL,
    write_csv     = NULL,
    progress      = interactive(),
    seed          = NULL
) {
  if (!requireNamespace("vegan", quietly = TRUE)) {
    stop("Package 'vegan' is required (for vegdist). Please install it.")
  }
  transform <- match.arg(transform)
  beta_dist <- match.arg(beta_dist)
  
  ## validate inputs
  if (ncol(data_wide) < 3L)
    stop("data_wide must have at least 3 columns: 1=id, 2=strata, 3..=species.")
  if (!is.numeric(trials) || length(trials) != 1L || is.na(trials) || trials < 1)
    stop("'trials' must be a single number >= 1.")
  trials <- as.integer(trials)
  if (!is.numeric(min_plots) || !is.numeric(max_plots))
    stop("'min_plots' and 'max_plots' must be numeric.")
  if (length(min_plots) != 1L || length(max_plots) != 1L ||
      is.na(min_plots) || is.na(max_plots) || min_plots < 1 || max_plots < min_plots)
    stop("Provide sensible 'min_plots' and 'max_plots' (min >= 1, max >= min).")
  if (!is.null(n_plots) && !is.na(n_plots) &&
      (!is.numeric(n_plots) || length(n_plots) != 1L || n_plots < 1))
    stop("'n_plots' must be a single number >= 1 (or NA).")
  if (!is.numeric(min_stratum_n) || length(min_stratum_n) != 1L || is.na(min_stratum_n) || min_stratum_n < 1)
    stop("'min_stratum_n' must be a single number >= 1.")
  if (!is.null(group_vec) && length(group_vec) != nrow(data_wide))
    stop("If provided, 'group_vec' must have length equal to nrow(data_wide).")
  if (!is.null(group_limits) && !all(c("min_plots","max_plots") %in% names(group_limits)))
    stop("'group_limits' must contain columns named 'min_plots' and 'max_plots' (first column = group names).")
  
  if (!is.null(seed)) set.seed(seed)
  
  # fix positions: 1 = id, 2 = strata
  dataset  <- data_wide[, c(1L, 2L, setdiff(seq_len(ncol(data_wide)), 1:2)), drop = FALSE]
  dataset1 <- dataset[, -c(1L, 2L), drop = FALSE]
  
  # coerce species to numeric
  non_num <- !vapply(dataset1, is.numeric, logical(1))
  if (any(non_num)) {
    suppressWarnings({
      dataset1[non_num] <- lapply(dataset1[non_num], function(z) as.numeric(as.character(z)))
    })
  }
  dataset1 <- as.matrix(dataset1)
  if (!is.numeric(dataset1)) stop("Species columns could not be coerced to numeric.")
  
  ## NA handling
  if (anyNA(dataset1)) {
    dataset1[is.na(dataset1)] <- 0
    message("NAs replaced by 0")
  }
  
  # transformations
  bin_flag <- FALSE
  if (transform == "binary") {
    dataset1 <- ifelse(dataset1 > 0, 1, 0)
    bin_flag <- TRUE
  } else if (transform == "sqrt") {
    dataset1 <- sqrt(pmax(dataset1, 0))
  } else if (transform == "log1p") {
    dataset1 <- log1p(pmax(dataset1, 0))
  }
  
  n <- nrow(dataset)
  strata   <- as.factor(dataset[[2L]])
  levs     <- levels(strata)
  selected <- rep(FALSE, n)
  
  # prepare group_limits lookup (first column = group names)
  have_group_limits <- !is.null(group_limits) &&
    ncol(group_limits) >= 3L &&
    all(c("min_plots","max_plots") %in% names(group_limits))
  if (have_group_limits) {
    gl_key       <- as.character(group_limits[[1L]])
    gl_min_plots <- suppressWarnings(as.numeric(group_limits[["min_plots"]]))
    gl_max_plots <- suppressWarnings(as.numeric(group_limits[["max_plots"]]))
  }
  
  # check if group limit definitions are provided
  if (!is.null(group_vec)) {
    if (!have_group_limits) {
      message("group_limits not provided; using global min_plots/max_plots for all groups.")
    } else {
      used_groups     <- unique(as.character(group_vec[!is.na(group_vec)]))
      defined_groups  <- unique(gl_key)  # first column of group_limits
      missing_groups  <- setdiff(used_groups, defined_groups)
      
      if (length(missing_groups) > 0) {
        message("Groups not defined in 'group_limits' (using global min/max): ",
                paste(sort(missing_groups), collapse = ", "))
      }
      if (any(is.na(group_vec))) {
        message("Some samples have NA in 'group_vec'; they will use global min/max.")
      }
    }
  }
  
  
  # initialise progress bar
  pb <- NULL
  if (isTRUE(progress)) {
    pb <- utils::txtProgressBar(min = 0, max = length(levs), style = 3)
    on.exit({ if (!is.null(pb)) close(pb) }, add = TRUE)
  }
  
  # helper: pick min/max for a stratum from group_limits
  pick_limits_for_stratum <- function(idx) {
    stratum_min_plots <- min_plots
    stratum_max_plots <- max_plots
    if (is.null(group_vec) || !have_group_limits) {
      return(list(min_plots = stratum_min_plots, max_plots = stratum_max_plots))
    }
    grp <- group_vec[idx][1]
    if (is.na(grp)) return(list(min_plots = stratum_min_plots, max_plots = stratum_max_plots))
    hit <- which(gl_key == as.character(grp))[1]
    if (length(hit) == 1 && !is.na(hit)) {
      if (!is.na(gl_min_plots[hit])) stratum_min_plots <- as.integer(gl_min_plots[hit])
      if (!is.na(gl_max_plots[hit])) stratum_max_plots <- as.integer(gl_max_plots[hit])
    }
    list(min_plots = stratum_min_plots, max_plots = stratum_max_plots)
  }
  
  ## loop BEGIN: iterating over all strata 
  for (k in seq_along(levs)) {
    lev <- levs[k]
    idx <- which(strata == lev)
    m   <- length(idx)
    
    # group-specific min/max (fallback to globals)
    lims  <- pick_limits_for_stratum(idx)
    stratum_min_plots <- lims$min_plots
    stratum_max_plots <- lims$max_plots
    
    ## precompute per-stratum distance matrices (much faster)
    # for beta_mean (can differ from score_dist)
    beta_mean <- 0
    if (m > 1L) {
      if (identical(beta_dist, score_dist)) {
        # compute D_score below and reuse it, set flag
        reuse_beta_from_score <- TRUE
      } else {
        reuse_beta_from_score <- FALSE
        D_beta <- vegan::vegdist(dataset1[idx, , drop = FALSE],
                                 method = beta_dist, binary = bin_flag, na.rm = TRUE)
        beta_mean <- mean(as.numeric(D_beta))
      }
    } else {
      reuse_beta_from_score <- FALSE
    }
    
    # precompute score distance matrix for this stratum
    if (m > 1L) {
      D_score <- vegan::vegdist(dataset1[idx, , drop = FALSE],
                                method = score_dist, binary = bin_flag, na.rm = TRUE)
      Dm <- as.matrix(D_score)
      if (reuse_beta_from_score) {
        beta_mean <- mean(D_score)  # mean of 'dist' object equals mean of upper triangle
      }
    }
    
    # decide plots for this stratum
    if (isTRUE(adaptive_n)) {
      # Compute the desired number of plots using a linear function
      # linear function where the slope is set to stratum_max_plots and the intercept is zero.
      # when beta_mean is zero, the result n_plots_est is 0
      n_plots_est   <- beta_mean * stratum_max_plots
      # ensure minimum and maximum thresholds
      target_n_plots <- as.integer(max(stratum_min_plots, min(stratum_max_plots, round(n_plots_est))))
      target_n_plots <- max( min(target_n_plots, m), min(2L, m) )
      min_stratum_size  <- target_n_plots
    } else {
      target_n_plots <- if (is.na(n_plots)) min(max_plots, m) else as.integer(n_plots)
      target_n_plots <- max( min(target_n_plots, m), min(2L, m) )
      min_stratum_size  <- as.integer(min_stratum_n)
    }
    
    # stratum smaller than minimum value: keep all; else run trials and score
    # trial_subsets: matrix that stores all random candidate subsets
    if (m < min_stratum_size  || target_n_plots < 2L) {
      selected[idx] <- TRUE
    } else {
      subset_mean_diss    <- numeric(trials)
      subset_var_diss <- numeric(trials)
      trial_subsets <- matrix(0L, nrow = trials, ncol = target_n_plots)
      
      # for each trial, it randomly samples target_n_plots relevés from those in the stratum 
      # and stores their indices in trial_subsets
      for (j in seq_len(trials)) {
        cand_local <- sample.int(m, size = target_n_plots, replace = FALSE)
        trial_subsets[j, ] <- cand_local
        # extract upper-triangle distances for the candidate subset
        dsub <- Dm[cand_local, cand_local, drop = FALSE]
        diss <- dsub[upper.tri(dsub, diag = FALSE)]
        subset_mean_diss[j]    <- mean(diss)
        subset_var_diss[j] <- if (length(diss) > 1L) stats::var(diss) else 0
        if (is.na(subset_var_diss[j])) subset_var_diss[j] <- 0
      }
      rank_mean_diss <- rank(subset_mean_diss, ties.method = "average")     # higher mean dissimilarity (heterogeneity) is better
      rank_var_diss <- rank(-subset_var_diss, ties.method = "average")     # lower variance is better
      subset_score    <- rank_mean_diss + rank_var_diss
      best_trial   <- which.max(subset_score)
      selected[idx[trial_subsets[best_trial, ]]] <- TRUE
    }
    
    if (!is.null(pb)) utils::setTxtProgressBar(pb, k)
  }
  ## Loop END 
  
  # create a matrix marker with two columns. 
  # first column is taken from the first column of dataset (identifier)
  # second column is set to 1 for all rows that were marked as selected.
  out <- data.frame(
    sample_id = dataset[[1L]],
    selected  = as.integer(selected),
    stringsAsFactors = FALSE
  )
  
  # optional file output
  if (!is.null(write_csv)) {
    dir.create(dirname(write_csv), showWarnings = FALSE, recursive = TRUE)
    utils::write.csv(out, file = write_csv, row.names = FALSE)
  }
  
  attr(out, "selected_rows") <- selected
  attr(out, "params") <- list(
    transform = transform, score_dist = score_dist, beta_dist = beta_dist,
    adaptive_n = adaptive_n, trials = trials,
    min_plots = min_plots, max_plots = max_plots, min_stratum_n = min_stratum_n,
    n_plots = n_plots, seed = seed
  )
  out
}
