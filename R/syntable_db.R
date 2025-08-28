#' Synoptic tables for long-format vegetation databases
#'
#' @description
#' `syntable_db` calculates synoptic tables for long-format vegetation data
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
#' @export
syntable_db <- function(db, cluster, abund = "percentage", type = "percfreq", digits = 0) {
  if (!all(c("RELEVE_NR", "TaxonName", "Cover_Perc") %in% names(db))) {
    stop("Input table must contain columns 'RELEVE_NR', 'TaxonName' and 'Cover_Perc'.")
  }

  narep <- FALSE
  if (any(db$Cover_Perc == "", na.rm = TRUE)) {
    db$Cover_Perc[db$Cover_Perc == ""] <- NA
    narep <- TRUE
  }
  if (any(is.na(db$Cover_Perc))) {
    db$Cover_Perc[is.na(db$Cover_Perc)] <- 0
    narep <- TRUE
  }
  if (narep) {
    print("NA and/or empty character values replaced by 0.")
  }

  if (any(is.na(cluster))) {
    stop("NA values in cluster not allowed.")
  }

  if (any(is.na(suppressWarnings(as.numeric(db$Cover_Perc))))) {
    warning("Non-numeric cover values transformed into 1. Using presence/absence scale.")
    db$Cover_Perc[db$Cover_Perc != 0] <- 1
    abund <- "pa"
  }
  db$Cover_Perc <- as.numeric(db$Cover_Perc)

  if (is.null(names(cluster))) {
    if (length(cluster) != length(unique(db$RELEVE_NR))) {
      stop("Cluster vector length must equal number of unique RELEVE_NR.")
    }
    names(cluster) <- unique(db$RELEVE_NR)
  }

  cluster <- factor(cluster)
  group <- levels(cluster)
  db$cluster <- factor(cluster[as.character(db$RELEVE_NR)], levels = group)
  if (any(is.na(db$cluster))) {
    stop("Cluster vector does not cover all RELEVE_NR.")
  }

  samplesize <- tapply(rep(1, length(cluster)), cluster, sum)

  if (!abund %in% c("percentage", "pa")) {
    stop("Argument 'abund' must be either percentages ('percentage') or presence/absence data ('pa')")
  }

  if (type == "totalfreq") {
    syntab <- xtabs((Cover_Perc > 0) ~ TaxonName + cluster, data = db)
    syntab <- as.data.frame.matrix(syntab)
    syntab <- syntab[, group, drop = FALSE]
    results <- list("syntable" = syntab, "samplesize" = samplesize)

  } else if (type == "percfreq") {
    freq <- xtabs((Cover_Perc > 0) ~ TaxonName + cluster, data = db)
    freq <- as.data.frame.matrix(freq)[, group, drop = FALSE]
    syntab <- sweep(freq, 2, samplesize, function(x, s) round(x * 100 / s, digits = digits))
    results <- list("syntable" = syntab, "samplesize" = samplesize)

  } else if (type == "mean" && abund == "percentage") {
    sumtab <- xtabs(Cover_Perc ~ TaxonName + cluster, data = db)
    sumtab <- as.data.frame.matrix(sumtab)[, group, drop = FALSE]
    syntab <- sweep(sumtab, 2, samplesize, function(x, s) round(x / s, digits = digits))
    results <- list("syntable" = syntab, "samplesize" = samplesize)

  } else if (type == "mean" && abund == "pa") {
    stop("Cannot calculate mean cover in clusters with presence/absence values.")

  } else if (type == "median" && abund == "percentage") {
    species <- sort(unique(db$TaxonName))
    syntab <- matrix(NA, nrow = length(species), ncol = length(group),
                     dimnames = list(species, group))
    for (g in group) {
      sampids <- names(cluster)[cluster == g]
      n_samp <- length(sampids)
      sub <- db[db$cluster == g, c("RELEVE_NR", "TaxonName", "Cover_Perc")]
      for (sp in species) {
        vals <- sub$Cover_Perc[sub$TaxonName == sp]
        if (length(vals) < n_samp) {
          vals <- c(vals, rep(0, n_samp - length(vals)))
        }
        syntab[sp, g] <- round(stats::median(vals), digits = digits)
      }
    }
    results <- list("syntable" = as.data.frame(syntab), "samplesize" = samplesize)

  } else if (type == "median" && abund == "pa") {
    stop("Cannot calculate median cover in clusters with presence/absence values.")

  } else if (type == "phi") {
    freq <- xtabs((Cover_Perc > 0) ~ TaxonName + cluster, data = db)
    freq <- as.data.frame.matrix(freq)[, group, drop = FALSE]
    N <- sum(samplesize)
    n <- rowSums(freq)
    phitab <- freq
    for (i in seq_along(group)) {
      phitab[, i] <- (N * freq[, i] - n * samplesize[i]) /
        sqrt(n * samplesize[i] * (N - n) * (N - samplesize[i]))
    }
    results <- list("syntable" = phitab, "samplesize" = samplesize)

  } else {
    stop("Cannot calculate synoptic table. Define correct type of species matrix values to use (abund = c('percentage', 'pa')). \nCheck correct type of synoptic table output type (type = c('totalfreq', 'percfreq', 'mean', 'median', 'phi')).")
  }

  return(invisible(results))
}
