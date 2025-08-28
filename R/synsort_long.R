#' Sorting functions for synoptic tables (long format)
#'
#' @description
#' `synsort_long` sorts synoptic tables produced by `syntable_long()` using
#' vegetation data supplied in long format. Only the "allspec" method is
#' implemented, omitting differential species analysis.
#'
#' @param syn1 Input synoptic table 1, a data frame with numerical data
#'   obtained from `syntable_long()`. The values of this table will be
#'   displayed in the final output table.
#' @param syn2 Optional second input table with additional numeric sorting
#'   criteria.
#' @param db Community data in long format. The first column must contain
#'   sample identifiers, the second column species names and the third column
#'   numeric values.
#' @param cluster Integer or character vector/factor with classification
#'   cluster identity. Ensure matching order of cluster identity and samples in
#'   `db` for correct allocation of cluster numbers to samples.
#' @param min1 Cluster-wise threshold minimum value for species shown in the
#'   final sorted synoptic table. Species below that minimum will be listed in
#'   the output (`$others` section).
#' @param min2 Threshold minimum value for considering species values of a
#'   numerical second input table `syn2`. Species below that minimum will not be
#'   displayed in the final synoptic table, but will be listed in the output
#'   (`$others` section).
#'
#' @return
#' An (invisible) list analogous to the output of [synsort()].
#'
#' @seealso [syntable_long()], [synsort()]
#'
#' @export
#' @importFrom stats xtabs
synsort_long <- function(syn1, syn2 = syn1, db, cluster, min1 = 0, min2 = 0) {
  if (ncol(db) < 3)
    stop("db must contain at least three columns: sample, species and value")

  wide <- xtabs(db[, 3] ~ db[, 1] + db[, 2])
  matrix <- as.data.frame.matrix(wide)

  synsort(syn1 = syn1, syn2 = syn2, matrix = matrix, cluster = cluster,
          method = "allspec", min1 = min1, min2 = min2)
}
