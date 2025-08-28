#' Sorting functions for synoptic tables
#'
#' @description
#' This function sorts synoptic tables from \code{\link{syntable}} function output. Sorting criteria
#' can be either numerical values in synoptic tables, such as cluster-wise frequencies or fidelity
#' measures, as well as combined criteria that also take into account differential character (according to
#' the criteria defined by Tsiripidis et al., 2009).
#'
#' The algorithm aims to sort species to blocked structure considering the defined criteria and input
#' tables, with the best characterizing species on the top of the block, followed by species with
#' descending importance for plant community description.
#'
#' @param syn1 Input synoptic table 1, a data frame with numerical data format, usually from
#' \code{\link{syntable}} function output. See Details for input table format.
#' The values of this table will be displayed in the final output table.
#' @param syn2 Optional second input table with additional numeric or differential character
#' sorting criteria.
#' @param matrix Optional species-sample matrix, already used for \code{\link{syntable}} function input; used only when calculating the
#'   sorted species-sample matrix (`samples = TRUE`).
#' @param cluster Integer or character vector/factor with classification cluster identity. Ensure matching order of
#' cluster identity and samples in matrix for correct allocation of cluster numbers to samples.
#' @param cluster_order Optional vector giving the desired order of cluster levels. If supplied,
#'   input tables and outputs are rearranged to follow this order.
#' @param cluster_nosort Optional vector of cluster names that should not be used for sorting.
#' @param method Sorting algorithm and synoptic table output options
#'   (\code{method = c("allspec", "alldiff", "totalfreq", "manual")}).
##' @param manual_order Optional character vector of species names (matching
#'   row names of \code{syn1}) to impose a manual row order when
#'   \code{method = "manual"}. Species not listed are appended and ordered by
#'   their overall frequency (row sums of \code{syn1}, descending).
#' @param min1 Cluster-wise threshold minimum value for species shown in the final sorted synoptic table.
#' Species below that minimum will be listed in the output (\code{$others} section).
#' @param min2 Threshold minimum value for considering species values of a numerical second input table \code{syn2}.
#' Species below that minimum will not be displayed in final synoptic table, but will be listed in the
#' output (\code{$others} section).
#' @param samples Logical; if \code{TRUE}, a sorted species-sample matrix is returned. Defaults to \code{FALSE}.
#'
#' @section Details:
#' Two types of sorted synoptic tables can be created with this function:
#'
#'  \itemize{
#'  \item \code{method = "allspec"} (\emph{default}) creates a sorted synoptic table basing on one or
#'  two numeric input tables, e.g. percentage or absolute frequencies, or phi fidelity values.
#'  Sorting criteria can be either given by only one input table by using only \code{syn1}
#'  argument, as well as by two input tables with specifying \code{syn2}, too.
#'  Thereby, only values of \code{syn1} will be shown in the final sorted table.
#'  \item \code{method = "alldiff"}: With including differential species
#'  character as sorting criteria, \code{syn1} must be numeric (e.g. percentage frequency) and
#'  \code{syn2} must contain information on differential character (output from \code{\link{syntable}}
#'  function with defined \code{type = "diffspec"}). The result table shows ALL diagnostic and
#'  non-diagnostic species, as long as they match the \code{min1} and \code{min2} thresholds.
#'  The algorithm detects highest cluster values of species calculated from
#'  \code{syn1} as base for sorting, but will consider differential character criterion
#'  from \code{syn2} as well. Species with high values in \code{syn1} AND
#'  positive differential character will then be listed on the top of a species block.
#'  Within such a block, the differentiating and high-abundant species are sorted in a way favoring
#'  species that are positive in only one or at least few clusters.
#'  \item \code{method = "totalfreq"}: Sorts species simply by their overall
#'    frequency in \code{syn1} (row sums) in descending order. Only species with
#'    at least one cluster value \eqn{\ge} \code{min1} are kept in the table; the
#'    rest are listed in \code{$others}.
#'    \item \code{method = "manual"}: The rows are ordered by \code{manual_order}
#'      (species not present are ignored). Any remaining species are appended,
#'      sorted by their overall frequency (row sums of \code{syn1}, descending).
#'      Species must still pass \code{min1} (max across clusters \eqn{\ge} \code{min1});
#'      the rest go to \code{$others}.
#'  }
#'
#' @return
#' Returns an (invisible) list composed of:
#' \itemize{
#'   \item \code{$output} Sorting method description
#'   \item \code{$species} Information to species included in the output table
#'   \item \code{$samplesize} Sample sizes in clusters
#'   \item \code{$syntable} Sorted synoptic table, with the numeric values of \code{syn1} in the left-side columns
#'   and differential character of species on the right-side of the output table. See Tsiripidis et al. (2009) for
#'   details and criteria for the assignment of a differential species as p = positive, n = negative,
#'   pn = positive/negative.
#'   \item \code{$others} Species that are omitted in Synoptic table due to their failing
#'   reaching the given threshold values \code{min1} and \code{min2}. Sorted alphabetically.
#'   \item \code{$samples} Sorted original species-sample matrix, with original Plot-IDs (as column
#'   names) and the cluster identity (Cluster_No as first row of output samples table) (only when `samples = TRUE`)
#'   }
#'
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
#' Tsiripidis, I., Bergmeier, E., Fotiadis, G. & Dimopoulos, P. (2009): A new algorithm for the
#' determination of differential taxa. \emph{Journal of Vegetation Science} \strong{20}: 233-240. \doi{https://doi.org/10.1111/j.1654-1103.2009.05273.x}
#'
#' @author Jenny Schellenberg (\email{jschell@gwdg.de}), Friedemann von Lampe
#' @seealso \code{\link{syntable}}

#' @examples
#' ### Synoptic table of Scheden vegetation data using syntable()-function:
#' # classification to create a vector of cluster identity
#' library(cluster)
#' pam1 <- pam(schedenveg, 4)
#'
#'
#' ### One input table for sorting:
#' ## Synoptic table with percentage frequency of species in clusters, all species
#' unordered <- syntable(schedenveg, pam1$clustering, abund = "percentage",
#'                       type = "percfreq")   # Unordered synoptic percentage frequency table
#' sorted <- synsort(syn1 = unordered$syntable, matrix = schedenveg,
#'                   cluster = pam1$clustering, method = "allspec", min1 = 0)
#' sorted             # view results
#' \dontrun{
#' # Export sorted synoptic table
#' write.csv(sorted$syntab, "syntab.csv")
#' # Export sorted species-sample matrix with original releve data for postprocessing
#' write.csv(sorted$samples, "output_species_sample.csv")}
#'
#' ## Synoptic table with only phi values
#' phi <- syntable(schedenveg, pam1$clustering, abund = "percentage",
#'                      type = "phi")         # calculates cluster-wise phi for each species
#' phi_table <- synsort(syn1 = phi$syntable, matrix = schedenveg, cluster = pam1$clustering,
#'                      method = "allspec", min1 = 0.3)
#' phi_table     # view results
#'
#' ## Synoptic table with total frequency (global ranking)
#' total <- synsort(syn1 = unordered$syntable,
#'                cluster = pam1$clustering,
#'                method = "totalfreq",
#'                min1 = 5)
#' total         # view results
#'
#' ### Two numerical tables for sorting:
#' ## Synoptic table showing percentage frequencies, but only for species with minimum phi-value
#' ## of 0.3 AND exclude species with less than 25% percentage frequency
#'
#' unordered <- syntable(schedenveg, pam1$clustering, abund = "percentage",
#'                       type = "percfreq")   # Unordered synoptic percentage frequency table
#' phitable <- syntable(schedenveg, pam1$clustering, abund = "percentage",
#'                      type = "phi")         # calculates cluster-wise phi for each species
#' # now sorting and arranging
#' phi_complete <- synsort(syn1 = unordered$syntable, syn2 = phitable$syntable,
#'                        matrix = schedenveg, cluster = pam1$clustering, method = "allspec",
#'                        min1 = 25, min2 = 0.3)
#' phi_complete      # view results
#'
#' ### Differential species analysis
#' differential <- syntable(schedenveg, pam1$clustering, abund = "percentage",
#'                          type = "diffspec")
#'
#' ## Synoptic table with percentage frequency (only species >25%) and
#' ## differential character.
#' complete <- synsort(syn1 = unordered$syntable, syn2 = differential$syntable,
#'                     matrix = schedenveg, cluster = pam1$clustering,
#'                     method = "alldiff", min1 = 25)
#' complete            # view result table
#' differential$differentials  # list differential species for clusters
#'
#'
#' @export
#' @importFrom Hmisc %nin%



synsort <- function(syn1, syn2 = syn1, matrix = NULL, cluster, cluster_order = NULL,
                    cluster_nosort = NULL, method = "allspec", min1 = 0,
                    min2 = 0, samples = FALSE, manual_order = NULL) {
  
  cluster <- as.factor(cluster)
  if (is.null(cluster_order)) {
    cluster_order <- levels(cluster)
  } else {
    if (!all(levels(cluster) %in% cluster_order))
      stop("cluster_order must contain all cluster levels")
  }
  if (!is.null(cluster_nosort)) {
    if (!all(cluster_nosort %in% cluster_order))
      stop("cluster_nosort must contain existing cluster levels")
    # cluster_order <- c(setdiff(cluster_order, cluster_nosort), cluster_nosort)
    # NOTE: no reordering here — nosort clusters remain where they are
  }
  
  cluster <- factor(cluster, levels = cluster_order)
  
  if (!all(cluster_order %in% colnames(syn1)))
    stop("cluster_order must match column names of syn1")
  
  syn1 <- syn1[, cluster_order, drop = FALSE]
  syn2 <- syn2[, cluster_order, drop = FALSE]
  
  names(cluster) <- row.names(matrix)   # Name cluster according to site names (if present)
  group <- cluster_order
  
  if (method == "allspec") {
    if (all(syn2 == syn1)) {
      frames <- list()
      
      all <- syn1[apply(syn1, 1, max) >= min1,]
      
      for (i in 1:length(group)) {
        frames[[i]] <- assign(paste0("frame",i), all[apply(all,1,max) == all[,i],])
        frames[[i]] <- frames[[i]][sort.list(frames[[i]][,i], decreasing=TRUE),] }
      
      for ( i in 2:length(group)) {
        duprows <- rownames(frames[[i]]) %in% rownames(frames[[1]])
        frames[[1]] <- rbind(frames[[1]], frames[[i]][!duprows,]) }
      allspec <- frames[[1]]
      if (samples) {
        if (is.null(matrix))
          stop("matrix must be provided when samples = TRUE")
        specsam <- data.frame(matrix(NA, nrow = nrow(allspec), ncol = nrow(matrix)))
        
        if (length(cluster) != nrow(matrix))
          stop("cluster must have the same length as samples in matrix")
        
        rownames(specsam) <- rownames(allspec)
        matrixplotnames <- names(sort(cluster))
        names(cluster) <- seq(1, length(cluster), 1)
        colnames(specsam) <- names(sort(cluster))
        
        for (k in 1:nrow(allspec))
          for (i in 1:nrow(specsam)) {{
            if(rownames(specsam)[i] == rownames(allspec)[k]) {specsam[i,] <-
              matrix[, names(matrix) == rownames(specsam)[i]][as.numeric(names(sort(cluster)))]}
            else {}
          }}
        
        names(specsam) <- matrixplotnames
        Cluster_No <- sort(cluster); Cluster_No <- paste0(Cluster_No, "L"); specsam <- rbind(Cluster_No = Cluster_No, specsam)
        specsam[1,] <- substr(specsam[1,],1,1)
      } else {
        specsam <- data.frame()
      }
      results <- list("output" = "Synoptic table sorted by numerical values of one input table",
                      "species" = paste0("species with minimum value =", min1, " in input table 1, others listet below"),
                      "samplesize" = tapply(rep(1,length(cluster)),cluster,sum),
                      "syntable" = allspec,
                      "others" = if (length(sort(rownames(syn1[apply(syn1,1,max) < min1,]))) == 0)
                      {"No species excluded from Synoptic table."
                      } else {sort(rownames(syn1[apply(syn1,1,max) < min1,]))},
                      "samples" = specsam)
    } else {
      if (is.numeric(unlist(syn2)) == TRUE)
      { all <-  syn1[rowSums(syn2) >= min2,]
      all <-  syn1[apply(syn1,1,max) >= min1,]
      } else {stop("check data format for syn2: must be numeric")}
      frames <- list()
      for ( i in 1:length(group)) {
        frames[[i]] <- assign(paste0("frame",i), all[apply(all,1,max) == all[,i],])
        frames[[i]] <- frames[[i]][sort.list(frames[[i]][,i], decreasing=TRUE),] }
      for ( i in 2:length(group)) {
        duprows <- rownames(frames[[i]]) %in% rownames(frames[[1]])
        frames[[1]] <- rbind(frames[[1]], frames[[i]][!duprows,]) }
      
      allspec <- frames[[1]]
      if (samples) {
        specsam <- data.frame(matrix(NA, nrow = nrow(allspec), ncol = nrow(matrix)))
        rownames(specsam) <- rownames(allspec)
        
        matrixplotnames <- names(sort(cluster))
        
        names(cluster) <- seq(1, length(cluster), 1)
        colnames(specsam) <- names(sort(cluster))
        
        for (k in 1:nrow(allspec))
          for (i in 1:nrow(specsam)) {{
            if(rownames(specsam)[i] == rownames(allspec)[k]) {specsam[i,] <-
              matrix[, names(matrix) ==
                       rownames(specsam)[i]][as.numeric(names(sort(cluster)))]
            } else {}
          }}
        names(specsam) <- matrixplotnames
        
        Cluster_No <- sort(cluster); Cluster_No <- paste0(Cluster_No, "L"); specsam <- rbind(Cluster_No = Cluster_No, specsam)
        specsam[1,] <- substr(specsam[1,],1,1)
      } else {
        specsam <- data.frame()
      }
      
      results <- list("output" = "synoptic table sorted by values of two numerical input tables",
                      "species" = paste0("species with minimum value = ", min1,
                                         " in input table 1 AND with minimum value =", min2,
                                         " in input table 2, others listet below"),
                      "samplesize" = tapply(rep(1,length(cluster)),cluster,sum),
                      "syntable" = allspec,
                      "others" = if (length(sort(rownames(syn1[apply(syn1,1,max) < min1,]))) == 0)
                      {"No species excluded from Synoptic table."
                      } else {sort(rownames(syn1[apply(syn1,1,max) < min1,]))},
                      "samples" = specsam)
    }
    
    
  }  else if (method == "totalfreq") {
    # keep species that reach min1 in any cluster (same semantics as "allspec")
    allspec <- syn1[apply(syn1, 1, max) >= min1, , drop = FALSE]
    
    # order by overall frequency (row-wise sum across clusters), descending
    ord <- order(rowSums(allspec, na.rm = TRUE), decreasing = TRUE)
    allspec <- allspec[ord, , drop = FALSE]
    
    # optional samples matrix
    if (samples) {
      if (is.null(matrix))
        stop("matrix must be provided when samples = TRUE")
      if (length(cluster) != nrow(matrix))
        stop("cluster must have the same length as samples in matrix")
      
      specsam <- data.frame(matrix(NA, nrow = nrow(allspec), ncol = nrow(matrix)))
      rownames(specsam) <- rownames(allspec)
      
      matrixplotnames <- names(sort(cluster))
      names(cluster) <- seq(1, length(cluster), 1)
      colnames(specsam) <- names(sort(cluster))
      
      for (k in 1:nrow(allspec))
        for (i in 1:nrow(specsam)) {{
          if (rownames(specsam)[i] == rownames(allspec)[k]) {
            specsam[i, ] <- matrix[, names(matrix) == rownames(specsam)[i]][as.numeric(names(sort(cluster)))]
          } else {}
        }}
      
      names(specsam) <- matrixplotnames
      Cluster_No <- sort(cluster); Cluster_No <- paste0(Cluster_No, "L")
      specsam <- rbind(Cluster_No = Cluster_No, specsam)
      specsam[1, ] <- substr(specsam[1, ], 1, 1)
    } else {
      specsam <- data.frame()
    }
    
    results <- list(
      "output" = "Synoptic table sorted by total (row-wise) frequency of input table 1",
      "species" = paste0("species with minimum value = ", min1,
                         " in at least one cluster; ordered by row sums of syn1 (descending)"),
      "samplesize" = tapply(rep(1, length(cluster)), cluster, sum),
      "syntable" = allspec,
      "others" = if (length(sort(rownames(syn1[apply(syn1, 1, max) < min1, , drop = FALSE]))) == 0) {
        "No species excluded from Synoptic table."
      } else {
        sort(rownames(syn1[apply(syn1, 1, max) < min1, , drop = FALSE]))
      },
      "samples" = specsam
    )
    
  } else if (method == "manual") {
    if (is.null(manual_order))
      stop("manual_order must be provided when method = 'manual'")
    
    # keep species that reach min1 in any cluster
    allspec <- syn1[apply(syn1, 1, max) >= min1, , drop = FALSE]
    
    # sanitize / evaluate the manual list
    wanted <- unique(as.character(manual_order))
    present_manual <- wanted[wanted %in% rownames(allspec)]
    missing_manual <- setdiff(wanted, rownames(allspec))
    if (length(missing_manual))
      warning("manual_order species not found (or filtered by min1) and will be ignored: ",
              paste(missing_manual, collapse = ", "))
    
    # remainder: species not specified manually
    remainder <- setdiff(rownames(allspec), present_manual)
    if (length(remainder)) {
      rem_ord <- order(rowSums(allspec[remainder, , drop = FALSE], na.rm = TRUE),
                       decreasing = TRUE)
      remainder <- remainder[rem_ord]
    }
    
    # final row order: manual first (in given order), then remainder by total freq
    final_rows <- c(present_manual, remainder)
    allspec <- allspec[final_rows, , drop = FALSE]
    
    # optional samples matrix
    if (samples) {
      if (is.null(matrix))
        stop("matrix must be provided when samples = TRUE")
      if (length(cluster) != nrow(matrix))
        stop("cluster must have the same length as samples in matrix")
      
      specsam <- data.frame(matrix(NA, nrow = nrow(allspec), ncol = nrow(matrix)))
      rownames(specsam) <- rownames(allspec)
      
      matrixplotnames <- names(sort(cluster))
      names(cluster) <- seq(1, length(cluster), 1)
      colnames(specsam) <- names(sort(cluster))
      
      for (k in 1:nrow(allspec))
        for (i in 1:nrow(specsam)) {{
          if (rownames(specsam)[i] == rownames(allspec)[k]) {
            specsam[i, ] <- matrix[, names(matrix) == rownames(specsam)[i]][as.numeric(names(sort(cluster)))]
          } else {}
        }}
      
      names(specsam) <- matrixplotnames
      Cluster_No <- sort(cluster); Cluster_No <- paste0(Cluster_No, "L")
      specsam <- rbind(Cluster_No = Cluster_No, specsam)
      specsam[1, ] <- substr(specsam[1, ], 1, 1)
    } else {
      specsam <- data.frame()
    }
    
    results <- list(
      "output" = "Synoptic table sorted by manual species order (then total frequency)",
      "species" = paste0("species with minimum value = ", min1,
                         "; manual order applied first; remaining species by row sums (descending)"),
      "samplesize" = tapply(rep(1, length(cluster)), cluster, sum),
      "syntable" = allspec,
      "others" = if (length(sort(rownames(syn1[apply(syn1, 1, max) < min1, , drop = FALSE]))) == 0) {
        "No species excluded from Synoptic table."
      } else {
        sort(rownames(syn1[apply(syn1, 1, max) < min1, , drop = FALSE]))
      },
      "samples" = specsam
    )
    
  }  else if (method == "alldiff") {
    # setup complete table
    syntab <- syn2[apply(syn1,1,max) >= min1,]
    syntab <- syntab[complete.cases(syntab),]
    all <-    syn1[apply(syn1,1,max) >= min1,]
    all <- all[complete.cases(all),]
    completetable <- merge(all, syntab, all.x=TRUE, by= "row.names", sort=F)
    rownames(completetable) <- completetable[,1]
    completetable <- completetable[,-1]
    name <- c("")
    for(i in 1:length(group)) {
      name[i] = paste0("perc ", sort(unique(cluster))[i])
      name[i+length(unique(cluster))] = paste0("diff ", sort(unique(cluster))[i]) }
    names(completetable) <- name
    
    frames <- list()
    for (i in 1:length(group)) {
      frames[[i]] <- assign(paste0("frame",i), completetable[apply(
        completetable[,1:length(group)],1,max) == completetable[,i],])
      frames[[i]] <- frames[[i]][sort.list(frames[[i]][,i], decreasing=TRUE),] }
    for ( i in 2:length(group)) {
      duprows <- rownames(frames[[i]]) %in% rownames(frames[[1]])
      frames[[1]] <- rbind(frames[[1]], frames[[i]][!duprows,]) }
    allspec <- frames[[1]]
    
    if (samples) {
      specsam <- data.frame(matrix(NA, nrow = nrow(allspec), ncol = nrow(matrix)))
      rownames(specsam) <- rownames(allspec)
      
      matrixplotnames <- names(sort(cluster))
      
      names(cluster) <- seq(1, length(cluster), 1)
      colnames(specsam) <- names(sort(cluster))
      
      for (k in 1:nrow(allspec))
        for (i in 1:nrow(specsam)) {{
          if(rownames(specsam)[i] == rownames(allspec)[k]) {specsam[i,] <-
            matrix[, names(matrix) == rownames(specsam)[i]][as.numeric(names(sort(cluster)))]}
          else {}
        }}
      
      names(specsam) <- matrixplotnames
      
      Cluster_No <- sort(cluster); Cluster_No <- paste0(Cluster_No, "L"); specsam <- rbind(Cluster_No = Cluster_No, specsam)
      specsam[1,] <- substr(specsam[1,],1,1)
    } else {
      specsam <- data.frame()
    }
    results <- list("output" = "complete synoptic table, sorted by values of numeric input table and differential species character",
                    "species" = paste0("species with minimum value of", sep=" ", min1,
                                       " and their differentiating character"),
                    "samplesize" = tapply(rep(1,length(cluster)),cluster,sum),
                    "syntable" = allspec,
                    "others" = if (length(sort(rownames(syn1[apply(syn1,1,max) < min1,]))) == 0)
                    {"No species excluded from Synoptic table."
                    } else {sort(rownames(syn1[apply(syn1,1,max) < min1,]))},
                    "samples" = specsam)
  } else {stop("Sorting of synoptic table failed: wrong method entry. Check correct formula input")}
  
  return(invisible(results))
}