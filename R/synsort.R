#' Sorting functions for synoptic tables
#'
#' @description Synoptic tables are a tool for the visualization and interpretation of previously
#' defined plant species groups (clusters), e.g. from cluster analysis or classification methods.
#' They help to determine characteristic patterning of species occurrences in plant communities
#' by calculating cluster-wise percentage or absolute frequencies, mean/median cover values, fidelity
#' (phi) or differential species character.
#'
#' This function sorts synoptic tables from \code{\link{syntable}} function output. Sorting criteria
#' can be either numerical values in synoptic tables, such as cluster-wise frequencies or fidelity
#' measures, as well as combined criteria with considering differential character, too (according to
#' the criteria defined by Tsiripidis et al., 2009).
#'
#' The algorithm aims to sort species to blocked structure considering the defined criteria and input
#' tables, with the best characterizing species on the top of the block, followed by species with
#' descending importance for plant community description.
#'
#' @param syn1 Input synoptic table 1, a dataframe with numerical data format, usually from
#' \code{\link{syntable}} function output. See Details for input table format.
#' The values of this table will be displayed in the final output table.
#' @param syn2 Optional second input table with additional numeric or differential character
#' sorting criteria.
#' @param original Species-sample matrix, already used for \code{syntable()} function input (there:
#' \code{spec})
#' @param cluster Vector with classification cluster identity, named with the unique plot IDs,
#' both in integer format. Ensure matching order of cluster identity in the cluster vector
#' and samples in used dataframe ( = original) for correct allocation of cluster numbers to samples.
#' @param method Sorting algorithm and synoptic table output options (\code{method = c("allspec", "alldiff")}).
#' See Details.
#' @param min1 Cluster-wise threshold minimum value for species shown in the final sorted synoptic table.
#' Species below that minimum will be listed in the  output (\code{$others} section).
#' @param min2 Threshold minimum value for considering species values of a numerical second input table \code{syn2}.
#' Species below that minimum will not be displayed in final synoptic table, but will be listed in the
#' output (\code{$others} section).
#'
#' @section Details:
#' Two types of sorted synoptic tables can be created with this function:
#'
#'  \itemize{
#'  \item{\code{method = "allspec"} (default)}{ creates a sorted synoptic table basing on one or
#'  two numeric input tables, e.g. percentage or absolute frequencies, or phi fidelity values.
#'  Sorting criteria can be either given by only one input table by using only \code{syn1}
#'  argument, as well as by two input tables with specifying \code{syn2}, too.
#'  Thereby, only values of \code{syn1} will be shown in the final sorted table.}
#'  \item{\code{method = "alldiff"}:}{ With including differential species
#'  character as sorting criteria, \code{syn1} must be numeric (e.g. percentage frequency) and
#'  \code{syn2} must contain information on differential character (output from \code{\link{syntable}}
#'  function with defined \code{type = "diffspec"}). The result table shows ALL diagnostic and
#'  non-diagnostic species, as long as they match the \code{min1} and \code{min2} thresholds.
#'  The algorithm detects highest cluster values of species calculated from
#'  \code{syn1} as base for sorting, but will consider differential character criterion
#'  from \code{syn2} as well. Species with high values in \code{syn1} AND
#'  positive differential character will then be listed on the top of a species block.
#'  Within such a block, the differentiating and high-abundant species are sorted in a way favoring
#'  species that are positive in only one or at least few clusters.} }
#'
#'
#'
#' @return
#' Returns a list composed of:
#'   \item{\code{$output}}{Sorting method description}
#'   \item{\code{$species}}{Information to species included in the output table}
#'   \item{\code{$samplesize}}{Sample sizes in clusters}
#'   \item{\code{$syntable}}{Sorted synoptic table, with the numeric values of code{syn1} in the left-side columns
#'   and differential character of species on the right-side of the output table. See Tsiripidis et al. (2009) for
#'   details and criteria for the assignment of a differential species as p = positive, n = negative,
#'   pn = positive/negative.}
#'   \item{\code{$others}}{Species that are omitted in Synoptic table due to their failing
#'   reaching the given threshold values \code{min1} and \code{min2}. Sorted alphabetically.}
#'   \item{\code{$samples}}{Sorted original species-sample matrix, with original Plot-IDs (as column
#'   names) and the cluster identity (Cluster_No as first row of output samples table)}
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
#' @author Jenny Schellenberg (\email{jschell@gwdg.de})
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
#' unordered <- syntable(schedenveg, pam1$clustering, abund = "perc",
#'                       type = "percfreq")   # Unordered synoptic percentage frequency table
#' sorted <- synsort(syn1 = unordered$syntable, original = schedenveg,
#'                   cluster = pam1$clustering, method = "allspec", min1 = 0)
#' sorted             # view results
#' \dontrun{
#' # export sorted synoptic table
#' write.csv(sorted$syntab, "syntab.csv")
#' # export sorted species-sample matrix with original releve data for postprocessing
#' write.csv(sorted$samples, "output_species_sample.csv")}
#'
#' ## Synoptic table with only phi values
#' phi <- syntable(schedenveg, pam1$clustering, abund = "perc",
#'                      type = "phi")         # calculates cluster-wise phi for each species
#' phi_table <- synsort(syn1 = phi$syntable, original = schedenveg, cluster = pam1$clustering,
#'                      method = "allspec", min1 = 0.3)
#' phi_table     # view results
#'
#'
#' ### Two numerical tables for sorting:
#' ## Synoptic table showing percentage frequencies, but only for species with minimum phi-value
#' ## of 0.3 AND exclude species with less than 25% percentage frequency
#'
#' unordered <- syntable(schedenveg, pam1$clustering, abund = "perc",
#'                       type = "percfreq")   # Unordered synoptic percentage frequency table
#' phitable <- syntable(schedenveg, pam1$clustering, abund = "perc",
#'                      type = "phi")         # calculates cluster-wise phi for each species
#' # now sorting and arranging
#' phi_complete <- synsort(syn1 = unordered$syntable, syn2 = phitable$syntable,
#'                        original = schedenveg, cluster = pam1$clustering, method = "allspec",
#'                        min1 = 25, min2 = 0.3)
#' phi_complete      # view results
#'
#' ### Differential species analysis
#' differential <- syntable(schedenveg, pam1$clustering, abund = "perc",
#'                          type = "diffspec")
#'
#' ## Synoptic table with percentage frequency (only species >25%) and
#' ## differential character.
#' complete <- synsort(syn1 = unordered$syntable, syn2 = differential$syntable,
#'                     original = schedenveg, cluster = pam1$clustering,
#'                     method = "alldiff", min1 = 25)
#' complete            # view result table
#' differential$differentials  # list differential species for clusters
#'
#' @export
#' @importFrom Hmisc %nin%




synsort <- function(syn1, syn2 = syn1 , original, cluster, method = "allspec", min1 = 0,
                    min2 = 0) {
  if (method == "allspec") {
      if (all(syn2 == syn1)) {
      frames <- list()
      all <- syn1[apply(syn1,1,max) >= min1,]
      for (i in 1:max(cluster)) {
        frames[[i]] <- assign(paste0("frame",i), all[apply(all,1,max) == all[,i],])
        frames[[i]] <- frames[[i]][sort.list(frames[[i]][,i], decreasing=TRUE),] }
      for ( i in 2:max(cluster)) {
        duprows <- rownames(frames[[i]]) %in% rownames(frames[[1]])
        frames[[1]] <- rbind(frames[[1]], frames[[i]][!duprows,]) }
      allspec <- frames[[1]]
      specsam <- data.frame(matrix(NA, nrow = nrow(allspec), ncol = nrow(original)))
      rownames(specsam) <- rownames(allspec)
      originalplotnames <- names(sort(cluster))
      names(cluster) <- seq(1, length(cluster), 1)
      colnames(specsam) <- names(sort(cluster))
      for (k in 1:nrow(allspec))
        for (i in 1:nrow(specsam)) {{
          if(rownames(specsam)[i] == rownames(allspec)[k]) {specsam[i,] <-
            original[, names(original) == rownames(specsam)[i]][as.numeric(names(sort(cluster)))]}
          else {}
        }}
      names(specsam) <- originalplotnames
      Cluster_No <- sort(cluster); Cluster_No <- paste0(Cluster_No, "L"); specsam <- rbind(Cluster_No = Cluster_No, specsam)
      specsam[1,] <- substr(specsam[1,],1,1)
      results <- list("output" = "Synoptic table sorted by numerical values of one input table",
                      "species" = paste0("species with minimum value =", min1, " in input table 1, others listet below"),
                      "samplesize" = tapply(rep(1,length(cluster)),cluster,sum),
                      "syntable" = allspec,
                      "others" = if (length(sort(rownames(syn1[apply(syn1,1,max)<=min1,]))) == 0)
                                  {"No species excluded from Synoptic table."
                                } else {sort(rownames(syn1[apply(syn1,1,max)<=min1,]))},
                      "samples" = specsam)
      } else {
          if (is.numeric(unlist(syn2)) == TRUE)
            { all <-  syn1[rowSums(syn2) >= min2,]
              all <-  syn1[apply(syn1,1,max) >= min1,]
            } else {stop("check data format for syn2: must be numeric")}
          frames <- list()
          for ( i in 1:max(cluster)) {
            frames[[i]] <- assign(paste0("frame",i), all[apply(all,1,max) == all[,i],])
            frames[[i]] <- frames[[i]][sort.list(frames[[i]][,i], decreasing=TRUE),] }
          for ( i in 2:max(cluster)) {
            duprows <- rownames(frames[[i]]) %in% rownames(frames[[1]])
            frames[[1]] <- rbind(frames[[1]], frames[[i]][!duprows,]) }
          allspec <- frames[[1]]
          specsam <- data.frame(matrix(NA, nrow = nrow(allspec), ncol = nrow(original)))
          rownames(specsam) <- rownames(allspec)
          originalplotnames <- names(sort(cluster))
          names(cluster) <- seq(1, length(cluster), 1)
          colnames(specsam) <- names(sort(cluster))
          for (k in 1:nrow(allspec))
            for (i in 1:nrow(specsam)) {{
              if(rownames(specsam)[i] == rownames(allspec)[k]) {specsam[i,] <-
                original[, names(original) ==
                           rownames(specsam)[i]][as.numeric(names(sort(cluster)))]
              } else {}
            }}
          names(specsam) <- originalplotnames
          Cluster_No <- sort(cluster); Cluster_No <- paste0(Cluster_No, "L"); specsam <- rbind(Cluster_No = Cluster_No, specsam)
          specsam[1,] <- substr(specsam[1,],1,1)
          results <- list("output" = "synoptic table sorted by values of two numerical input tables",
                        "species" = paste0("species with minimum value = ", min1,
                                           " in input table 1 AND with minimum value =", min2,
                                           " in input table 2, others listet below"),
                        "samplesize" = tapply(rep(1,length(cluster)),cluster,sum),
                        "syntable" = allspec,
                        "others" = if (length(sort(rownames(syn1[apply(syn1,1,max)<=min1,]))) == 0)
                        {"No species excluded from Synoptic table."
                        } else {sort(rownames(syn1[apply(syn1,1,max)<=min1,]))},
                        "samples" = specsam)
      }
    } else if (method == "alldiff") {
      # setup complete table
      syntab <- syn2[apply(syn1,1,max) >= min1,]
      syntab <- syntab[complete.cases(syntab),]
      all <-    syn1[apply(syn1,1,max) >= min1,]
      all <- all[complete.cases(all),]
      completetable <- merge(all, syntab, all.x=TRUE, by= "row.names", sort=F)
      rownames(completetable) <- completetable[,1]
      completetable <- completetable[,-1]
      name <- c("")
      for(i in 1:length(unique(cluster))) {
        name[i] = paste0("perc ", sort(unique(cluster))[i])
        name[i+length(unique(cluster))] = paste0("diff ", sort(unique(cluster))[i]) }
      names(completetable) <- name
      frames <- list()
      for (i in 1:max(cluster)) {
        frames[[i]] <- assign(paste0("frame",i), completetable[apply(
          completetable[,1:max(cluster)],1,max) == completetable[,i],])
        frames[[i]] <- frames[[i]][sort.list(frames[[i]][,i], decreasing=TRUE),] }
      for ( i in 2:max(cluster)) {
        duprows <- rownames(frames[[i]]) %in% rownames(frames[[1]])
        frames[[1]] <- rbind(frames[[1]], frames[[i]][!duprows,]) }
      allspec <- frames[[1]]
      specsam <- data.frame(matrix(NA, nrow = nrow(allspec), ncol = nrow(original)))
      rownames(specsam) <- rownames(allspec)
      originalplotnames <- names(sort(cluster))
      names(cluster) <- seq(1, length(cluster), 1)
      colnames(specsam) <- names(sort(cluster))
      for (k in 1:nrow(allspec))
        for (i in 1:nrow(specsam)) {{
          if(rownames(specsam)[i] == rownames(allspec)[k]) {specsam[i,] <-
            original[, names(original) == rownames(specsam)[i]][as.numeric(names(sort(cluster)))]}
          else {}
          }}
      names(specsam) <- originalplotnames
      Cluster_No <- sort(cluster); Cluster_No <- paste0(Cluster_No, "L"); specsam <- rbind(Cluster_No = Cluster_No, specsam)
      specsam[1,] <- substr(specsam[1,],1,1)
      results <- list("output" = "complete synoptic table, sorted by values of numeric input table and differential species character",
                    "species" = paste0("species with minimum value of", sep=" ", min1,
                                       " and their differentiating character"),
                    "samplesize" = tapply(rep(1,length(cluster)),cluster,sum),
                    "syntable" = allspec,
                    "others" = if (length(sort(rownames(syn1[apply(syn1,1,max)<=min1,]))) == 0)
                    {"No species excluded from Synoptic table."
                    } else {sort(rownames(syn1[apply(syn1,1,max)<=min1,]))},
                    "samples" = specsam)
    } else {stop("Sorting of synoptic table failed: wrong method entry. Check correct formula input")}

}



