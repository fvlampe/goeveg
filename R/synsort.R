#' Sorting functions for synoptic tables
#'
#' @description Synoptic tables are a tool for interpretation of cluster species composition.
#' This function provides sorting options for synoptic tables, sorting criteria can be either
#' values in synoptic tables, such as frequencies, as well as combined criteria with considering
#' differential character, too.
#' Sorting algorithm aims to sort species in given cluster column order to blocked structure.
#' Thereby, species with high frequencies and/or differential character are displayed blocked
#' for each cluster or several neighbouring clusters.
#'
#' @param syn1 Input synoptic table 1 (as dataframe) with priority entries for sorting.
#' Usually dataframe from \code{\link{syntable}} function output,
#' but function should work with every synoptic table input, as long as formats are
#' appropriate. The values of this table will be displayed in the final output table.
#' @param syn2 Optional second input table with additional sorting criteria. Note that
#' values of second input table will be considered in sorting, but not be displayed
#' in final synoptic table with \code{method = "allspec"}.
#' @param cluster Integer vector with classification cluster identity. Ensure matching order
#' of cluster identity and samples in dataframe for correct allocation of cluster numbers to samples.
#' @param method Sorting algorithm (\code{method = c("allspec", "p_diff", "n_diff", "pn_diff", "accspec", "all_diff")}).
#' See Details.
#' @param min1 Treshold minimum value for considering species of \code{syn1} in ordering algorithm.
#' Species below that minimum will neither be considered in algorithm nor displayed in final
#' synoptic table, but will be listed in the \code{$others} output.
#' @param min2 Treshold minimum value for considering species of \code{syn2} in ordering algorithm.
#' Species below that minimum will neither be considered in algorithm nor displayed in final
#' synoptic table, but will be listed in the \code{$others} vector.
#' @param relate2 Specifies relation of given second table minimum values to either related to
#' entire dataset (default) or to each cluster only (\code{relate2 = c("entire", "cluster"))}.
#'
#' @section Details:
#' Six types of synoptic tables can be created with this function.
#'
#'  \describe{\item{\code{method = "allspec"}}{creates a sorted synoptic table basing on numeric input tables,
#'  e.g. common percentage frequency tables. Sorting criteria can be either given by one input table (\code{syn1}), as well
#'  as by two input tables (\code{syn1, syn2}). Thereby, only values of \code{syn1} will be shown in the final sorted table,
#'  but values of \code{syn2} will be considered, too. The second minimum treshold (\code{min2}) for values in \code{syn2}
#'  will be either applied on single clusters (e.g. recommended for phi values with specifying \code{relate2 = "cluster"})
#'  for having minimum phi/cluster or e.g. on total frequencies in entire dataset for excluding rare species from
#'  synoptic table, applying the minimum treshold on the entire dataset (\code{relate2 = "entire"}).}}
#'
#'  With including differential species character as sorting criterion (\code{method = c("p_diff", "n_diff", "pn_diff", "accspec", "all_diff")}),
#'  input table \code{syn1} must be numeric, the second one with information on differential character (output from \code{\link{syntable}} function with
#' \code{type="diffspec"}). Again, algorithm detects highest cluster values of species in \code{syn1} as base for sorting,
#' but will sort them considering differentiating character criterion (from second
#' input table \code{syn2}). Species with high values in \code{syn1} AND differential character will then be listed
#' on the top of a species block. Within differentiating species, prevalence of diagnostic character
#' is considered by favoring positive and/or cluster-specific differential character. Available types are:
#' \describe{
#'  \item{\code{method = "p_diff"}}{creates a synoptic table of diagnostic species with numerical values of input table
#'  \code{syn1}}
#'  \item{\code{method = "accspec"}}{creates a synoptic table of non-diagnostic species with numerical values of input table
#'  \code{syn1}}
#'  \item{\code{method = "all_diff"}}{showing all diagnostic and non-diagnostic species}}
#'
#' @return
#' Returns a list composed of:
#'   \item{\code{$output}}{sorting method description}
#'   \item{\code{$species}}{species sorting criteria}
#'   \item{\code{$samplesize}}{sample sizes in clusters}
#'   \item{\code{$syntable}}{sorted synoptic table}
#'   \item{\code{$others}}{species that failed to be included in the final table due to
#' treshold values given by min1 and min2}
#'    \item{\code{$differential}}{In case of combined sorting with considering differential species character,
#'    a table with differential character of species.}
#'
#' @references
#' Bruelheide, H. (2000): A new measure of fidelity and its application to defining species groups. - \emph{Journal
#' of Vegetation Science} 11: 167-178.
#'
#' Chytry, M., Tichy, L., Holt, J., Botta-Dukat, Z. (2002): Determination of diagnostic species with statistical fidelity measures.
#' \emph{Journal of Vegetation Science} 13: 79-90.
#'
#' Sokal, R.R. & Rohlf, F.J. (1995): Biometry. 3rd edition Freemann, New York.
#'
#' Tsiripidis, I., Bergmeier, E., Fotiadis, G. & Dimopoulos, P. (2009): A new algorithm for the determination
#' of differential taxa. - \emph{Journal of Vegetation Science} 20: 233-240.
#'
#' @author Jenny Schellenberg (\email{jschell@gwdg.de})
#' @seealso \code{\link{syntable}}
#' @examples
#' ## Synoptic table of Scheden vegetation data:
#' library(cluster)
#' pam1 <- pam(schedenveg, 4)
#'
#' ## Unordered synoptic tables
#' # Unordered synoptiv percentage frequency table
#' unordered <- syntable(schedenveg, pam1$clustering, abund = "perc",
#'                       type = "percfreq")
#' # Differential species analysis
#' differential <- syntable(schedenveg, pam1$clustering, abund = "perc",
#'                          type = "diffspec")
#' # Fidelity phi
#' phitable <- syntable(schedenveg, pam1$clustering, abund = "perc",
#'                      type = "phi")
#'
#' ## Common complete synoptic table: sort by percentage frequency,
#' ## show all species
#' sorted <- synsort(syn1 = unordered$syntable, cluster = pam1$clustering,
#'                   method = "allspec", min1 = 0)
#' sorted             # view results
#'
#' ## Synoptic table, with only positive differentiating species with
#' ## minimum 25% frequency in table
#' positive <- synsort(syn1 = unordered$syntable, syn2 = differential$syntable,
#'                     cluster = pam1$clustering, method = "p_diff", min1 = 25)
#' positive           # view results
#'
#' ## Synoptic table, with percentage frequency (only species >25%) and
#' ## differential character.
#' complete <- synsort(syn1 = unordered$syntable, syn2 = differential$syntable,
#'                     cluster = pam1$clustering, method = "all_diff", min1 = 25)
#' complete
#'
#' ## Synoptic table, species with minimum phi-value of 0.3, show
#' ## percentage frequency
#' phi_complete <- synsort(syn1 = unordered$syntable, syn2 = phitable$syntable,
#'                         cluster = pam1$clustering, method = "allspec", min1 = 25, min2 = 0.3)
#' phi_complete
#'
#' ## Synoptic table with only phi values
#' phi_table <- synsort(syn1 = phitable$syntable, cluster = pam1$clustering,
#'                      method = "allspec", min1 = 0.3)
#' phitable
#'
#' ## Synoptic table showing diagnostic phi phi-values (>=0.3) and
#' ## differential character
#' phi_diff <- synsort(syn1 = phitable$syntable, syn2 = differential$syntable,
#'                     cluster = pam1$clustering, method = "all_diff", min1 = 0.3)
#' phi_diff
#' @export
#' @importFrom Hmisc %nin%




synsort <- function(syn1, syn2 = syn1 , cluster, method = "allspec", min1 = 0, min2 = 0, relate2 = "entire") {

if (method == "allspec") {
      if (all(syn2 == syn1)) {
      frames <- list()
      all <- syn1[apply(syn1,1,max) >= min1,]
      for (i in 1:max(cluster)) {
        frames[[i]] <- assign(paste0("frame",i), all[apply(all,1,max) == all[,i],])
        frames[[i]] <- frames[[i]][order(frames[[i]][,i], decreasing=TRUE),] }
      for ( i in 2:max(cluster)) {
        duprows <- rownames(frames[[i]]) %in% rownames(frames[[1]])
        frames[[1]] <- rbind(frames[[1]], frames[[i]][!duprows,]) }
      allspec <- frames[[1]]
      results <- list("output" = "synoptic table sorted by values of one input table",
                      "species" = paste0("species with minimum value =", min1, " in input table 1, others listet below"),
                      "samplesize" = tapply(rep(1,length(cluster)),cluster,sum),
                      "syntable" = allspec,
                      "others" = sort(rownames(syn1[apply(syn1,1,max)<=min1,])))
      return(results)

      } else {
          if (relate2 == "entire") {
            if (is.numeric(unlist(syn2)) == TRUE)
            { all <-  syn1[rowSums(syn2) >= min2,]
              all <-  syn1[apply(syn1,1,max) >= min1,]
            } else if (is.character(unlist(syn2)) == TRUE)
              {  all <-  syn1[apply(syn1,1,max) >= min1,]
              } else {stop("check data entry type for syn2: must be numeric or character")}

          frames <- list()
          for ( i in 1:max(cluster)) {
            frames[[i]] <- assign(paste0("frame",i), all[apply(all,1,max) == all[,i],])
            frames[[i]] <- frames[[i]][order(frames[[i]][,i], decreasing=TRUE),] }
          for ( i in 2:max(cluster)) {
            duprows <- rownames(frames[[i]]) %in% rownames(frames[[1]])
            frames[[1]] <- rbind(frames[[1]], frames[[i]][!duprows,]) }
          allspec <- frames[[1]]
          results <- list("output" = "synoptic table sorted by values of one input table",
                        "species" = paste0("species with minimum value = ", min1,
                                           " in input table 1 AND with minimum value =", min2, " in input table 2,others listet below"),
                        "samplesize" = tapply(rep(1,length(cluster)),cluster,sum),
                        "syntable" = allspec,
                        "others" = sort(rownames(syn1[apply(syn1,1,max)<=min1,])))
          return(results)

          } else {
            all = syn1[apply(syn2,1,min) >= min2,]
            all = syn1[apply(syn1,1,max) >= min1,]
            for ( i in 1:max(cluster)) {
              frames[[i]] <- assign(paste0("frame",i), all[apply(all,1,max) == all[,i],])
              frames[[i]] <- frames[[i]][order(frames[[i]][,i], decreasing=TRUE),] }
            for ( i in 2:max(cluster)) {
              duprows <- rownames(frames[[i]]) %in% rownames(frames[[1]])
              frames[[1]] <- rbind(frames[[1]], frames[[i]][!duprows,]) }
            allspec <- frames[[1]]
            results <- list("output" = "synoptic table sorted by values of one input table",
                            "species" = paste0("species with minimum value = ", min1,
                                               " in input table 1 AND with minimum value =", min2, " in input table 2,others listet below"),
                            "samplesize" = tapply(rep(1,length(cluster)),cluster,sum),
                            "syntable" = allspec,
                            "others" = sort(rownames(syn1[apply(syn1,1,max)<=min1,])))
            return(results)
          }
      }

} else if (method =="p_diff" | method == "accspec" | method == "all_diff") {
    # setup complete table
    syntab <- syn2[apply(syn1,1,max) >= min1,]
    syntab <- syntab[complete.cases(syntab),]
    all <-    syn1[apply(syn1,1,max) >= min1,]
    all <- all[complete.cases(all),]
    completetable <- merge(syntab, all, all.x=TRUE, by= "row.names", sort=F)
    rownames(completetable) <- completetable[,1]
    completetable <- completetable[,-1]
    name <- c("")
    for(i in 1:length(unique(cluster))) {
      name[i] = paste0("diff ", sort(unique(cluster))[i])
      name[i+length(unique(cluster))] = paste0("cluster ", sort(unique(cluster))[i]) }
    names(completetable) <- name


  # posidiff: a vector with sums of clusters "p"-differentiating
  # by the species, respectively
    pos_synord <- data.frame(matrix(NA, ncol = max(cluster), nrow = length(syntab[,1])))
    for ( i in 1:length(syntab[1,])) {
        for ( k in 1:length(syntab[,1])) {
          if (syntab[k,i] == "p") {pos_synord[k,i] <- 1}
          else {pos_synord[k,i] <- 0}   }   }
    rownames(pos_synord) <- rownames(syntab)
    pos_synord <- data.frame(lapply(pos_synord, as.numeric))
    rownames(pos_synord) <- rownames(all)
    names(pos_synord) <- seq(1,length(unique(cluster)),1)
    posidiff <- apply(pos_synord,1,sum)

  # accdiff: a vector with sums of clusters "-" non-differentiating
    acc_synord <- data.frame(matrix(NA, ncol = max(cluster), nrow = length(syntab[,1])))
    for ( i in 1:length(syntab[1,])) {
      for ( k in 1:length(syntab[,1])) {
        if (syntab[k,i] == "-") {acc_synord[k,i] <- 1}
        else {acc_synord[k,i] <- 0}   }   }

    acc_synord <- data.frame(lapply(acc_synord, as.numeric))
    rownames(acc_synord) <- rownames(syntab)
    names(acc_synord) <- seq(1,length(unique(cluster)),1)
    accdiff <- apply(acc_synord,1,sum) # vector with sums of
    # clusters "-"-not differentiated by the species, respectively

    # selections from all non-differentianting species from the whole dataset
    diffacc <- syntab[accdiff == max(accdiff),]
    synacc <- all[accdiff == max(accdiff),]

    # for method = "accspec": show percentage frequency of
    # all non-differentiating species with values (syn1) > min1
    acctable <- merge(synacc, diffacc, by= "row.names", all.x = TRUE, sort=F)
    rownames(acctable) <- acctable[,1]
    acctable <- acctable[,-1]
    frames <- list()
    for (i in 1:length(unique(cluster))) {
      frames[[i]] <- assign(paste0("frame",i),
                     acctable[apply(acctable[, min(sort(unique(cluster))):
                                               max(sort(unique(cluster)))],
                                               1, max) == acctable[,(i)],])
      frames[[i]] <- frames[[i]][order(frames[[i]][,i], decreasing=TRUE),]
      frames[[i]] <- frames[[i]] [complete.cases(frames[[i]]),] }
    acctable <- do.call("rbind", frames)
    acctable <- acctable[substr(rownames(acctable), nchar(rownames(acctable)), nchar(rownames(acctable))) != 1 &
                    substr(rownames(acctable), nchar(rownames(acctable)), nchar(rownames(acctable))) != 2 &
                    substr(rownames(acctable), nchar(rownames(acctable)), nchar(rownames(acctable))) != 3 &
                    substr(rownames(acctable), nchar(rownames(acctable)), nchar(rownames(acctable))) != 4 &
                    substr(rownames(acctable), nchar(rownames(acctable)), nchar(rownames(acctable))) != 5 &
                    substr(rownames(acctable), nchar(rownames(acctable)), nchar(rownames(acctable))) != 6 &
                    substr(rownames(acctable), nchar(rownames(acctable)), nchar(rownames(acctable))) != 7 &
                    substr(rownames(acctable), nchar(rownames(acctable)), nchar(rownames(acctable))) != 8 &
                    substr(rownames(acctable), nchar(rownames(acctable)), nchar(rownames(acctable))) != 9,]
    names(acctable) <- rep(c(seq(1,max(cluster), 1)),2)
    for ( i in 1:max(cluster)) {
          names(acctable)[max(cluster) + i] <- paste0("diff", names(acctable)[i])
        }

    if (method == "accspec") {
        results <- list(
          "output" = "synoptic table sorted by values of numeric input table and non-differential, accompanying species",
          "species" = paste0("species with minimim value of", sep=" ", min1, " in input numeric table and that are ACCOMPANYING"),
          "samplesize" = tapply(rep(1,length(cluster)),cluster,sum),
          "syntable" = acctable,
          "diagnostic species" = sort(rownames(syn1)[Hmisc::`%nin%` (rownames(syn1), rownames(synacc))]),
          "others" = sort(rownames(syn1)[Hmisc::`%nin%` (rownames(syn1), rownames(acctable))]))
        return(results)

    } else {                                          # for method  == "alldiff"
        # sort for each column for blocked structure
        frames1 <- list()
        in1 <- completetable[posidiff == 1,]   # for species only differentiating for one
        for (i in 1:length(unique(cluster))) {
          frames1[[i]] = assign(paste0("frame1_",i), in1[in1[,i] == "p",])
          frames1[[i]] = frames1[[i]][order(frames1[[i]][,i + length(unique(cluster))], decreasing=TRUE),]
        if (length(rownames(frames1[[i]])) == 0) {
          } else { rownames(frames1[[i]]) <- paste0(rownames(frames1[[i]]), "_")  }  }
        pos_in1 <- do.call("rbind", frames1)


        in2 <- completetable[posidiff == 2,]
        maxi <- apply(in2[, (1 + length(unique(cluster))) : (2*length(unique(cluster)))], 1, max)
        frames2 <- list()
        for (i in 1:length(unique(cluster))) {
           frames2[[i]] = assign(paste0("frame2_",i), in2[in2[, i + length(unique(cluster))] == maxi, ])
           frames2[[i]] = frames2[[i]][order(frames2[[i]][, i + length(unique(cluster))], decreasing =T),]
           if (length(rownames(frames2[[i]])) == 0) {
             } else { rownames(frames2[[i]]) <- paste0(rownames(frames2[[i]]), "_")  }   }
        pos_in2 <- unique(do.call("rbind", frames2))

        in3 <- completetable[posidiff == 3,]
        maxi <- apply(in3[, (1 + length(unique(cluster))) : (2*length(unique(cluster)))], 1, max)
        frames3 <- list()
        for (i in 1:length(unique(cluster))) {
          frames3[[i]] = assign(paste0("frame3_",i), in3[in3[, i + length(unique(cluster))] == maxi, ])
          frames3[[i]] = frames3[[i]][order(frames3[[i]][, i + length(unique(cluster))], decreasing=TRUE),]
          if (length(rownames(frames3[[i]])) == 0) {
            } else { rownames(frames3[[i]]) <- paste0(rownames(frames3[[i]]), "_")  }   }
        pos_in3 <- unique(do.call("rbind", frames3))

        in4 <- completetable[posidiff == 4,]
        maxi <- apply(in4[, (1 + length(unique(cluster))) : (2*length(unique(cluster)))], 1, max)
        frames4 <- list()
        for (i in 1:length(unique(cluster))) {
          frames4[[i]] = assign(paste0("frame4_",i), in4[in4[, i + length(unique(cluster))] == maxi, ])
          frames4[[i]] = frames4[[i]][order(frames4[[i]][,i + length(unique(cluster))], decreasing=TRUE),]
        if (length(rownames(frames4[[i]])) == 0) {
        } else { rownames(frames4[[i]]) <- paste0(rownames(frames4[[i]]), "_")  }   }
        pos_in4 <- unique(do.call("rbind", frames4))

        in5 <- completetable[posidiff == 5,]
        maxi <- apply(in5[, (1 + length(unique(cluster))) : (2*length(unique(cluster)))], 1, max)
        frames5 <- list()
        for (i in 1:length(unique(cluster))) {
          frames5[[i]] = assign(paste0("frame5_",i), in5[in5[, i + length(unique(cluster))] == maxi, ])
          frames5[[i]] = frames5[[i]][order(frames5[[i]][,i + length(unique(cluster))], decreasing=TRUE),]
        if (length(rownames(frames5[[i]])) == 0) {
          } else { rownames(frames5[[i]]) <- paste0(rownames(frames5[[i]]), "_")  }   }
        pos_in5 <- unique(do.call("rbind", frames5))

        in6 <- completetable[posidiff == 6,]
        maxi <- apply(in6[, (1 + length(unique(cluster))) : (2*length(unique(cluster)))], 1, max)
        frames6 <- list()
        for (i in 1:length(unique(cluster))) {
          frames6[[i]] = assign(paste0("frame6_",i), in6[in6[, i + length(unique(cluster))] == maxi, ])
          frames6[[i]] = frames6[[i]][order(frames6[[i]][,i + length(unique(cluster))], decreasing=TRUE),]
        if (length(rownames(frames6[[i]])) == 0) {
          } else { rownames(frames6[[i]]) <- paste0(rownames(frames6[[i]]), "_")  }   }
        pos_in6 <- unique(do.call("rbind", frames6))

        in7 <- completetable[posidiff == 7,]
        maxi <- apply(in7[, (1 + length(unique(cluster))) : (2*length(unique(cluster)))], 1, max)
        frames7 <- list()
        for (i in 1:length(unique(cluster))) {
          frames7[[i]] = assign(paste0("frame7_",i), in7[in7[, i + length(unique(cluster))] == maxi, ])
          frames7[[i]] = frames7[[i]][order(frames7[[i]][,i + length(unique(cluster))], decreasing=TRUE),]
          if (length(rownames(frames7)) == 0) {
          } else { rownames(frames7[[i]]) <- paste0(rownames(frames7[[i]]), "_")  }   }
        pos_in7 <- unique(do.call("rbind", frames7))

        in8 <- completetable[posidiff >= 8,]
        frames8 <- list()
        maxi <- apply(in8[, (1 + length(unique(cluster))) : (2*length(unique(cluster)))], 1, max)
        for (i in 1:length(unique(cluster))) {
          frames8[[i]] = assign(paste0("frame8_",i), in8[in8[, i + length(unique(cluster))] == maxi, ])
          frames8[[i]] = frames8[[i]][order(frames8[[i]][,i + length(unique(cluster))], decreasing=TRUE),]
        if (length(rownames(frames8)) == 0) {
          } else { rownames(frames8[[i]]) <- paste0(rownames(frames8[[i]]), "_")  }   }
        pos_in8 <- unique(do.call("rbind", frames8))

    completetable <- rbind(pos_in1, pos_in2, pos_in3, pos_in4, pos_in5, pos_in6, pos_in7, pos_in8)
    completetable <- completetable[substr(rownames(completetable), nchar(rownames(completetable)), nchar(rownames(completetable))) != 1 &
                                   substr(rownames(completetable), nchar(rownames(completetable)), nchar(rownames(completetable))) != 2 &
                                   substr(rownames(completetable), nchar(rownames(completetable)), nchar(rownames(completetable))) != 3 &
                                   substr(rownames(completetable), nchar(rownames(completetable)), nchar(rownames(completetable))) != 4 &
                                   substr(rownames(completetable), nchar(rownames(completetable)), nchar(rownames(completetable))) != 5 &
                                   substr(rownames(completetable), nchar(rownames(completetable)), nchar(rownames(completetable))) != 6 &
                                   substr(rownames(completetable), nchar(rownames(completetable)), nchar(rownames(completetable))) != 7 &
                                   substr(rownames(completetable), nchar(rownames(completetable)), nchar(rownames(completetable))) != 8 &
                                   substr(rownames(completetable), nchar(rownames(completetable)), nchar(rownames(completetable))) != 9,]
    rownames(completetable) <- substr(rownames(completetable),
                                      nchar(rownames(completetable))-nchar(rownames(completetable))+1,
                                      nchar(rownames(completetable))-1)
    # add the accompagnying species
    acctable_c <- cbind(acctable[ , (max(cluster)+1) : (max(cluster)*2)],
                    acctable[ , 1 : (max(cluster)) ] )
    names(acctable_c) <- names(completetable)
    completetable_end <- rbind(completetable, acctable_c)


    results <- list("output" = "complete synoptic table, sorted by values of numeric input table and differential species character",
                    "species" = paste0("species with minimim value of", sep=" ", min1, " and their differentiating character"),
                    "samplesize" = tapply(rep(1,length(cluster)),cluster,sum),
                    "syntable" = completetable_end,
                    "others" = sort(rownames(syn1)[Hmisc::`%nin%` (rownames(syn1), rownames(completetable))]))
    }

  } else {stop("Sorting of synoptic table failed: wrong method entry. Check correct formula input")}

}



