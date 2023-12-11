#' Merge taxa with identical names
#' @description
#' The function offers a simple way to merge taxa with identical names in a vegetation table, e.g. due to necessary harmonization of the taxon level, to combine taxa of different layers or to remove duplicates.
#' The original cover-abundance scales are maintained.
#'
#' @param vegtable Data frame with samples in columns and taxa in rows. Taxon names must be in the first column.
#' @param scale Cover-abundance scale(s) of data. Default is percentage cover (values between 0 and 100) (\code{"percentage"}). Alternatively it can be one of the included scales in \code{\link{scale_tabs}} (e.g. \code{"braun.blanquet"}) or a custom scale defined as data frame following the same format.
#' You can also provide a vector of the same length as the number of samples, to define individual scales (from \code{\link{scale_tabs}}) for each sample.
#' @param layers A logical evaluation to \code{TRUE} or \code{FALSE} indicating whether vegetation layers are to be included with layer information stored in the second column.
#' @param method Choice of method to combine cover. \code{"independent"} (= default) or \code{"exclusive"}. See details for methods.
#'
#' @section Details:
#' The format required for this function is a data frame with samples in columns and taxa in rows, which corresponds to the export format of vegetation tables from Turboveg (Hennekens & Schaminee 2001).
#' Taxon names must be in the first column of the table (not row names as these do not allow duplicates).
#' If vegetation layers are to be included, layer information must be stored in the second column of the table. Taxa will then be merged only within the defined layers.
#'
#' If a cover-abundance scale from \code{\link{scale_tabs}} is defined, all cover-abundance values will be transformed into percentage cover for the merging process, and then back-transformed into the original cover-abundance scale.
#'
#' When combining cover values there are two possibilities following Fischer (2014) and Tichý & Holt (2011).
#' \description{
#'      \item[\code{method = "independent" }]{assumes that covers can overlap and that they do so independently of each other (e.g. individuals of the lower layer are growing beneath individuals of the upper layer).
#'      In case of two layers the sum of the cover of layer 1 and layer 2 is the sum of the covers minus the overlap. ($p1 + p2- p1 * p2)
#'      Usually the choice when merging the same taxon from different (sub-)layers.}
#'      \item[\code{method = "exclusive" }]{assumes the covers are mutually exclusive, cover values will be simply summed up (e.g. individuals grow side by side).
#'      Usually the chouce when merging within layers, e.g. aggregating distinct taxa.}
#'      }
#' Percentage cover values will eventually be truncated to 100%.
#'
#' @return
#' A data frame based on \code{vegtable} with merged taxa.
#'
#' @examples
#' ## Merge taxa with identical names without any layer information
#' # Transpose table to required format
#' schedenveg.t <- data.frame(species = names(schedenveg), t(schedenveg))
#' # Add two duplicated taxa
#' schedenveg.t <- rbind(schedenveg.t, schedenveg.t[c(55, 61), ])
#'
#' # Merge duplicated taxa using default 'independent' method
#' merge_taxa(schedenveg.t, scale = "percentage")
#'
#' @references
#' Fischer, H. S. (2015): On the combination of species cover values from different vegetation layers.
#' \emph{Applied Vegetation Science}, \strong{18}: 169–170. \doi{https://doi.org/10.1111/avsc.12130}
#'
#' Tichý, L. & Holt, J. (2011): JUICE. Program for management, analysis and classification of ecological data. Vegetation Science Group, Masaryk University Brno, CZ.
#' @author Friedemann von Lampe (\email{fvonlampe@uni-goettingen.de})
#' @export
#' @import stats utils




merge_taxa <- function(vegtable, scale = "percentage", layers = "FALSE", method = "independent") {

  # Transform cover values to percentage values ----
  vegcover <- data.frame(0)

  if(layers == TRUE) {
    vegcover <- vegtable[-c(1, 2)]
  } else {
    vegcover <- vegtable[-1]
  }

  percentage = FALSE

  if(!is.data.frame(scale)) {

    if (length(scale) == 1) {
      if (scale == "percentage") {
        # Check percentage values
        percentage = TRUE

        # Check if all values are numeric and between 0 and 100
        if(any(is.na(as.numeric(unlist(vegcover))))) {
          stop("Non-numeric percentage cover values.")
        }
        if(any(unlist(vegcover) < 0) |
           any(unlist(vegcover) > 100)) {
          stop("Percentage cover values need to be between 0 and 100.")
        }
      } else {
        # Convert non-percentage values based on scale table
        vegcover <- cov2per(vegcover, scale = scale)
      }

    } else if(length(scale) == ncol(vegcover)) {
      # When vector of length
      for(i in 1:ncol(vegcover)) {
        vegcover[, i] <- cov2per(vegcover[, i], scale = scale[i])
      }
    } else {
      stop("Vector length for cover-abundance scale does not match number of samples.")
    }

  } else if(is.data.frame(scale)) {
    # Convert non-percentage values based on dataframe (same function)
    vegcover <- cov2per(vegcover, scale = scale)
  }  else {
    stop("Unknown cover-abundance scale.")
  }

  if(layers == TRUE) {
    vegtable <- cbind(vegtable[,c(1:2)], vegcover)
  } else {
    vegtable <- cbind(vegtable[,1], vegcover)
  }

  # Convert zero to NA to keep only existing taxa in long table
  if(0 %in% as.matrix(vegtable)) {
    vegtable[vegtable == 0] <- NA
  }

  # Conversion to long table  ----
  if(layers == TRUE) {
    tab_long <- na.omit(cbind(vegtable[1:2], stack(vegtable[3:ncol(vegtable)]), row.names = NULL))
    names(tab_long) <- c("taxon", "layer", "cover", "site")
  } else {
    tab_long <- na.omit(cbind(vegtable[1], stack(vegtable[2:ncol(vegtable)]), row.names = NULL))
    names(tab_long) <- c("taxon", "cover", "site")
  }

  # Split long table to list with releves
  tab_long_list <- split(tab_long, tab_long$site)


  # Var A: Merge identical taxa (without layers) ----
  if(layers == FALSE) {

    tax_dups = character(0)

    # Run through each releve in list
    for(l in 1:length(tab_long_list)) {

      # Check each releve for duplicated taxa
      tax_dups <- tab_long_list[[l]][duplicated(tab_long_list[[l]]$taxon), ]$taxon

      # Run loop for every duplicated taxon
      if(length(tax_dups) != 0) {

        for(n in 1:length(tax_dups)) {

          # Select duplicated taxa
          tax_list <- tab_long_list[[l]][tab_long_list[[l]]$taxon == tax_dups[n], ]

          # Calculate combined cover values based on selected method
          if(method == "independent") {
            values <- tax_list$cover/100
            cov_comb <- (1 - prod(1 - values))*100
            cov_comb
          } else if(method == "exclusive") {
            cov_comb <- sum(tax_list$cover)
          } else {
            stop("Unknown method")
          }

          # Remove duplicated taxa from list item
          tab.a <- tab_long_list[[l]][tab_long_list[[l]]$taxon != tax_dups[n], ]
          tab.a

          # Add merged taxa to list item
          tab.b <- data.frame(
            taxon = tax_dups[n],
            cover = cov_comb,
            site = tab_long_list[[l]]$site[1])
          tab.b
          tab_long_list[[l]] <- rbind(tab.a, tab.b)
        }
        print(paste0(length(tax_dups), " taxa combined in sample ", names(tab_long_list[l]), ": ", toString(tax_dups)))
      }
    }
  }


  # Var B: Merge identical taxa (within layers) ----
  if(layers == TRUE) {

    tax_dups = character(0)

    # Run through each releve in list
    for(l in 1:length(tab_long_list)) {

      tab_long_list[[l]]$specieslayer <- paste(tab_long_list[[l]]$taxon, tab_long_list[[l]]$layer)

      # Check each releve for duplicated taxon-layers combinations
      tax_dups <- tab_long_list[[l]][duplicated(tab_long_list[[l]]$specieslayer), ]$specieslayer

      # Run loop for every duplicated taxon
      if(length(tax_dups) != 0) {

        for(n in 1:length(tax_dups)) {

          # Select duplicated taxa
          tax_list <- tab_long_list[[l]][tab_long_list[[l]]$specieslayer == tax_dups[n], ]

          # Calculate combined cover values based on selected method
          if(method == "independent") {
            values <- tax_list$cover/100
            cov_comb <- (1 - prod(1 - values))*100
            cov_comb
          } else if(method == "exclusive") {
            cov_comb <- sum(tax_list$cover)
          } else {
            stop("Unknown method")
          }

          # Remove duplicated taxa from list item
          tab.a <- tab_long_list[[l]][tab_long_list[[l]]$specieslayer != tax_dups[n], ]

          # Add merged taxa to list item
          tab.b <- data.frame(
            taxon = tax_list$taxon[1],
            layer = tax_list$layer[1],
            cover = cov_comb,
            site = tax_list$site[1],
            specieslayer = tax_list$specieslayer[1])

          tab_long_list[[l]] <- rbind(tab.a, tab.b)
        }
        print(paste0(length(tax_dups), " taxa combined in sample ", names(tab_long_list[l]), ": ", toString(tax_dups)))
      }

      # Remove combined taxon-layer colum
      tab_long_list[[l]] <- tab_long_list[[l]][, -5]
    }
  }


  # Backtransform into original cover abundance scale ----

  # Unlist releves
  tab_long_new <- do.call(rbind, tab_long_list)

  # Truncate cover values > 100
  tab_long_new$cover[tab_long_new$cover > 100] <- 100

  # Conversion to wide table ----
  if(layers == TRUE) {
    tab_new <- reshape(tab_long_new, timevar = 'site',idvar = c('taxon', 'layer'), direction = 'wide',
                       v.names = 'cover')
  } else {
    tab_new <- reshape(tab_long_new, timevar = 'site',idvar = 'taxon', direction = 'wide',
                       v.names = 'cover')
  }

  names(tab_new) <- gsub("cover.", "", names(tab_new))
  row.names(tab_new) <- 1:nrow(tab_new)

  # Convert NA to zero
  tab_new[is.na(tab_new)] <- 0

  # Back-transform values
  if(percentage == FALSE) {
    if(layers == TRUE) {
      vegcover <- tab_new[-c(1, 2)]
    } else {
      vegcover <- tab_new[-1]
    }

    if(!is.data.frame(scale)) {

      if (length(scale) == 1) {
        vegcover <- per2cov(vegcover, scale = scale)
      } else if(length(scale) > 1) {
        for(i in 1:ncol(vegcover)) {
          vegcover[, i] <- per2cov(vegcover[, i], scale = scale[i])
        }
      }

    } else if(is.data.frame(scale)) {
      vegcover <- per2cov(vegcover, scale = scale)
    }

    if(layers == TRUE) {
      tab_new <- cbind(tab_new[,c(1:2)], vegcover)
      names(tab_new)[1:2] <- c("Taxon name", "Layer")
    } else {
      tab_new <- cbind(tab_new[,1], vegcover)
      names(tab_new)[1] <- "Taxon name"
    }
  }


  return(tab_new)
}
