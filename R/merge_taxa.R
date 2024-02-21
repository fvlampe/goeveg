#' Merge taxa with identical names
#' @description
#' The function offers a simple way to merge taxa with identical names in a vegetation table, e.g. due to necessary harmonization of the taxon level, to combine taxa of different layers or to remove duplicates.
#' The original cover-abundance scales are maintained. 
#'
#' @param vegtable Data frame with samples in columns and taxa in rows. Taxon names must be in the first column.
#' @param scale Cover-abundance scale(s) of data. Default is percentage cover (values between 0 and 100) (\code{"percentage"}). Alternatively it can be one of the included scales in \code{\link{scale_tabs}} (e.g. \code{"braun.blanquet"}) or a custom scale defined as data frame following the same format.
#' You can also provide a vector of the same length as the number of samples, to define individual scales (from \code{\link{scale_tabs}}) for each sample.
#' @param layers A logical evaluation to \code{TRUE} or \code{FALSE} \emph{(default)} indicating whether vegetation layers are to be included with layer information stored in the second column.
#' @param method Choice of method to combine cover. \code{"independent"} \emph{(default)} or \code{"exclusive"}. See details for methods.
#' @param clean_matrix A logical evaluation to \code{TRUE} or \code{FALSE} \emph{(default)} indicating whether taxa that do not occur (i.e. sum of abundances is zero) or samples that are empty (i.e. do not have any taxa present) are to be removed.
#' @param backtransform A logical evaluation to \code{TRUE} \emph{(default)} or \code{FALSE} indicating whether cover-abundance values should be back-transformed into their original cover-abundance scale or left as percentage cover.
#' 
#' @section Details:
#' The format required for this function is a data frame with samples in columns and taxa in rows, which corresponds to the export format of vegetation tables from Turboveg (Hennekens & Schaminee 2001).
#' Taxon names must be in the first column of the table (not row names as these do not allow duplicates).
#' If vegetation layers are to be included, layer information must be stored in the second column of the table. Taxa will then be merged only within the defined layers.
#'
#' If a cover-abundance scale from \code{\link{scale_tabs}} is defined, all cover-abundance values will be transformed into percentage cover for the merging process, and then back-transformed into the original cover-abundance scale.
#'
#' When combining cover values there are two possibilities following Fischer (2014) and Tichý & Holt (2011).
#' \itemize{
#'      \item \code{method = "independent" } \emph{(default)} assumes that covers can overlap and that they do so independently of each other (e.g. individuals of the lower layer are growing beneath individuals of the upper layer).
#'      In case of two layers the sum of the cover of layer 1 and layer 2 is the sum of the covers minus the overlap. ($p1 + p2- p1 * p2)
#'      Usually the choice when merging the same taxon from different (sub-)layers.
#'      \item \code{method = "exclusive" } assumes the covers are mutually exclusive, cover values will be simply summed up (e.g. individuals grow side by side).
#'      Usually the chouce when merging within layers, e.g. aggregating distinct taxa.
#'      }
#' Percentage cover values will eventually be truncated to 100%.
#' 
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
#' schedenveg.merged <- merge_taxa(schedenveg.t, scale = "percentage")
#'
#' @references
#' Fischer, H. S. (2015): On the combination of species cover values from different vegetation layers.
#' \emph{Applied Vegetation Science}, \strong{18}: 169–170. \doi{https://doi.org/10.1111/avsc.12130}
#'
#' Tichý, L. & Holt, J. (2011): JUICE. Program for management, analysis and classification of ecological data. Vegetation Science Group, Masaryk University Brno, CZ.
#' @author Friedemann von Lampe (\email{fvonlampe@uni-goettingen.de})
#' @export
#' @import stats utils


merge_taxa <- function(vegtable, scale = "percentage", layers = "FALSE", method = "independent", clean_matrix = FALSE, backtransform = TRUE) {
  
  narep <- FALSE
  # Check for "" values and replace by NA
  if (length(vegtable[vegtable == ""]) != 0) {
    vegtable[vegtable == ""] <- NA
    narep <- TRUE
  }
  # Check for NA values
  if (any(is.na(vegtable))) {
    vegtable[is.na(vegtable)] <- 0
    narep <- TRUE
  }
  if(narep == TRUE) print("NA and/or empty character values replaced by 0.")
  
  # create progress bar
  pb <- txtProgressBar(min = 0, max = 5, style = 3)
  
  # Transform cover values to percentage values ----
  vegcover <- data.frame(0)
  
  if(layers == TRUE) {
    vegcover <- vegtable[-c(1, 2)]
  } else {
    vegcover <- vegtable[-1]
  }
  
  percentage = FALSE
  
  setTxtProgressBar(pb, 1)
  
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
      # When vector of length of columns
      for(i in 1:ncol(vegcover)) {
        vegcover[, i] <- cov2per(vegcover[, i], scale = as.character(scale[i]))
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
  
  setTxtProgressBar(pb, 2)
  
  if(layers == TRUE) {
    vegtable <- data.frame(vegtable[,c(1:2)], vegcover)
    names(vegtable)[1:2] <- c("taxon", "layer")
  } else {
    vegtable <- data.frame(vegtable[ ,1], vegcover)
    names(vegtable)[1] <- c("taxon")
  }
  
  setTxtProgressBar(pb, 3)
  
  # Store duplicated taxa
  tax_dups = vegtable$taxon[duplicated(vegtable$taxon)]
  
  
  # Var A: Merge identical taxa (without layers) ----
  if(layers == FALSE) {
    
    # Calculate combined cover values based on selected method
    if(method == "independent") {
      
      vegtable[ , 2: ncol(vegtable)] <- vegtable[ , 2: ncol(vegtable)]/100
      vegtable <- aggregate(.~taxon, vegtable, function(x) (1 - prod(1 - x))*100)
      
    } else if(method == "exclusive") {
      vegtable <- aggregate(.~taxon, vegtable, sum)
    } else {
      stop("Unknown method")
    }
    
    # Truncate cover values > 100
    values <- vegtable[2:ncol(vegtable)]
    values[values > 100] <- 100
    vegtable[2:ncol(vegtable)] <- values
    
    setTxtProgressBar(pb, 4)
  }
  
  # Var B: Merge identical taxa (within layers) ----
  if(layers == TRUE) {
    
    # Calculate combined cover values based on selected method
    if(method == "independent") {
      
      vegtable[ , 3: ncol(vegtable)] <- vegtable[ , 3: ncol(vegtable)]/100
      vegtable <- aggregate(.~taxon + layer, vegtable, function(x) (1 - prod(1 - x))*100)
      
    } else if(method == "exclusive") {
      vegtable <- aggregate(.~taxon +  layer, vegtable, sum)
    } else {
      stop("Unknown method")
    }
    
    # Truncate cover values > 100
    values <- vegtable[3:ncol(vegtable)]
    values[values > 100] <- 100
    vegtable[3:ncol(vegtable)] <- values
    
    setTxtProgressBar(pb, 4)
  }
  
  
  if(layers == TRUE) {
    vegcover <- vegtable[-c(1, 2)]
  } else {
    vegcover <- vegtable[-1]
  }
  
  # Back-transform values
  if(percentage == FALSE && backtransform == TRUE) {
    
    if(!is.data.frame(scale)) {
      
      if (length(scale) == 1) {
        vegcover <- per2cov(vegcover, scale = scale)
      } else if(length(scale) > 1) {
        for(i in 1:ncol(vegcover)) {
          vegcover[, i] <- per2cov(vegcover[, i], scale = as.character(scale[i]))
        }
      }
      
    } else if(is.data.frame(scale)) {
      vegcover <- per2cov(vegcover, scale = scale)
    }
  }
  if(clean_matrix == TRUE) {
    # Calculate frequencies and richness
    freq <- apply(vegcover > 0, 1, sum)
    freq
    freq.0 <- length(freq[freq == 0])  
    
    rich <- apply(vegcover > 0, 2, sum)
    rich
    rich.0 <- length(rich[rich == 0]) 
    
    # Remove empty species and samples
    vegcover <- vegcover[freq > 0, rich > 0]
    vegtable <- vegtable[freq > 0, ]
  }
  
  if(layers == TRUE) {
    vegtable <- data.frame(vegtable[,c(1:2)], vegcover)
    names(vegtable)[1:2] <- c("taxon", "layer")
  } else {
    vegtable <- data.frame(vegtable[,1], vegcover)
    names(vegtable)[1] <- c("taxon")
  }
  
  setTxtProgressBar(pb, 5)
  
  # Close progress bar
  close(pb)
  
  
  if(length(tax_dups) > 0) {
    print(paste0(length(tax_dups), " duplicated taxa merged: ", toString(tax_dups)))
  }  else {
    print("Samples: No duplicated taxa within samples to combine.")
  }
  
  if(clean_matrix == TRUE) {
    print(paste0("Clean matrix: Removed ", freq.0, " species without occurence and ", 
                 rich.0, " sample(s) without any species"))
  }
  
  return(vegtable)
}