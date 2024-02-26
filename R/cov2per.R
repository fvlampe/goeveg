#' Conversion between cover-abundance codes and percentage cover
#' @description
#' These functions perform a conversion between cover-abundance codes from different survey scales and percentage cover. They can be applied on a matrix-like object or a single vector.
#'
#' \code{cov2per} performs conversion from cover-abundance to percentage values
#'
#' \code{per2cov} performs conversion from percentage to cover-abundance values
#'
#' You may chose between a set of scales with pre-defined conversion values, in \code{\link{scale_tabs}} or define your own conversion table following the same format.
#'
#' @param matrix Community data, a vector or matrix-like object with cover-abundance values
#' @param scale Cover-abundance scale (from \code{\link{scale_tabs}}) or dataframe with custom conversion table following the same format.
#' @param multiscale A logical evaluation to \code{FALSE} \emph{(default)} or \code{TRUE} indicating whether individual scales per sample should be used. For individual scales, you must provide a vector of the same length as the number of samples to the \code{scale} argument, defining the scale (from \code{\link{scale_tabs}}) for each sample (with samples as rows in matrix).
#'
#' @section Details:
#' If scales are not only cover-based but also abundance-based (e.g. Braun-Blanquet, Kohler) there are often no unique definitions about their conversion
#' into percentage cover. Therefore it is necessary to define and give reference to the applied conversion table.
#'
#' Cover-abundance codes are transformed into the mean percentage cover of their class.
#' For the conversion of percentage cover to cover-abundance codes, all values between the lower and upper border of the class are transformed into the corresponding code.
#'
#' The included cover-abundance scales and the associated conversion tables with references are explained in \code{\link{scale_tabs}}. On this site you also find definitions for defining a custom table.
#' 
#'
#' @examples
#' ## Conversion of species matrix with percentage cover to Braun-Blanquet values
#' schedenveg.bb <- per2cov(schedenveg)
#'
#' ## Conversion of only 10 samples to Londo values
#' schedenveg.londo <- per2cov(schedenveg[, 1:10], scale = "londo")
#'
#' ## Conversion of species matrix with Braun-Blanquet values to percentage cover
#' schedenveg.per <- cov2per(schedenveg.bb)
#'
#' @seealso \code{\link{scale_tabs}} for explanation and references of included conversion tables
#'
#' @author Friedemann von Lampe (\email{fvonlampe@uni-goettingen.de})
#' @import stats
#' @export
#' @returns A dataframe or vector containing the transformed data


cov2per <- function(matrix, scale = "braun.blanquet", multiscale = FALSE) {
  
  # Check for NA values
  if (any(is.na(matrix))) {
    matrix[is.na(matrix)] <- 0
    warning("NA values in matrix transformed into 0.")
  }
  
  if(multiscale == FALSE) {
    
    # Create cover objects
    cover <- matrix
    cover_new <- cover
    
    # Load cover-abundance scale
    if (is.data.frame(scale)) {
      
      if(ncol(scale) != 3) {
        stop("Custom conversion table need exactly 3 columns: code, cov_mean & cov_max.")
      }
      
      if(nrow(scale) < 1) {
        stop("No entries in custom conversion table.")
      }
      
      scale_tab <- rbind(c(0,0,0), scale)
      names(scale_tab) <- c("code", "cov_mean", "cov_max")
      
    } else if (any(scale == names(goeveg::scale_tabs))) {
      
      # Selection of pre-defined conversion tables
      scale_tab <- goeveg::scale_tabs[[scale]]
    } else {
      stop("Unknown cover-abundance scale.")
    }
    
    # Check if all values of matrix are defined
    testmatrix <- as.matrix(matrix)
    notdef <- length(testmatrix[!testmatrix %in% scale_tab$code])
    
    if(notdef > 0) {
      stop(paste0(notdef, " value(s) in matrix not defined in cover-abundance scale. E.g.: ",
                  testmatrix[!testmatrix %in% scale_tab$code][1]))
    }
    
    # Convert into percentage values
    for (i in 1:length(scale_tab$cov_mean)) {
      cover_new[cover == scale_tab$code[i]] <- scale_tab$cov_mean[i]
    }
    
    # Convert into numeric values
    if(is.null(dim(matrix))) {
      matrix <- as.numeric(cover_new)
    } else {
      matrix <- apply(cover_new, 2,  function(x) as.numeric(as.character(x)))
      # Following names are only assigned when still 2-dimensional
      if(!is.null(dim(matrix))) {
        matrix <- data.frame(matrix, row.names = row.names(cover))
        names(matrix) <- names(cover)
      }
    }
    
  } else if(multiscale == TRUE) {
    
    # Test if length of scale matches length of samples (rows)
    if(length(scale) != nrow(matrix)) {
      stop("Vector length for cover-abundance scales does not match number of samples (rows).")
    }
    
    # create progress bar
    pb <- txtProgressBar(min = 0, max = nrow(matrix), style = 3)
    
    # Make loop per sample (row)
    for(i in 1:nrow(matrix)) {
      
      setTxtProgressBar(pb, i)
      
      # Create cover objects
      cover <- matrix[i ,]
      cover_new <- cover
      
      # Load cover-abundance scale
      if (any(as.character(scale[i]) == names(goeveg::scale_tabs))) {
        
        # Selection of pre-defined conversion tables
        scale_tab <- goeveg::scale_tabs[[as.character(scale[i])]]
      } else {
        stop(paste0("Unknown cover-abundance scale in row ", i, ": ", as.character(scale[i])))
      }
      
      # Check if all values of matrix are defined
      testmatrix <- as.matrix(matrix[i ,])
      notdef <- length(testmatrix[!testmatrix %in% scale_tab$code])
      
      if(notdef > 0) {
        stop(paste0(notdef, " value(s) in sample ", i ," not defined in cover-abundance scale. E.g.: ",
                    testmatrix[!testmatrix %in% scale_tab$code][1]))
      }
      
      # Convert into percentage values
      for (l in 1:length(scale_tab$cov_mean)) {
        cover_new[cover == scale_tab$code[l]] <- scale_tab$cov_mean[l]
      }
      
      # Convert into numeric values
      matrix[i ,] <- apply(cover_new, 2,  function(x) as.numeric(as.character(x)))
      matrix[i ,] <- data.frame(matrix[i ,], row.names = row.names(cover))
      names(matrix[i ,]) <- names(cover)
      
    }
    
    # Close bar
    close(pb)
    
  } else {
    stop("Type unknown.")
  }
  
  # Check if still NA values (should not be possible)
  if (any(is.na(matrix))) {
    warning("NA values in final matrix. Check scales and values.")
  }
  
  return(matrix)
}

#' @export
#' @rdname cov2per






per2cov <- function(matrix, scale = "braun.blanquet", multiscale = FALSE) {
  
  # Check for NA values
  if (any(is.na(matrix))) {
    matrix[is.na(matrix)] <- 0
    warning("NA values in matrix transformed into 0.")
  }
  
  # tranform matrix to numeric; creates NA for non-numeric values
  if(is.null(dim(matrix))) {
    matrix <- as.numeric(matrix)
  } else {
    matrix <- apply(matrix, 2,  function(x) as.numeric(as.character(x)))
  }
  
  # Check if all values are numeric and between 0 and 100
  if(any(is.na(matrix))) {
    stop("Non-numeric percentage cover values.")
  }
  
  if(any(matrix < 0) | any(matrix > 100)) {
    stop("Percentage cover values need to be between 0 and 100.")
  }
  
  if(multiscale == FALSE) {
    
    # Create cover objects
    cover <- matrix
    cover_new <- cover
    
    # Load cover-abundance scale
    if (is.data.frame(scale)) {
      
      if(ncol(scale_tab) != 3) {
        stop("Custom conversion table need exactly 3 columns: code, cov_mean & cov_max.")
      }
      
      if(nrow(scale_tab) < 1) {
        stop("No entries in custom conversion table.")
      }
      
      scale_tab <- rbind(c(0,0,0), scale_tab)
      names(scale_tab) <- c("code", "cov_mean", "cov_max")
      
    } else if (any(scale == names(goeveg::scale_tabs))) {
      
      # Selection of pre-defined conversion tables
      scale_tab <- goeveg::scale_tabs[[scale]]
    } else {
      stop("Unknown cover-abundance scale")
    }
    
    
    # Clean scale table to keep only values for back-transformation
    scale_tab_clean <- na.omit(scale_tab)
    
    # Calculate back-transformed cover values
    # according to the class limits of the corresponding scale table
    for(l in 1:length(scale_tab_clean$cov_max)-1) {
      cover_new[(cover > scale_tab_clean$cov_max[l]) &
                  (cover <= scale_tab_clean$cov_max[l+1])] <- scale_tab_clean$code[l+1]
    }
    
    if(is.null(dim(cover_new))) {
      matrix <- cover_new
    } else {
      matrix <- as.data.frame(cover_new)
    }
    
    
  } else if(multiscale == TRUE) {
    
    # Test if length of scale matches length of samples (rows)
    if(length(scale) != nrow(matrix)) {
      stop("Vector length for cover-abundance scale does not match number of samples (rows).")
    }
    
    matrix_new <- matrix
    
    # create progress bar
    pb <- txtProgressBar(min = 0, max = nrow(matrix), style = 3)
    
    # Make loop per sample (row)
    for(i in 1:nrow(matrix)) {
      
      setTxtProgressBar(pb, i)
      
      # Create cover objects
      cover <- matrix[i, ]
      cover_new <- cover
      
      # Load cover-abundance scale
      if (any(as.character(scale[i]) == names(goeveg::scale_tabs))) {
        
        # Selection of pre-defined conversion tables
        scale_tab <- goeveg::scale_tabs[[as.character(scale[i])]]
      } else {
        stop(paste0("Unknown cover-abundance scale in row ", i, ": ", as.character(scale[i])))
      }
      
      # Clean scale table to keep only values for back-transformation
      scale_tab_clean <- na.omit(scale_tab)
      
      # Calculate back-transformed cover values
      # according to the class limits of the corresponding scale table
      for(l in 1:length(scale_tab_clean$cov_max)-1) {
        cover_new[(cover > scale_tab_clean$cov_max[l]) &
                    (cover <= scale_tab_clean$cov_max[l+1])] <- scale_tab_clean$code[l+1]
      }
      
      matrix_new[i, ] <- cover_new
    }
    
    matrix <- as.data.frame(matrix_new)
    
    # Close bar
    close(pb)
    
  } else {
    stop("Type unknown.")
  }
  
  # Check if still NA values (should not be possible)
  if (any(is.na(matrix))) {
    warning("NA values in final matrix. Check scales and values.")
  }
  
  return(matrix)
}






