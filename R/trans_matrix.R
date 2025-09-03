#' Transpose species matrix
#' @description
#' The function transposes a species matrix, while preserving correct species and sample names. 
#' The new column names must be stored as row names of the data frame. They may also be stored in the first column, when chosing the argument \code{row.names = F}. 
#' 
#' @param matrix Community data, a data frame. 
#' @param row.names A logical evaluation indicating whether the new column names are stored as row names of the data frame \code{TRUE} \emph{(default)} or in the first column \code{FALSE}.
#' @param rmchar A logical evaluation indicating whether the first character of the original column names should be removed (default: \code{FALSE}).
#' 
#' @return
#' A transposed data frame.
#' 
#' @section Details:
#' Sometimes vegetation data is organized as a data frame with samples in columns and taxa in rows, with taxon names stored in the first column, e.g. as result of the function \code{\link{merge_taxa}}.
#' In this case you can use \code{row.names = F} to directly convert this species matrix into a statistically analyzable format, e.g. with \code{vegan}.
#' 
#' If your dataframe contains prepended “X” to each header due to numbered samples, you can use \code{rmchar = TRUE} to remove the first character of the column names during transposing. 
#' (You may also avoid this problem at all by using \code{check.names = FALSE} when loading the data in \code{\link[utils]{read.table}})
#'
#' @examples
#' # Transpose species matrix
#' schedenveg.trans <- trans_matrix(schedenveg)
#' 
#' @author Friedemann von Lampe
#' @export
#' @import stats utils


trans_matrix <- function(matrix, row.names = T, rmchar = FALSE) {
  
  if(row.names == T) {
    colnames <- names(matrix)
    if(rmchar == TRUE) colnames <- sub('.', '', colnames) 
    rownames <- row.names(matrix)                 # Row names as names
    
    matrix.trans <- data.frame(t(matrix), row.names = colnames)
    names(matrix.trans) <- rownames
  } else {
    colnames <- names(matrix[-1])
    if(rmchar == TRUE) colnames <- sub('.', '', colnames) 
    rownames <- matrix[[1]]                    # First column as names
    
    matrix.trans <- data.frame(t(matrix)[-1,], row.names = colnames)
    names(matrix.trans) <- rownames
  }
  
  # Replace empty character values
  matrix.trans[matrix.trans == ""] <- NA
  
  return(matrix.trans)
}