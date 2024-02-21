#' Transpose species matrix
#' @description
#' The function transposes a species matrix, while preserving correct species and sample names. 
#' The new column names must be stored as row names of the data frame. They may also be stored in the first column, when chosing the argument \code{row.names = F}. 
#'
#' @param matrix Community data, a data frame. 
#' @param row.names A logical evaluation indicating whether the new column names are stored as row names of the data frame \code{TRUE} \emph{(default)} or in the first column \code{FALSE}.
#'
#' @return
#' A transposed data frame.
#'
#' @examples
#' # Transpose species matrix
#' schedenveg.trans <- trans_matrix(schedenveg)
#' 
#' @author Friedemann von Lampe (\email{fvonlampe@uni-goettingen.de})
#' @export
#' @import stats utils


trans_matrix <- function(matrix, row.names = T) {
  if(row.names == T) {
    
    colnames <- names(matrix)
    rownames <- row.names(matrix)                 # Row names as names
    
    matrix.trans <- data.frame(t(matrix), row.names = colnames)
    names(matrix.trans) <- rownames
  } else {
    colnames <- names(matrix[-1])
    rownames <- matrix[[1]]                    # First column as names
    
    matrix.trans <- data.frame(t(matrix)[-1,], row.names = colnames)
    names(matrix.trans) <- rownames
  }
  return(matrix.trans)
}