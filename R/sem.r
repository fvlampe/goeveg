#' Standard error of the mean (SEM)
#' @description  Compute the standard error of the mean (SEM). The SEM is the standard deviation of the sample-mean's estimate of a population mean.
#' It therefore describes the accuracy of the calculation of a sample's mean.
#' If \code{na.rm} is \code{TRUE} then missing values are removed before computation proceeds.
#' @param x a numeric vector
#' @param na.rm logical. Should missing values be removed?
#' @section Details:
#' The SEM of a zero-length vector (after removal of \code{NA}s if \code{na.rm = TRUE}) is not defined and gives an error.
#' The SEM of a length-one vector is \code{NA}.
#' @return
#' A numeric scalar -- the standard error of the mean.
#' @examples
#' ## Calculate mean and SEM for variable soil depth
#' mean(schedenenv$soil_depth)
#' sem(schedenenv$soil_depth)
#' @seealso \code{\link[stats]{sd}}
#' @export
#' @import stats

sem <- function(x, na.rm = FALSE)  {
  	sd(x, na.rm = na.rm)/sqrt(length(x))
}



