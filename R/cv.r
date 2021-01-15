#' Coefficient of variation (CV)
#' @description
#' Compute the coefficient of variation (CV). The CV, also known as relative standard deviation (RSD), is a standardized measure of dispersion of a probability distribution or frequency distribution.
#' It is defined as the ratio of the standard deviation to the mean and is often expressed as a percentage.
#' In contrast to the standard deviation, it enables comparison between datasets as the CV is independent of the unit in which the measurement has been taken.
#' If \code{na.rm} is \code{TRUE} then missing values are removed before computation proceeds.
#' @param x a numeric vector
#' @param na.rm logical. Should missing values be removed?
#' @section Details:
#' The coefficient of variation (CV) should be computed only for data measured on a ratio scale, as these are the measurements that can only take non-negative values.
#' The CV may not have any meaning for data on an interval scale.
#'
#' According to \cite{Dormann 2013} CV-values below 0.05 (5\%) indicate very high precision of the data, values above 0.2 (20\%) low precision.
#' However, this is considered as a rule of thumb. In studies of highly variable systems (e.g. some ecological studies) CV values above 1 may occur.
#'
#' The CV of a zero-length vector (after removal of \code{NA}s if \code{na.rm = TRUE}) is not defined and gives an error.
#' If there is only a single value, \code{sd} is \code{NA} and \code{cv} returns \code{NA}.
#' @examples
#' ## Calculate CV for variable soil depth
#' cv(schedenenv$soil_depth)
#' @seealso \code{\link[stats]{sd}}
#' @references Dormann, C. (2013). Parametrische Statistik. Verteilungen, maximum likelihood und GLM in R. \emph{Springer}.
#' @references "What is the difference between ordinal, interval and ratio variables? Why should I care?" \emph{GraphPad Software Inc}. \url{https://www.graphpad.com/support/faqid/1089/}.
#' @export
#' @import stats

cv <- function(x, na.rm = FALSE)  {
  sd(x, na.rm = na.rm)/mean(x, na.rm = na.rm)
}

