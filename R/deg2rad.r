#' Conversion between degrees and radians
#' @description
#' \code{deg2rad} performs conversion from degrees to radians.
#'
#' \code{rad2deg} performs conversion from radians to degrees.
#' @section Details:
#' Radians and degrees are both units used for measuring angles.
#'
#' A degree is a measure of angle equal to 1/360th of a revolution, or circle.
#' A radian is the measurement of angle equal to the length of an arc divided by the radius of the circle or arc.
#' A circle is comprised of 2*pi radians, which is the equivalent of 360 degrees.
#'
#' A common application in ecological studies is the conversion of measured exposition (in degrees) of plots into statistically meaningful measures, such as the north value or the east value.
#' For this, the cosine (for northness) or sine (for eastness) is applied to the radian of the exposition.
#'
#' @examples
#' ## Covert the value pi to degrees
#' rad2deg(pi)
#'
#' # Calculate north and east values based on exposition measured in degrees
#' north <- cos(deg2rad(schedenenv$exp))
#' east <- sin(deg2rad(schedenenv$exp))
#'
#'
#' @references BIPM (2019): The International System of Units (SI). Bureau international des poids et mesures, ninth edition. \url{http://www.bipm.org/en/si/si_brochure/}, ISBN 978-92-822-2272-0
#' @import stats
#' @export
#' @param x a numeric vector
#' @returns a numeric vector the same length as \code{x}

deg2rad <- function(x) {(x * pi) / (180)}

#' @export
#' @rdname deg2rad
rad2deg <- function(x) {(x * 180) / (pi)}



