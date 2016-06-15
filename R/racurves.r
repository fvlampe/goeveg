#' Multiple rank-abundance curves
#' @description This function draws multiple rank-abundance curves for selected samples into one diagram.
#' If you wish to draw a simple rank-abundance curve see \link{racurve}.
#' @param matrix Community data, a matrix-like object with samples in rows and species in columns. Rank-abundance curves are drawn for all selected rows (samples).
#' @param main The main title.
#' @param bw If set on \code{FALSE} the lines will be drawn in colours instead of black/white lines with different line types.
#' @section Details:
#' Rank abundance curves or Whittaker plots (see \cite{Whittaker 1965}) are used to display relative species abundance as biodiversity component.
#' They are a means to visualise species richness and species evenness.
#'
#' The axes of the diagram will be scaled according automatically.
#' As the line type is used to differenciate between samples, a maximum of 6 curves per diagram is feasible in black/white mode.
#' @examples
#' ## Draw multiple rank-abundance curves for selected samples
#' racurves(schedenveg[c(1,7,20,25), ])
#'
#' ## Draw multiple rank-abundance curves for selected samples with coloured lines
#' racurves(schedenveg[c(1,7,20,25), ], bw = FALSE)
#' @seealso \code{\link{racurve}}
#' @references Whittaker, R. H. (1965). Dominance and Diversity in Land Plant Communities: Numerical relations of species express the importance of competition in community function and evolution. \emph{Science} \strong{147 :} 250-260.
#' @author Friedemann Goral (\email{fgoral@gwdg.de})
#' @export

racurves <-  function(matrix, main = "Rank-abundance diagram", bw = TRUE) {
  if(!is.data.frame(matrix)) {
    matrix <- data.frame(matrix)
  }

  if(nrow(matrix) >= 1) {

    abund <- apply(matrix, 2, sum)
    sum(abund)
    rel.abund.plot <- prop.table(as.matrix(matrix), margin=1)
    rang.abund <- 1:ncol(matrix)
    rang.abund <- data.frame(rang.abund)

    ls <- nrow(matrix)

    for(i in 1:ls) {
      rang.abund[,i] <- sort(rel.abund.plot[i,], decreasing = T)
    }

    rang.abund[rang.abund==0] <- NA
    complete <- 0

    for(i in 1:ls) {
      complete[i] <- sum(complete.cases(rang.abund[,i]))
    }

    plot(rang.abund[,1], type="n", xlab="Rank", ylab="Relative Abundance",
         main = main,
         xlim=c(0,max(complete)), ylim =c(0,max(rang.abund, na.rm=T)))

    if(bw == TRUE) {

      for(i in 1:ls) {
        lines(rang.abund[,i], lty=i)
      }
      legend("topright", inset = 0.02, legend = row.names(matrix), lty = 1:ls,
             cex = 0.85, bty="n")
    } else {

       for(i in 1:ls) {
         lines(rang.abund[,i], col=1+i)
      }
      legend("topright", inset = 0.02, legend = row.names(matrix), lty = 1, col = 2:(ls+1),
             cex = 0.85, bty="n")

    }

  } else {
    print("FATAL: No samples in matrix")
  }
}
