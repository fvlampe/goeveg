#' Synoptic tables and calculation of cluster-wise frequencies, fidelity and
#' differential species character
#'
#' @description Synoptic tables are a tool for the visualization and interpretation of previously
#' defined plant species groups (clusters), e.g. from cluster analysis or classification methods.
#' They help to determine characteristic patterning of species occurrences in plant communities
#' by calculating cluster-wise percentage or absolute frequencies, mean/median cover values, fidelity
#' (phi) or differential species character.
#' \code{syntable} function calculates an unordered synoptic table for plant community analysis, using
#' an input species-sample dataframe and a numeric vector of cluster identity input.
#' The unordered output table can be sorted automatically with \code{\link[goeveg]{synsort}} function
#' in this package.
#'
#' @param spec Species matrix or dataframe with species in columns and samples in rows.
#' Values must be numeric, with point "." as decimal character, or integer.  Missing values, NA or NaN
#' are not allowed. Species and sample names must be defined as column- and rownames, respectively.
#' @param cluster Integer vector with classification cluster identity. Ensure matching order of
#' cluster identity and samples in dataframe for correct allocation of cluster numbers to samples.
#' @param abund Data input type. Define whether input species matrix or dataframe is presence/absence
#' data (\code{abund = "freq"}) or percentage cover (\code{abund = "perc"}, default).
#' @param type Type of synoptic table output \code{type = c("percfreq", "totalfreq", "mean",
#' "median", "diffspec", "phi")}. See Details.
#'
#' @section Details:
#' For synoptic table calculation, six types are available.
#'   \itemize{\item{\code{type = "percfreq" }}{Default, creates a percentage frequency table}
#'   \item{\code{type = "totalfreq" }}{Creates an absolute frequency table}
#'   \item{\code{type = "mean" }}{Calculates mean of species values given in \code{spec} per cluster}
#'   \item{\code{type = "median" }}{Calculates median of species values given in \code{spec} per
#'    cluster}
#'   \item{\code{type = "diffspec" }}{Calculates differential character of species according to
#'    Tsiripidis et al. 2009, with resulting character p = positive, n = negative, pn = positive-
#'    negative or no differential character (-). Consider that differential character is always
#'    restricted to some and not necessarily all of the other units, thus considering percentage
#'    frequency is essential for correct interpretation of the diagnostic species character.}
#'   \item{\code{type = "phi" }}{Calculates fidelity measure phi (algorithm basing on Sokal & Rohlf
#'    1995, Bruelheide 2000). Values are ranging between -1 and 1 with high values near 1 indicating
#'    high fidelity.}}
#'
#' For sorting the output synoptic table, use \code{\link{synsort}} function, providing several
#' options.
#'
#' @return
#' The function returns a list of result components.
#'   \item{\code{$syntable }}{unordered synoptic table for given species and clusters}
#'   \item{\code{$samplesize }}{total samples in clusters}
#'
#'Additionally for differential species character calculation:
#'   \item{\code{$onlydiff }}{Synoptic table only with differential species}
#'   \item{\code{$others }}{List of non-differential species}
#'   \item{\code{$differentials }}{Lists differential species for each cluster}
#'
#'
#' @references
#' Bruelheide, H. (2000): A new measure of fidelity and its application to defining species groups.
#'  \emph{Journal of Vegetation Science} \strong{11}: 167-178. \doi{https://doi.org/10.2307/3236796}
#'
#' Chytry, M., Tichy, L., Holt, J., Botta-Dukat, Z. (2002): Determination of diagnostic species with
#'  statistical fidelity measures. \emph{Journal of Vegetation Science} \strong{13}: 79-90. \doi{https://doi.org/10.1111/j.1654-1103.2002.tb02025.x}
#'
#' Sokal, R.R. & Rohlf, F.J. (1995): Biometry. 3rd edition Freemann, New York.
#'
#' Tsiripidis, I., Bergmeier, E., Fotiadis, G. & Dimopoulos, P. (2009): A new algorithm for the
#' determination of differential taxa. \emph{Journal of Vegetation Science} \strong{20}: 233-240. \doi{https://doi.org/10.1111/j.1654-1103.2009.05273.x}
#'
#' @author Jenny Schellenberg (\email{jschell@gwdg.de})
#' @seealso \code{\link{synsort}}
#' @examples
#' ## Synoptic table of Scheden vegetation data
#' library(cluster)
#' pam1 <- pam(schedenveg, 4)  # PAM clustering with 4 clusters output
#'
#' ## 1) unordered synoptic percentage frequency table
#' unordered <- syntable(schedenveg, pam1$clustering, abund = "perc",
#'                       type = "percfreq")
#' unordered                   # view results
#'
#' ## 2) differential species analysis
#' differential <- syntable(schedenveg, pam1$clustering, abund = "perc",
#'                          type = "diffspec")
#' # show complete table with differential character of species
#' differential$syntable
#' # list differential species for second cluster
#' differential$differentials$group2
#'
#' ## 3) Synoptic table with phi fidelity
#' phitable <- syntable(schedenveg, pam1$clustering, abund = "perc",
#'                      type = "phi")
#' phitable
#' @export


syntable <- function(spec, cluster, abund = "perc", type = "percfreq") {

  spec <- as.data.frame(apply(spec,2,as.numeric))
  group <- sort(unique(cluster))
  syntab <- data.frame(matrix(NA, nrow=length(names(spec)), ncol=max(cluster)))
  names(syntab) <- c(seq(1, max(cluster),1))
  rownames(syntab) <- names(spec)

  if (abund != "perc" & abund != "freq") {
    stop("Check entry. Argument 'abund' must be either percentages ('perc') or
         presence/absence data ('freq')")
  } else {
    if (type=="totalfreq")   {
     for (i in 1:max(cluster)) {syntab[,i] <- apply(spec[cluster==i,]>0,2,sum)}
     samplesize <- tapply(rep(1,length(cluster)),cluster,sum)
     result <- list("syntable" = syntab,
                   "samplesize" = samplesize)

     } else if (type=="percfreq") {
    samplesize <- tapply(rep(1,length(cluster)),cluster,sum)
    for (i in 1:max(cluster)) {
      syntab[,i] <- round(apply(spec[cluster==i,]>0,2,sum)*(100/samplesize[i]), digits=0)}
    result <- list("syntable" = syntab, "samplesize" = samplesize)

  } else if (type=="mean" & abund=="perc") {
    samplesize <- tapply(rep(1,length(cluster)),cluster,sum)
    for (i in 1:max(cluster)) { syntab[,i] <- round(apply(spec[cluster==i,],2,mean), digits=0) }
    result <- list("syntable" = syntab, "samplesize" = samplesize)

  } else if (type=="mean" & abund=="freq") {
    stop("Can?t calculate mean cover in clusters with frequency values.")

  } else if (type=="median" & abund=="cover") {
    for (i in 1:max(cluster)) { syntab[,i] <- round(apply(spec[cluster==i,],2,median), digits=0)   }
    result <- list("syntable" = syntab, "samplesize" = samplesize)

  } else if (type=="median" & abund=="freq") {
    stop("Can?t calculate median cover in clusters with frequency values.")

  } else if (type=="median" & abund=="perc") {
    samplesize <- tapply(rep(1,length(cluster)),cluster,sum)
    for (i in 1:max(cluster)) { syntab[,i] <- round(apply(spec[cluster==i,],2,median), digits=0) }
    result <- list("syntable" = syntab, "samplesize" = samplesize)

  } else if (type=="diffspec") {
    samplesize <- tapply(rep(1,length(cluster)),cluster,sum)
    for (i in 1:max(cluster)) {
      syntab[,i] <- round(apply(spec[cluster==i,]>0,2,sum)*(100/samplesize[i]), digits=0)}
    syn <- syntab

    small <- data.frame(matrix(NA, nrow=length(rownames(syn)), ncol=length(group)))
    names(small) <- group
    for (i in 1:length(rownames(syn))) { small[i,] <- sort(as.numeric(syn[i,])) }

    Avsmall <- data.frame(matrix(NA, nrow=length(rownames(syn)), ncol=length(group)))
    Avsmall[,1] <- small[,1]
    for (i in 1:length(rownames(syn))) { for ( k in 2:length(group)) {
      Avsmall[i,k] <- mean(c(as.numeric(small[i,1:k]))) }   }

    CondAvsmall <- data.frame(matrix(NA, nrow=length(rownames(syn)), ncol=length(group)))
    CondAvsmall[,1] <- rep("", length(rownames(syn)))
    for ( k in 1:length(rownames(syn))) {
      if (small[k,2] >= 2*Avsmall[k,1]+20) {CondAvsmall[k,2] <- "1"}
      else {CondAvsmall[k,2] <- "0"}   }
    for ( i in 3:length(group)) { for ( k in 1:length(rownames(syn))) {
      if (CondAvsmall[k,i-1]=="1") {CondAvsmall[k,i] <- "1"}
      else if (small[k,i] >= 2*Avsmall[k,i-1]+20) {CondAvsmall[k,i] <- "1"}
      else {CondAvsmall[k,i] <- "0"}  } }

    Avsmall2 <- data.frame(matrix(NA, nrow=length(rownames(syn)), ncol=length(group)))
    Avsmall2[,1] <- Avsmall[,1]
    for ( i in 2:length(group)) { for ( k in 1:length(rownames(syn))) {
      if (CondAvsmall[k,i]=="1") {Avsmall2[k,i] <- Avsmall2[k,i-1]}
      else {Avsmall2[k,i] <- Avsmall[k,i]}  } }

    large <- data.frame(matrix(NA, nrow=length(rownames(syn)), ncol=length(group)))
    names(large) <- group
    for (i in 1:length(rownames(syn))) { large[i,] <- sort(as.numeric(syn[i,]), decreasing = T) }

    Avlarge <- data.frame(matrix(NA, nrow=length(rownames(syn)), ncol=length(group)))
    Avlarge[,1] <- large[,1]
    for (i in 1:length(rownames(syn))) { for ( k in 2:length(group)) {
      Avlarge[i,k] <- mean(c(as.numeric(large[i,1:k]))) }   }

    CondAvlarge <- data.frame(matrix(NA, nrow=length(rownames(syn)), ncol=length(group)))
    CondAvlarge[,1] <- rep("", length(rownames(syn)))
    for ( k in 1:length(rownames(syn))) {
      if (large[k,2] <= (Avlarge[k,1]/2)-10) {CondAvlarge[k,2] <- "1"}
      else {CondAvlarge[k,2] <- "0"}    }
    for ( i in 3:length(group)) { for ( k in 1:length(rownames(syn))) {
      if (CondAvlarge[k,i-1]=="1") {CondAvlarge[k,i] <- "1"}
      else if (large[k,i] <= (Avlarge[k,i-1]/2)-10) {CondAvlarge[k,i] <- "1"}
      else {CondAvlarge[k,i] <- "0"}   } }

    Avlarge2 <- data.frame(matrix(NA, nrow=length(rownames(syn)), ncol=length(group)))
    Avlarge2[,1] <- Avlarge[,1]
    for ( i in 2:length(group)) {  for ( k in 1:length(rownames(syn))) {
      if (CondAvlarge[k,i]=="1") {Avlarge2[k,i] <- Avlarge2[k,i-1]}
      else {Avlarge2[k,i] <- Avlarge[k,i]}    } }

    CondAvsmall2 <- data.frame(matrix(NA, nrow=length(rownames(syn)), ncol=length(group)))
    CondAvsmall2[,1] <- rep("", length(rownames(syn)))
    for ( i in 2:length(group)) {for ( k in 1:length(rownames(syn))) {
      if (small[k,i] > (Avlarge2[k,i-1]/2)-10) {CondAvsmall2[k,i] <- "1"}
      else {CondAvsmall2[k,i] <- "0"}  } }

    Avsmall3 <- data.frame(matrix(NA, nrow=length(rownames(syn)), ncol=length(group)))
    Avsmall3[,1] <- Avsmall2[,1]
    for ( i in 2:length(group)) { for ( k in 1:length(rownames(syn))) {
      if (CondAvsmall2[k,i]=="1") {Avsmall3[k,i] <- Avsmall3[k,i-1]}
      else {Avsmall3[k,i] <- Avsmall2[k,i]}  } }

    CondAvlarge2 <- data.frame(matrix(NA, nrow=length(rownames(syn)), ncol=length(group)))
    CondAvlarge2[,1] <- rep("", length(rownames(syn)))
    for ( i in 2:length(group)) { for ( k in 1:length(rownames(syn))) {
      if (large[k,i] < 2*Avsmall2[k,i-1]+20) {CondAvlarge2[k,i] <- "1"}
      else {CondAvlarge2[k,i] <- "0"} } }

    Avlarge3 <- data.frame(matrix(NA, nrow=length(rownames(syn)), ncol=length(group)))
    Avlarge3[,1] <- Avlarge2[,1]
    for ( i in 2:length(group)) {  for ( k in 1:length(rownames(syn))) {
      if (CondAvlarge2[k,i]=="1") {Avlarge3[k,i] <- Avlarge3[k,i-1]}
      else {Avlarge3[k,i] <- Avlarge2[k,i]}  } }

    Small_Avsmall3 <- data.frame(matrix(NA, nrow=length(rownames(syn)), ncol=length(group)))
    for ( i in 1:length(group)) { for ( k in 1:length(rownames(syn))) {
      if (syn[k,1] == "") {Small_Avsmall3[k,i] <- ""}
      else if (syn[k,i]>=2*max(Avsmall3[k,])+20) {Small_Avsmall3[k,i] <- "p"}
      else if (syn[k,i] < 2*max(Avsmall3[k,])+20) {Small_Avsmall3[k,i] <- "n"}
      else {stop("Calculation of differential species failed: Check synoptic table input")}  } }

    Large_AvLarge3 <- data.frame(matrix(NA, nrow=length(rownames(syn)), ncol=length(group)))
    for ( i in 1:length(group)) {  for ( k in 1:length(rownames(syn))) {
      if (syn[k,1] == "") {Large_AvLarge3[k,i] <- ""}
      else if (syn[k,i] >  min(Avlarge3[k,])/2-10) {Large_AvLarge3[k,i] <- "p"}
      else if (syn[k,i] <= min(Avlarge3[k,])/2-10) {Large_AvLarge3[k,i] <- "n"}
      else {stop("Calculation of differential species failed: Check synoptic table input")}   } }

    Sum_cond <-  data.frame(matrix(NA, nrow=length(rownames(syn)), ncol=length(group)))
    for ( i in 1:length(group)) { for ( k in 1:length(rownames(syn))) {
      Sum_cond[k,i] <- paste0(Small_Avsmall3[k,i], Large_AvLarge3[k,i])   } }

    Results <- data.frame(matrix(NA, nrow=length(rownames(syn)), ncol=length(group)))
    for ( i in 1:length(group)) { for ( k in 1:length(rownames(syn))) {
      if (Sum_cond[k,i] == "pp") {Results[k,i] <- "p"}
      else if (Sum_cond[k,i] == "nn") {Results[k,i] <- "n"}
      else if (Sum_cond[k,i] == "pn") {Results[k,i] <- "pn"}
      else {Results[k,i] <- "-"} } }
    rownames(Results) <- rownames(syn)
    names(Results) <- group

    subresult <- data.frame(matrix(NA, nrow=length(rownames(syn)), ncol=length(group)))
    for ( i in 1:length(group)) { for ( k in 1:length(rownames(Results))) {
      if (Results[k,i] == "-") {subresult[k,i] <- "0"}
      else if (Results[k,i] == "n") {subresult[k,i] <- "1"}
      else if (Results[k,i] == "p") {subresult[k,i] <- "1"}
      else if (Results[k,i] == "pn") {subresult[k,i] <- "1"}
      else {subresult[k,i] <- "0"} } }
    subresult <- as.data.frame(apply(subresult,2,as.numeric))

    onlydiff <- Results[apply(subresult,1,sum)>0,]
    others <- rownames(Results[apply(subresult,1,sum)==0,])

    positivespec <- data.frame(as.matrix(NA, nrow=500, ncol=200))
    for ( i in 1:length(group)) {
      if (length(rownames(Results[Results[,i]=="p",]) > 0)) {
        a <- rownames(Results[Results[,i]=="p",])
        for ( k in 1:length(a)) {
          positivespec[i,k] <- rownames(Results[Results[,i]=="p",])[k] }
      } else {positivespec[i,1] <- NA} }
    names(positivespec) <- rep("", length(positivespec[1,]))

    negativespec <- data.frame(as.matrix(NA, nrow=500, ncol=200))
    for ( i in 1:length(group)) {
      if (length(rownames(Results[Results[,i]=="n",]) > 0)) {
        a <- rownames(Results[Results[,i]=="n",])
        for ( k in 1:length(a)) {
          negativespec[i,k] <- rownames(Results[Results[,i]=="n",])[k] }
      }  else {negativespec[i,1] <- NA}   }
    names(negativespec) <- rep("", length(negativespec[1,]))

    positivenegativespec <- data.frame(as.matrix(NA, nrow=500, ncol=200))
    for ( i in 1:length(group)) {
      if (length(rownames(Results[Results[,i]=="pn",]) > 0)) {
        a <- rownames(Results[Results[,i]=="pn",])
        for ( k in 1:length(a)) {
          positivenegativespec[i,k] <- rownames(Results[Results[,i]=="pn",])[k] }
      } else {positivenegativespec[i,1] <- NA}   }
    names(positivenegativespec) <- rep("", length(positivenegativespec[1,]))

    diffspeclist <- list()
    for ( i in 1:max(cluster)) {
      pos <- positivespec[i,]
      pos <- pos[!is.na(pos)]
      if (length(pos) > 0) {pos=pos}
      else {pos <- noquote("no positive diff species")}
      neg <- negativespec[i,]
      neg <- neg[!is.na(neg)]
      if (length(neg) > 0) {neg=neg}
      else {neg <- noquote("no negative diff species")}
      posneg <- positivenegativespec[i,]
      posneg <- posneg[!is.na(posneg)]
      if (length(posneg) > 0) {posneg=posneg}
      else {posneg <- noquote("no positive/negative diff species")}
      name <- paste("group", i, sep="")
      tmp <- list("positive diff" = pos, "negative diff" = neg, "positive/negative diff" = posneg)
      diffspeclist[[name]] <- tmp }

    result <- list("syntable" = Results,
                   "onlydiff" = onlydiff,
                   "others" = others,
                   "samplesize" = samplesize,
                   "differentials" = diffspeclist)

    } else if (type == "phi") {
      for (i in 1:max(cluster)) {syntab[,i] <- apply(spec[cluster==i,]>0,2,sum)}

      samplesize <- tapply(rep(1,length(cluster)),cluster,sum)
      N = length(spec[,1])
      n = rowSums(syntab)

      phitab <-  data.frame(matrix(NA, nrow=length(names(spec)), ncol=max(cluster)))
      for (i in 1:max(cluster)) {
      for (k in 1:length(spec[1,])) {
       phitab[k,i] <- (N*syntab[k,i] - n[k]*samplesize[i]) / sqrt(n[k]*samplesize[i]*(N-n[k])*(N-samplesize[i]))
      }}

      names(phitab) <- c(seq(1, max(cluster),1))
      rownames(phitab) <- names(spec)
      result <- list("syntable" = phitab,
                     "samplesize" = samplesize)

  } else {stop("Cannot calculate synoptic table. Define correct type of species matrix values to use (abund = c('cover', 'freq')). \nCheck correct type of synoptic table output type (type = c('totalfreq', 'percfreq', 'mean', 'median', 'diffspec')).")
  }
  }
}
