#' Synoptic tables and calculation of cluster-wise frequencies, fidelity and
#' differential species character
#'
#' @keywords internal
#' @noRd


syntable_wide <- function(matrix, cluster, abund = "percentage", type = "percfreq", digits = 0) {

  narep <- FALSE
  # Check for "" values and replace by NA
  if (length(matrix[matrix == ""]) != 0) {
    matrix[matrix == ""] <- NA
    narep <- TRUE
  }
  # Check for NA values
  if (any(is.na(matrix))) {
    matrix[is.na(matrix)] <- 0
    narep <- TRUE
  }
  if(narep == TRUE) print("NA and/or empty character values replaced by 0.")

  if (any(is.na(cluster))) {
    stop("NA values in cluster not allowed.")
  }

  # Check if values are numeric
  if(any(is.na(suppressWarnings(as.numeric(unlist(matrix)))))) {
    warning("Non-numeric cover values transformed into 1. Using presence/absence scale.")
    matrix[matrix != 0] <- 1
    abund = "pa"
  }

  # Conversion to numeric dataframe
  matrix <- as.data.frame(apply(matrix, 2, as.numeric))

  # Define groups
  cluster <- as.factor(cluster)
  group <- levels(cluster)


  # Create empty synoptic table
  syntab <- data.frame(matrix(NA, nrow = length(names(matrix)), ncol = length(group)))
  names(syntab) <- group
  rownames(syntab) <- names(matrix)

  # Calculate results depending on choices

  if (abund != "percentage" & abund != "pa") {
    stop("Argument 'abund' must be either percentages ('percentage') or
         presence/absence data ('pa')")
  } else {
    if (type=="totalfreq")   {
      for (i in cluster) {syntab[,i] <- apply(matrix[cluster==i,]>0, 2 ,sum)}
      samplesize <- tapply(rep(1, length(cluster)), cluster, sum)
      results <- list("syntable" = syntab,
                      "samplesize" = samplesize)

    } else if (type=="percfreq") {
      samplesize <- tapply(rep(1,length(cluster)), cluster ,sum)
      for (i in cluster) {
        syntab[,i] <- round(apply(matrix[cluster==i,]>0, 2, sum) * (100/samplesize[i]), digits=digits)}
      results <- list("syntable" = syntab, "samplesize" = samplesize)

    } else if (type=="mean" & abund=="percentage") {
      samplesize <- tapply(rep(1,length(cluster)),cluster,sum)
      for (i in cluster) { syntab[,i] <- round(apply(matrix[cluster==i,],2,mean), digits=digits) }
      results <- list("syntable" = syntab, "samplesize" = samplesize)

    } else if (type=="mean" & abund=="pa") {
      stop("Cannot calculate mean cover in clusters with presence/absence values.")

    } else if (type=="median" & abund=="percentage") {
      samplesize <- tapply(rep(1,length(cluster)),cluster,sum)
      for (i in cluster) { syntab[,i] <- round(apply(matrix[cluster==i,],2,median), digits=digits) }
      results <- list("syntable" = syntab, "samplesize" = samplesize)

    } else if (type=="median" & abund=="pa") {
      stop("Cannot calculate median cover in clusters with presence/absence values.")

    } else if (type=="diffspec") {
      # create progress bar
      pb <- txtProgressBar(min = 1, max = 23, style = 3)

      #Calculate samplesize and percentage frequency
      samplesize <- tapply(rep(1,length(cluster)),cluster,sum)
      for (i in cluster) {
        syntab[,i] <- round(apply(matrix[cluster==i,]>0, 2, sum) * (100/samplesize[i]), digits=digits)}
      syn <- syntab
      setTxtProgressBar(pb, 2)

      small <- data.frame(matrix(NA, nrow=length(rownames(syn)), ncol=length(group)))
      names(small) <- group
      for (i in 1:length(rownames(syn))) { small[i,] <- sort(as.numeric(syn[i,])) }
      setTxtProgressBar(pb, 3)

      Avsmall <- data.frame(matrix(NA, nrow=length(rownames(syn)), ncol=length(group)))
      Avsmall[,1] <- small[,1]
      for (i in 1:length(rownames(syn))) { for ( k in 2:length(group)) {
        Avsmall[i,k] <- mean(c(as.numeric(small[i,1:k]))) }   }
      setTxtProgressBar(pb, 4)

      CondAvsmall <- data.frame(matrix(NA, nrow=length(rownames(syn)), ncol=length(group)))
      CondAvsmall[,1] <- rep("", length(rownames(syn)))
      for ( k in 1:length(rownames(syn))) {
        if (small[k,2] >= 2*Avsmall[k,1]+20) {CondAvsmall[k,2] <- "1"}
        else {CondAvsmall[k,2] <- "0"}   }
      for ( i in 3:length(group)) { for ( k in 1:length(rownames(syn))) {
        if (CondAvsmall[k,i-1]=="1") {CondAvsmall[k,i] <- "1"}
        else if (small[k,i] >= 2*Avsmall[k,i-1]+20) {CondAvsmall[k,i] <- "1"}
        else {CondAvsmall[k,i] <- "0"}  } }
      setTxtProgressBar(pb, 5)

      Avsmall2 <- data.frame(matrix(NA, nrow=length(rownames(syn)), ncol=length(group)))
      Avsmall2[,1] <- Avsmall[,1]
      for ( i in 2:length(group)) { for ( k in 1:length(rownames(syn))) {
        if (CondAvsmall[k,i]=="1") {Avsmall2[k,i] <- Avsmall2[k,i-1]}
        else {Avsmall2[k,i] <- Avsmall[k,i]}  } }
      setTxtProgressBar(pb, 6)

      large <- data.frame(matrix(NA, nrow=length(rownames(syn)), ncol=length(group)))
      names(large) <- group
      for (i in 1:length(rownames(syn))) { large[i,] <- sort(as.numeric(syn[i,]), decreasing = T) }
      setTxtProgressBar(pb, 7)

      Avlarge <- data.frame(matrix(NA, nrow=length(rownames(syn)), ncol=length(group)))
      Avlarge[,1] <- large[,1]
      for (i in 1:length(rownames(syn))) { for ( k in 2:length(group)) {
        Avlarge[i,k] <- mean(c(as.numeric(large[i,1:k]))) }   }
      setTxtProgressBar(pb, 8)

      CondAvlarge <- data.frame(matrix(NA, nrow=length(rownames(syn)), ncol=length(group)))
      CondAvlarge[,1] <- rep("", length(rownames(syn)))
      for ( k in 1:length(rownames(syn))) {
        if (large[k,2] <= (Avlarge[k,1]/2)-10) {CondAvlarge[k,2] <- "1"}
        else {CondAvlarge[k,2] <- "0"}    }
      for ( i in 3:length(group)) { for ( k in 1:length(rownames(syn))) {
        if (CondAvlarge[k,i-1]=="1") {CondAvlarge[k,i] <- "1"}
        else if (large[k,i] <= (Avlarge[k,i-1]/2)-10) {CondAvlarge[k,i] <- "1"}
        else {CondAvlarge[k,i] <- "0"}   } }
      setTxtProgressBar(pb, 9)

      Avlarge2 <- data.frame(matrix(NA, nrow=length(rownames(syn)), ncol=length(group)))
      Avlarge2[,1] <- Avlarge[,1]
      for ( i in 2:length(group)) {  for ( k in 1:length(rownames(syn))) {
        if (CondAvlarge[k,i]=="1") {Avlarge2[k,i] <- Avlarge2[k,i-1]}
        else {Avlarge2[k,i] <- Avlarge[k,i]}    } }
      setTxtProgressBar(pb, 10)

      CondAvsmall2 <- data.frame(matrix(NA, nrow=length(rownames(syn)), ncol=length(group)))
      CondAvsmall2[,1] <- rep("", length(rownames(syn)))
      for ( i in 2:length(group)) {for ( k in 1:length(rownames(syn))) {
        if (small[k,i] > (Avlarge2[k,i-1]/2)-10) {CondAvsmall2[k,i] <- "1"}
        else {CondAvsmall2[k,i] <- "0"}  } }
      setTxtProgressBar(pb, 11)

      Avsmall3 <- data.frame(matrix(NA, nrow=length(rownames(syn)), ncol=length(group)))
      Avsmall3[,1] <- Avsmall2[,1]
      for ( i in 2:length(group)) { for ( k in 1:length(rownames(syn))) {
        if (CondAvsmall2[k,i]=="1") {Avsmall3[k,i] <- Avsmall3[k,i-1]}
        else {Avsmall3[k,i] <- Avsmall2[k,i]}  } }
      setTxtProgressBar(pb, 12)

      CondAvlarge2 <- data.frame(matrix(NA, nrow=length(rownames(syn)), ncol=length(group)))
      CondAvlarge2[,1] <- rep("", length(rownames(syn)))
      for ( i in 2:length(group)) { for ( k in 1:length(rownames(syn))) {
        if (large[k,i] < 2*Avsmall2[k,i-1]+20) {CondAvlarge2[k,i] <- "1"}
        else {CondAvlarge2[k,i] <- "0"} } }
      setTxtProgressBar(pb, 13)

      Avlarge3 <- data.frame(matrix(NA, nrow=length(rownames(syn)), ncol=length(group)))
      Avlarge3[,1] <- Avlarge2[,1]
      for ( i in 2:length(group)) {  for ( k in 1:length(rownames(syn))) {
        if (CondAvlarge2[k,i]=="1") {Avlarge3[k,i] <- Avlarge3[k,i-1]}
        else {Avlarge3[k,i] <- Avlarge2[k,i]}  } }
      setTxtProgressBar(pb, 14)

      Small_Avsmall3 <- data.frame(matrix(NA, nrow=length(rownames(syn)), ncol=length(group)))
      for ( i in 1:length(group)) { for ( k in 1:length(rownames(syn))) {
        if (syn[k,1] == "") {Small_Avsmall3[k,i] <- ""}
        else if (syn[k,i]>=2*max(Avsmall3[k,])+20) {Small_Avsmall3[k,i] <- "p"}
        else if (syn[k,i] < 2*max(Avsmall3[k,])+20) {Small_Avsmall3[k,i] <- "n"}
        else {stop("Calculation of differential species failed: Check synoptic table input")}  } }
      setTxtProgressBar(pb, 15)

      Large_AvLarge3 <- data.frame(matrix(NA, nrow=length(rownames(syn)), ncol=length(group)))
      for ( i in 1:length(group)) {  for ( k in 1:length(rownames(syn))) {
        if (syn[k,1] == "") {Large_AvLarge3[k,i] <- ""}
        else if (syn[k,i] >  min(Avlarge3[k,])/2-10) {Large_AvLarge3[k,i] <- "p"}
        else if (syn[k,i] <= min(Avlarge3[k,])/2-10) {Large_AvLarge3[k,i] <- "n"}
        else {stop("Calculation of differential species failed: Check synoptic table input")}   } }
      setTxtProgressBar(pb, 16)

      Sum_cond <-  data.frame(matrix(NA, nrow=length(rownames(syn)), ncol=length(group)))
      for ( i in 1:length(group)) { for ( k in 1:length(rownames(syn))) {
        Sum_cond[k,i] <- paste0(Small_Avsmall3[k,i], Large_AvLarge3[k,i])   } }
      setTxtProgressBar(pb, 17)

      Results <- data.frame(matrix(NA, nrow=length(rownames(syn)), ncol=length(group)))
      for ( i in 1:length(group)) { for ( k in 1:length(rownames(syn))) {
        if (Sum_cond[k,i] == "pp") {Results[k,i] <- "p"}
        else if (Sum_cond[k,i] == "nn") {Results[k,i] <- "n"}
        else if (Sum_cond[k,i] == "pn") {Results[k,i] <- "pn"}
        else {Results[k,i] <- "-"} } }
      rownames(Results) <- rownames(syn)
      names(Results) <- group
      setTxtProgressBar(pb, 18)

      subresult <- data.frame(matrix(NA, nrow=length(rownames(syn)), ncol=length(group)))
      for ( i in 1:length(group)) { for ( k in 1:length(rownames(Results))) {
        if (Results[k,i] == "-") {subresult[k,i] <- "0"}
        else if (Results[k,i] == "n") {subresult[k,i] <- "1"}
        else if (Results[k,i] == "p") {subresult[k,i] <- "1"}
        else if (Results[k,i] == "pn") {subresult[k,i] <- "1"}
        else {subresult[k,i] <- "0"} } }
      subresult <- as.data.frame(apply(subresult,2,as.numeric))
      setTxtProgressBar(pb, 19)

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
      setTxtProgressBar(pb, 20)

      negativespec <- data.frame(as.matrix(NA, nrow=500, ncol=200))
      for ( i in 1:length(group)) {
        if (length(rownames(Results[Results[,i]=="n",]) > 0)) {
          a <- rownames(Results[Results[,i]=="n",])
          for ( k in 1:length(a)) {
            negativespec[i,k] <- rownames(Results[Results[,i]=="n",])[k] }
        }  else {negativespec[i,1] <- NA}   }
      names(negativespec) <- rep("", length(negativespec[1,]))
      setTxtProgressBar(pb, 21)

      positivenegativespec <- data.frame(as.matrix(NA, nrow=500, ncol=200))
      for ( i in 1:length(group)) {
        if (length(rownames(Results[Results[,i]=="pn",]) > 0)) {
          a <- rownames(Results[Results[,i]=="pn",])
          for ( k in 1:length(a)) {
            positivenegativespec[i,k] <- rownames(Results[Results[,i]=="pn",])[k] }
        } else {positivenegativespec[i,1] <- NA}   }
      names(positivenegativespec) <- rep("", length(positivenegativespec[1,]))
      setTxtProgressBar(pb, 22)

      diffspeclist <- list()
      for (i in 1:length(group)) {
        pos <- positivespec[i, ]
        pos <- pos[!is.na(pos)]
        if (length(pos) > 0) {pos=pos}  else {pos <- noquote("no positive diff species")}
        neg <- negativespec[i,]
        neg <- neg[!is.na(neg)]
        if (length(neg) > 0) {neg=neg}  else {neg <- noquote("no negative diff species")}
        posneg <- positivenegativespec[i,]
        posneg <- posneg[!is.na(posneg)]
        if (length(posneg) > 0) {posneg=posneg} else {posneg <- noquote("no positive/negative diff species")}
        name <- group[i]
        tmp <- list("positive diff" = pos, "negative diff" = neg, "positive/negative diff" = posneg)
        diffspeclist[[name]] <- tmp }
      setTxtProgressBar(pb, 23)

      # Close progress bar
      close(pb)

      results <- list("syntable" = Results,
                      "onlydiff" = onlydiff,
                      "others" = others,
                      "samplesize" = samplesize,
                      "differentials" = diffspeclist)

    } else if (type == "phi") {
      # Calculate total frequency and sample size
      for (i in cluster) {syntab[,i] <- apply(matrix[cluster == i,]>0, 2 ,sum)}

      samplesize <- tapply(rep(1, length(cluster)), cluster, sum)
      N = length(matrix[,1])
      n = rowSums(syntab)

      phitab <-  data.frame(matrix(NA, nrow=length(names(matrix)), ncol=length(group)))
      for (i in 1:length(group)) {
        for (k in 1:length(matrix[1,])) {
          phitab[k,i] <- (N * syntab[k,i] - n[k] * samplesize[i]) / sqrt(n[k] * samplesize[i]*(N-n[k]) * (N-samplesize[i]))
        }}

      names(phitab) <- group
      rownames(phitab) <- names(matrix)
      results <- list("syntable" = phitab,
                      "samplesize" = samplesize)

    } else {stop("Cannot calculate synoptic table. Define correct type of species matrix values to use (abund = c('percentage', 'pa')). \nCheck correct type of synoptic table output type (type = c('totalfreq', 'percfreq', 'mean', 'median', 'diffspec')).")
    }
  }

  return(invisible(results))

}
