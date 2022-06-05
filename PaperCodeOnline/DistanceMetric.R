Scaled.Frob <- function(D1,D2){
    scaled.diff <- (D1-D2)*log(1+2*D1)
    return(sqrt(sum(scaled.diff^2)/2))
}


Normalized.Scaled.Frob <- function(D1,D2){
    d.max <- max(D1)
    scaled.diff <- (D1-D2)*log(1+2*D1)/log(1+d.max)
    return(sqrt(sum(scaled.diff^2)/2))
}


Scaled.FP <- function(G1,G2,D1){
    upper.index <- which(upper.tri(G1))
    g1 <- G1[upper.index]
    g2 <- G2[upper.index]
    d1 <- D1[upper.index]
    same.clust.index <- which(g2==1)
    scaled.FP <- sum((1-g1[same.clust.index])*log(2*d1[same.clust.index]+1))
    scaled.FP <- scaled.FP/length(same.clust.index)
    return(scaled.FP)
}

Standard.FP <- function(G1,G2,D1){
    upper.index <- which(upper.tri(G1))
    g1 <- G1[upper.index]
    g2 <- G2[upper.index]
    d1 <- D1[upper.index]
    same.clust.index <- which(g2==1)
    scaled.FP <- sum((1-g1[same.clust.index]))
    scaled.FP <- scaled.FP/length(same.clust.index)
    return(scaled.FP)
}



Scaled.FN <- function(G1,G2,D2){
    upper.index <- which(upper.tri(G1))
    g1 <- G1[upper.index]
    g2 <- G2[upper.index]
    d2 <- D2[upper.index]
    diff.clust.index <- which(g2==0)
    scaled.FN <- sum((g1[diff.clust.index])*exp(2*d2[diff.clust.index]+1))
    scaled.FN <- scaled.FN/length(diff.clust.index)
    return(scaled.FN)
}



Binary.Similarity <- function(s1,s2){
    n <- min(length(s1),length(s2))
    min(which(s1[1:n]!=s2[1:n]))
}
