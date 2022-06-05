### Implementation of Wang and Bickel 2015, by spectral clustering, for both SBM and DCSBM
## A: adjacency matrix
## Kmax: selecting from 1:Kmax
source("SP.R")
WangLRT <- function(A,Kmax,lambda=NULL){
    n <- nrow(A)
    if(is.null(lambda)) lambda <- 1/n
    d <- colSums(A)
    SBM.result <- DCSBM.result <- rep(0,Kmax)
    sbm.p <- sum(A)/(n*(n-1))
    upper.index <- which(upper.tri(A))
    a <- A[upper.index]
    sbm.ll <- sum(a*log(sbm.p)) + sum((1-a)*log(1-sbm.p))

    SBM.result[1] <- sbm.ll-lambda*n*log(n)
    SBM.clust.list <- Arash.reg.SP.path(A,K=Kmax,tau=0.25,lap=TRUE)
    for(K in 2:Kmax){
        ##  SBM
        SBM.clust <- SBM.clust.list[[K]]
        #DCSBM.clust <- Arash.reg.SSP(A,K=K,tau=0.25,lap=TRUE)
        g <- SBM.clust$cluster
        n.K <- as.numeric(table(g))
        Pi <- n.K/n
        B <- matrix(0,K,K)
        for(j in 1:(K)){
            for(k in j:K){
                j.index <- which(g==j)
                k.index <- which(g==k)
                if(j!=k){
                    B[j,k] <- mean(A[j.index,k.index])
                    B[k,j] <- B[j,k]
                }else{
                    B[j,j] <- sum(A[j.index,j.index])/(length(j.index)^2-length(j.index))
                }
            }
        }
        Z <- matrix(0,n,K)
        Z[cbind(1:n,g)] <- 1
        SBM.P <- Z%*%B%*%t(Z)
        sbm.p <- SBM.P[upper.index]
        sbm.p[sbm.p>(1-1e-6)] <- 1-1e-6
        sbm.p[sbm.p<1e-6] <- 1e-6
        sbm.ll <- sum(a*log(sbm.p)) + sum((1-a)*log(1-sbm.p)) + sum(n.K*log(Pi))
        SBM.result[K] <- sbm.ll-lambda*(K)*(K+1)*n*log(n)/2

    }
    SBM.K <- which.max(SBM.result)
    return(list(SBM.K=SBM.K,SBM.BIC=SBM.result))
}



