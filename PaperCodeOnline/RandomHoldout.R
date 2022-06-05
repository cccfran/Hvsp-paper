#library(pROC)
library(AUC)
#library(softImpute)
library(Matrix)
source("SimonFunk.R")
#### Randomly holdout 1/K of the edges to evaluate number of communities
EdgeCV <- function(A,max.K,cv=NULL,B=3,holdout.p=0.1,soft=TRUE,tau=0,dc.est=2,fast=FALSE,kappa=NULL){
    n <- nrow(A)
    edge.index <- which(upper.tri(A))
    edge.n <- length(edge.index)
    holdout.index.list <- list()
    if(is.null(cv)){
        holdout.n <- floor(holdout.p*edge.n)

        for(j in 1:B){
            holdout.index.list[[j]] <- sample(x=edge.n,size=holdout.n)
        }
    }else{
        sample.index <- sample.int(edge.n)
        max.fold.num <- ceiling(edge.n/cv)
        fold.index <- rep(1:cv,each=max.fold.num)[edge.n]
        cv.index <- fold.index[sample.index]
        B <- cv
        for(j in 1:B){
            holdout.index.list[[j]] <- which(cv.index==j)
        }
    }
    #print(fast)
    result <- lapply(holdout.index.list,holdout.evaluation,A=A,max.K=max.K,soft=soft,tau=tau,dc.est=dc.est,fast=fast,p.sample=1-holdout.p,kappa=kappa)
    dc.block.err.mat <- dc.loglike.mat <- bin.dev.mat <- roc.auc.mat <- impute.err.mat <- block.err.mat <- loglike.mat <- matrix(0,nrow=B,ncol=max.K)
    no.edge.seq <- rep(0,B)
    Omega.list <- A.list <- Imputed.A.list <- list()
    for(b in 1:B){
        impute.err.mat[b,] <- result[[b]]$impute.sq.err
        block.err.mat[b,] <- result[[b]]$block.sq.err
        loglike.mat[b,] <- result[[b]]$loglike
        roc.auc.mat[b,] <- result[[b]]$roc.auc
        bin.dev.mat[b,] <- result[[b]]$bin.dev
        no.edge.seq[b] <- result[[b]]$no.edge
        dc.block.err.mat[b,] <- result[[b]]$dc.block.sq.err
        dc.loglike.mat[b,] <- result[[b]]$dc.loglike

    }

    return(list(impute.err=colMeans(impute.err.mat),block.err=colMeans(block.err.mat),loglike=colSums(loglike.mat),auc=colMeans(roc.auc.mat),no.edge.seq=no.edge.seq,dc.block.err=colMeans(dc.block.err.mat),dc.loglike=colSums(dc.loglike.mat),bin= colMeans(bin.dev.mat),auc.mat=roc.auc.mat,dev.mat=loglike.mat,l2.mat=block.err.mat,SSE.mat=impute.err.mat,dc.dev.mat=dc.loglike.mat,dc.l2.mat=dc.block.err.mat))
}



holdout.evaluation <- function(holdout.index,A,max.K,soft=TRUE,tau=0,dc.est=1,fast=FALSE,p.sample=1,kappa=NULL){
    n <- nrow(A)
    edge.index <- which(upper.tri(A))
    edge.n <- length(edge.index)
    A.new <- matrix(0,n,n)
    A.new[upper.tri(A.new)] <- A[edge.index]
    A.new[edge.index[holdout.index]] <- NA
    A.new <- A.new + t(A.new)
    degrees <- colSums(A.new,na.rm=TRUE)
    no.edge <- 0
    no.edge <- sum(degrees==0)

    Omega <- which(is.na(A.new))
    non.miss <- which(!is.na(A.new))
    #A.new[non.miss] <- A.new[non.miss] + 0.5
    dc.block.sq.err <-  dc.loglike <- roc.auc <- bin.dev <- block.sq.err <- impute.sq.err <- loglike <- rep(0,max.K)

    for(k in 1:max.K){
        print(k)
        #print(fast)
        if(k==1){
            if(soft){
                SVD <- softImpute(x=as.matrix(A.new),rank.max=k,maxit=300,type="svd")
                tmp.est <- list(A=matrix(SVD$u,ncol=k)%*%t(matrix(SVD$v,ncol=k))*as.numeric(SVD$d),SVD=SVD)

            }else{
                tmp.est <- iter.SVD.core(A.new,K=k,fast=fast,p.sample=p.sample,kappa=kappa)#SF.SVD(A.new,K=k,verbose=FALSE,nstart=5)

            }
        }else{
            if(soft){
                SVD <- softImpute(x=as.matrix(A.new),rank.max=k,warm.start=tmp.est$SVD,maxit=300)
                tmp.est <- list(A=SVD$u%*%diag(SVD$d)%*%t(SVD$v),SVD=SVD)

            }else{
                tmp.est <- iter.SVD.core(A.new,K=k,init=tmp.est$SVD,fast=fast,p.sample=p.sample,kappa=kappa)#SF.SVD(A.new,K=k,verbose=FALSE,nstart=5)
            }
        }
        A.approx <- tmp.est$A
        impute.sq.err[k] <- sum((A.approx[Omega]-A[Omega])^2)
        response <- A[Omega]
        predictors <- A.approx[Omega]
        print("AUC claculation")
        #print(system.time(tmp.roc <- pROC::roc(response=response,predictor=predictors)))
        #print(length(unique(predictors)))
        aa <- AUC::roc(predictions=predictors,labels=factor(response))
        #tmp.roc.smooth <- smooth(tmp.roc,method="binormal")
        roc.auc[k] <- auc(aa)#as.numeric(tmp.roc$auc)
        #print(tmp.roc$auc)
        #print(auc(aa))
        #roc.auc[k] <- as.numeric(tmp.roc.smooth$auc)
        trunc.predictors <- predictors
        trunc.predictors[predictors>(1-1e-6)] <- 1-1e-6
        trunc.predictors[predictors<1e-6] <- 1e-6
        bin.dev[k] <- sum((response-trunc.predictors)^2)#-sum(response*log(trunc.predictors)) - sum((1-response)*log(1-trunc.predictors))
        if(k==1){
            pb <- (sum(A.new,na.rm=TRUE)+1)/(sum(!is.na(A.new)) -sum(!is.na(diag(A.new)))+1)
            if(pb < 1e-6) pb <- 1e-6
            if(pb > 1-1e-6) pb <- 1-1e-6
            A.Omega <- A[Omega]
            block.sq.err[k] <- sum((pb-A[Omega])^2)
            loglike[k] <- -sum(A.Omega*log(pb)) - sum((1-A.Omega)*log(1-pb))

        }

        #U.approx <- eigen(A.approx)$vectors[,1:k]
        print("SBM calculation")
        ptm <- proc.time()
        if(k==1) {U.approx <- matrix(tmp.est$SVD$v,ncol=k)}else{
            U.approx <- tmp.est$SVD$v[,1:k]
            if(tau>0){
            A.approx <- A.approx + tau*mean(colSums(A.approx))/n
            d.approx <- colSums(A.approx)
            L.approx <- diag(1/sqrt(d.approx))%*%A.approx%*%diag(1/sqrt(d.approx))
            A.approx.svd <- irlba(L.approx,nu=k,nv=k)
            U.approx <- A.approx.svd$v[,1:k]
            }
        }
        km <- kmeans(U.approx,centers=k,nstart=30,iter.max=30)
        B <- matrix(0,k,k)
        Theta <- matrix(0,n,k)
        for(i in 1:k){
            for(j in i:k){
                N.i <- which(km$cluster==i)
                N.j <- which(km$cluster==j)
                if(i!=j){
                    B[i,j] <- B[j,i] <- (sum(A.new[N.i,N.j],na.rm=TRUE)+1)/(sum(!is.na(A.new[N.i,N.j]))+1)
                } else{
                   B[i,j] <- B[j,i] <- (sum(A.new[N.i,N.j],na.rm=TRUE)+1)/(sum(!is.na(A.new[N.i,N.j])) -sum(!is.na(diag(A.new[N.i,N.j])))+1)
                }

            }
            Theta[N.i,i] <- 1
        }
        P.hat <- Theta%*%B%*%t(Theta)
        diag(P.hat) <- 0
        block.sq.err[k] <- sum((P.hat[Omega]-A[Omega])^2)
        P.hat.Omega <- P.hat[Omega]
        A.Omega <- A[Omega]
        P.hat.Omega[P.hat.Omega < 1e-6] <- 1e-6
        P.hat.Omega[P.hat.Omega > (1-1e-6)] <- 1-1e-6
        loglike[k] <- -sum(A.Omega*log(P.hat.Omega)) - sum((1-A.Omega)*log(1-P.hat.Omega))
        print(proc.time() - ptm)
#### Degree correct model
        V <- U.approx
        print("DCSBM calculation")
        ptm <- proc.time()
        #V.norms <- apply(V,1,function(x) sqrt(sum(x^2)))
        if(k==1) {V.norms <- as.numeric(abs(V))}else{
            V.norms <- apply(V,1,function(x) sqrt(sum(x^2)))
        }

        iso.index <- which(V.norms==0)
        Psi <- V.norms
        Psi <- Psi / max(V.norms)
        inv.V.norms <- 1/V.norms
        inv.V.norms[iso.index] <- 1

        V.normalized <- diag(as.numeric(inv.V.norms))%*%V

        #V.norms <- apply(V,1,function(x) sqrt(sum(x^2)))
        #Psi <- V.norms
        #Psi <- Psi / max(V.norms)
        #V.normalized <- diag(1/V.norms)%*%V
        #Psi.outer <- outer(Psi,Psi)
        if(k==1){
        if(dc.est>1){
            B <- sum(A.new,na.rm=TRUE)+0.01

            partial.d <- colSums(A.new,na.rm=TRUE)
            partial.gd <- B
            phi <- rep(0,n)
            B.g <- partial.gd
            phi <- as.numeric(partial.d/B.g)
            B <- B/p.sample
            P.hat <- t(t(matrix(B,n,n)*phi)*phi)
            #P.hat <- diag(phi)%*%matrix(B,n,n)%*%diag(phi)
            diag(P.hat) <- 0
        }
            dc.block.sq.err[k] <- sum((pb-A[Omega])^2)
            P.hat.Omega <- P.hat[Omega]
            A.Omega <- A[Omega]
            P.hat.Omega[P.hat.Omega < 1e-6] <- 1e-6
            P.hat.Omega[P.hat.Omega > (1-1e-6)] <- 1-1e-6

            dc.loglike[k] <- -sum(A.Omega*log(P.hat.Omega)) - sum((1-A.Omega)*log(1-P.hat.Omega))


        }else{
        km <- kmeans(V.normalized,centers=k,nstart=30,iter.max=30)
        if(dc.est>1){
            B <- matrix(0,k,k)
            Theta <- matrix(0,n,k)
            for(i in 1:k){
                for(j in 1:k){
                N.i <- which(km$cluster==i)
                N.j <- which(km$cluster==j)
                B[i,j] <- sum(A.new[N.i,N.j],na.rm=TRUE)+0.01
                }
                Theta[N.i,i] <- 1
            }
            Theta <- Matrix(Theta,sparse=TRUE)
            partial.d <- colSums(A.new,na.rm=TRUE)
            partial.gd <- colSums(B)
            phi <- rep(0,n)
            B.g <- Theta%*%partial.gd
            phi <- as.numeric(partial.d/B.g)
            B <- B/p.sample
            tmp.int.mat <- Theta*phi
            P.hat <-as.matrix(tmp.int.mat%*%B%*%t(tmp.int.mat))
            #P.hat <- diag(phi)%*%Theta%*%B%*%t(Theta)%*%diag(phi)
            diag(P.hat) <- 0
        }
        dc.block.sq.err[k] <- sum((P.hat[Omega]-A[Omega])^2)
        P.hat.Omega <- P.hat[Omega]
        A.Omega <- A[Omega]
        P.hat.Omega[P.hat.Omega < 1e-6] <- 1e-6
        P.hat.Omega[P.hat.Omega > (1-1e-6)] <- 1-1e-6
        dc.loglike[k] <- -sum(A.Omega*log(P.hat.Omega)) - sum((1-A.Omega)*log(1-P.hat.Omega))
    }
        print(proc.time() - ptm)



    }
    return(list(impute.sq.err=impute.sq.err,block.sq.err=block.sq.err,loglike=loglike,roc.auc=roc.auc,no.edge=no.edge,dc.block.sq.err=dc.block.sq.err,dc.loglike=dc.loglike,bin.dev=bin.dev))
}



#### Using AUC to evaluate intrinsic rank
AUC.Rank <- function(A,max.K,B=3,holdout.p=0.1,soft=FALSE,fast=TRUE){
    n <- nrow(A)
    edge.index <- which(upper.tri(A))
    edge.n <- length(edge.index)
    holdout.index.list <- list()

    holdout.n <- floor(holdout.p*edge.n)

    for(j in 1:B){
        holdout.index.list[[j]] <- sample(x=edge.n,size=holdout.n)
    }
    result <- lapply(holdout.index.list,missing.AUC,A=A,max.K=max.K,fast=fast,sample.p=1-holdout.p)
    roc.auc.mat <- matrix(0,nrow=B,ncol=max.K)

    for(b in 1:B){
        roc.auc.mat[b,] <- result[[b]]$roc.auc
    }

    auc.seq <- colMeans(roc.auc.mat)
    auc.sd <- apply(roc.auc.mat,2,sd)/sqrt(B)
    return(list(rank=which.max(auc.seq),auc=auc.seq,auc.sd=auc.sd,rank.1sd=min(which(auc.seq>(max(auc.seq)-auc.sd[which.max(auc.seq)])))))
    #return(list(rank=which.max(auc.seq),auc=auc.seq))
}



missing.AUC <- function(holdout.index,A,max.K,soft=FALSE,fast=fast,sample.p=sample.p){
    n <- nrow(A)
    edge.index <- which(upper.tri(A))
    edge.n <- length(edge.index)
    A.new <- matrix(0,n,n)
    A.new[upper.tri(A.new)] <- A[edge.index]
    A.new[edge.index[holdout.index]] <- NA
    A.new <- A.new + t(A.new)
    degrees <- colSums(A.new,na.rm=TRUE)
    no.edge <- 0
    no.edge <- sum(degrees==0)

    Omega <- which(is.na(A.new))
    imputed.A <- list()
    roc.auc <- rep(0,max.K)

    for(k in max.K:1){
        print(k)
        if(k==max.K){
            if(soft){
                SVD <- softImpute(x=as.matrix(A.new),rank.max=k)
                tmp.est <- list(A=SVD$u%*%diag(SVD$d)%*%t(SVD$v),SVD=SVD)

            }else{
                tmp.est <- iter.SVD.core(A.new,K=k)#SF.SVD(A.new,K=k,verbose=FALSE,nstart=5)
            }
        }else{
            if(soft){
                SVD <- softImpute(x=as.matrix(A.new),rank.max=k,warm.start=tmp.est$SVD)
                tmp.est <- list(A=SVD$u%*%diag(SVD$d)%*%t(SVD$v),SVD=SVD)

            }else{
                tmp.est <- iter.SVD.core(A.new,K=k,init=tmp.est$SVD)#SF.SVD(A.new,K=k,verbose=FALSE,nstart=5)
            }
        }
        A.approx <- tmp.est$A
        response <- A[Omega]
        predictors <- A.approx[Omega]
        tmp.roc <- roc(response=response,predictor=predictors)
        roc.auc[k] <- as.numeric(tmp.roc$auc)
        imputed.A[[k]] <- A.approx
    }
    return(list(roc.auc=roc.auc,imputed.A=imputed.A,Omega=Omega))
}







#### Using AUC to evaluate intrinsic rank
ECV.directed.Rank <- function(A,max.K,B=3,holdout.p=0.1,soft=FALSE){
    n <- nrow(A)
    edge.index <- 1:n^2
    edge.n <- length(edge.index)
    #edge.index <- which(upper.tri(A))
    #edge.n <- length(edge.index)

    holdout.index.list <- list()

    holdout.n <- floor(holdout.p*edge.n)

    for(j in 1:B){
        holdout.index.list[[j]] <- sample(x=edge.n,size=holdout.n)
    }
    result <- lapply(holdout.index.list,missing.Rank.fast.all,A=A,max.K=max.K,p.sample=1-holdout.p)
    sse.mat <- roc.auc.mat <- matrix(0,nrow=B,ncol=max.K)

    for(b in 1:B){
        roc.auc.mat[b,] <- result[[b]]$roc.auc
        sse.mat[b,] <- result[[b]]$sse
    }

    auc.seq <- colMeans(roc.auc.mat)
    auc.sd <- apply(roc.auc.mat,2,sd)/sqrt(B)
    sse.seq <- colMeans(sse.mat)
    sse.sd <- apply(sse.mat,2,sd)/sqrt(B)
    return(list(rank=which.max(auc.seq),auc=auc.seq,auc.sd=auc.sd,sse=sse.seq,sse.sd=sse.sd))
    #return(list(rank=which.max(auc.seq),auc=auc.seq))
}



missing.Rank <- function(holdout.index,A,max.K,soft=FALSE){
    n <- nrow(A)
    A.new <- A
    A.new[holdout.index] <- NA
    degrees <- colSums(A.new,na.rm=TRUE)
    no.edge <- 0
    no.edge <- sum(degrees==0)

    Omega <- which(is.na(A.new))
    imputed.A <- list()
    sse <- roc.auc <- rep(0,max.K)

    for(k in 1:max.K){
        print(k)
        if(k==1){
            if(soft){
                SVD <- softImpute(x=as.matrix(A.new),rank.max=k)
                tmp.est <- list(A=SVD$u%*%diag(SVD$d)%*%t(SVD$v),SVD=SVD)

            }else{
                tmp.est <- iter.SVD.core(A.new,K=k)#SF.SVD(A.new,K=k,verbose=FALSE,nstart=5)
            }
        }else{
            if(soft){
                SVD <- softImpute(x=as.matrix(A.new),rank.max=k,warm.start=tmp.est$SVD)
                tmp.est <- list(A=SVD$u%*%diag(SVD$d)%*%t(SVD$v),SVD=SVD)

            }else{
                tmp.est <- iter.SVD.core(A.new,K=k,init=tmp.est$SVD)#SF.SVD(A.new,K=k,verbose=FALSE,nstart=5)
            }
        }
        A.approx <- tmp.est$A
        response <- A[Omega]
        predictors <- A.approx[Omega]
        aa <- AUC::roc(predictions=predictors,labels=factor(response))
        roc.auc[k] <- auc(aa)
        sse[k] <- mean((response-predictors)^2)
        imputed.A[[k]] <- A.approx
    }
    return(list(roc.auc=roc.auc,imputed.A=imputed.A,Omega=Omega,sse=sse))
}






missing.Rank.fast.all <- function(holdout.index,A,max.K,p.sample=NULL){
    n <- nrow(A)
    A.new <- A
    A.new[holdout.index] <- NA
    if(is.null(p.sample)){
        p.sample <- 1-length(holdout.index)/n^2
    }
    degrees <- colSums(A.new,na.rm=TRUE)
    no.edge <- 0
    no.edge <- sum(degrees==0)

    Omega <- which(is.na(A.new))
    imputed.A <- list()
    sse <- roc.auc <- rep(0,max.K)
    SVD.result <- iter.SVD.core.fast.all(A.new,max.K,fast=TRUE,p.sample=p.sample)

    for(k in 1:max.K){
        print(k)
        tmp.est <- SVD.result[[k]]
        A.approx <- tmp.est$A.thr
        response <- A[Omega]
        predictors <- A.approx[Omega]
        aa <- AUC::roc(predictions=predictors,labels=factor(response))
        roc.auc[k] <- auc(aa)
        sse[k] <- mean((response-predictors)^2)
        imputed.A[[k]] <- A.approx
    }
    return(list(roc.auc=roc.auc,imputed.A=imputed.A,Omega=Omega,sse=sse))
}








#### Using AUC to evaluate intrinsic rank
ECV.undirected.Rank.weighted <- function(A,max.K,B=3,holdout.p=0.1,soft=FALSE,fast=fast){
    n <- nrow(A)
    #edge.index <- 1:n^2
    #edge.n <- length(edge.index)
    edge.index <- which(upper.tri(A))
    edge.n <- length(edge.index)

    holdout.index.list <- list()

    holdout.n <- floor(holdout.p*edge.n)

    for(j in 1:B){
        holdout.index.list[[j]] <- sample(x=edge.n,size=holdout.n)
    }
    result <- lapply(holdout.index.list,missing.undirected.Rank.weighted.fast.all,A=A,max.K=max.K,soft=soft,fast=fast,p.sample=1-holdout.p)
    sse.mat <- roc.auc.mat <- matrix(0,nrow=B,ncol=max.K)

    for(b in 1:B){
        #roc.auc.mat[b,] <- result[[b]]$roc.auc
        sse.mat[b,] <- result[[b]]$sse
    }

    #auc.seq <- colMeans(roc.auc.mat)
    #auc.sd <- apply(roc.auc.mat,2,sd)/sqrt(B)
    sse.seq <- colMeans(sse.mat)
    sse.sd <- apply(sse.mat,2,sd)/sqrt(B)
    return(list(sse=sse.seq,sse.sd=sse.sd))
    #return(list(rank=which.max(auc.seq),auc=auc.seq))
}



missing.undirected.Rank.weighted <- function(holdout.index,A,max.K,soft=FALSE,fast=fast,p.sample=1){
    n <- nrow(A)
    #A.new <- A
    #A.new[holdout.index] <- NA
    edge.index <- which(upper.tri(A))
    edge.n <- length(edge.index)
    A.new <- matrix(0,n,n)
    A.new[upper.tri(A.new)] <- A[edge.index]
    A.new[edge.index[holdout.index]] <- NA
    A.new <- A.new + t(A.new)
    diag(A.new) <- diag(A)
    degrees <- colSums(A.new,na.rm=TRUE)
    no.edge <- 0
    no.edge <- sum(degrees==0)

    Omega <- which(is.na(A.new))
    imputed.A <- list()
    sse <- roc.auc <- rep(0,max.K)

    for(k in 1:max.K){
        print(k)
        if(k==1){
            if(soft){
                SVD <- softImpute(x=as.matrix(A.new),rank.max=k)
                tmp.est <- list(A=SVD$u%*%diag(SVD$d)%*%t(SVD$v),SVD=SVD)

            }else{
                tmp.est <- iter.SVD.core(A.new,K=k,fast=fast,p.sample=p.sample)#SF.SVD(A.new,K=k,verbose=FALSE,nstart=5)
            }
        }else{
            if(soft){
                SVD <- softImpute(x=as.matrix(A.new),rank.max=k,warm.start=tmp.est$SVD)
                tmp.est <- list(A=SVD$u%*%diag(SVD$d)%*%t(SVD$v),SVD=SVD)

            }else{
                tmp.est <- iter.SVD.core(A.new,K=k,init=tmp.est$SVD,fast=fast,p.sample=p.sample)#SF.SVD(A.new,K=k,verbose=FALSE,nstart=5)
            }
        }
        if(k==1){
        A.approx <- matrix(tmp.est$SVD$u,ncol=1)%*%t(matrix(tmp.est$SVD$v,ncol=1))*tmp.est$SVD$d[1]
        }else{
            A.approx <- tmp.est$SVD$u%*%t(tmp.est$SVD$v*tmp.est$SVD$d)
        }
        response <- A[Omega]
        predictors <- A.approx[Omega]
        #aa <- AUC::roc(predictions=predictors,labels=factor(response))
        #roc.auc[k] <- auc(aa)
        sse[k] <- mean((response-predictors)^2)
        imputed.A[[k]] <- A.approx
    }
    return(list(imputed.A=imputed.A,Omega=Omega,sse=sse))
}









missing.undirected.Rank.weighted.fast.all <- function(holdout.index,A,max.K,soft=FALSE,fast=fast,p.sample=1){
    n <- nrow(A)
    #A.new <- A
    #A.new[holdout.index] <- NA
    edge.index <- which(upper.tri(A))
    edge.n <- length(edge.index)
    A.new <- matrix(0,n,n)
    A.new[upper.tri(A.new)] <- A[edge.index]
    A.new[edge.index[holdout.index]] <- NA
    A.new <- A.new + t(A.new)
    diag(A.new) <- diag(A)
    degrees <- colSums(A.new,na.rm=TRUE)
    no.edge <- 0
    no.edge <- sum(degrees==0)

    Omega <- which(is.na(A.new))
    imputed.A <- list()
    sse <- roc.auc <- rep(0,max.K)
    SVD.result <- iter.SVD.core.fast.all(A.new,max.K,fast=TRUE,p.sample=p.sample)
    for(k in 1:max.K){
        print(k)
        tmp.est <- SVD.result[[k]]
        #if(k==1){
        #A.approx <- matrix(tmp.est$SVD$u,ncol=1)%*%t(matrix(tmp.est$SVD$v,ncol=1))*tmp.est$SVD$d[1]
        #}else{
         #   A.approx <- tmp.est$SVD$u%*%t(tmp.est$SVD$v*tmp.est$SVD$d)
        #}
        A.approx <- tmp.est$A
        response <- A[Omega]
        predictors <- A.approx[Omega]
        #aa <- AUC::roc(predictions=predictors,labels=factor(response))
        #roc.auc[k] <- auc(aa)
        sse[k] <- mean((response-predictors)^2)
        imputed.A[[k]] <- A.approx
    }
    return(list(imputed.A=imputed.A,Omega=Omega,sse=sse))
}










EdgeCV.fast.all <- function(A,max.K,cv=NULL,B=3,holdout.p=0.1,soft=TRUE,tau=0,dc.est=2,fast=FALSE,kappa=NULL){
    n <- nrow(A)
    edge.index <- which(upper.tri(A))
    edge.n <- length(edge.index)
    holdout.index.list <- list()
    if(is.null(cv)){
        holdout.n <- floor(holdout.p*edge.n)

        for(j in 1:B){
            holdout.index.list[[j]] <- sample(x=edge.n,size=holdout.n)
        }
    }else{
        sample.index <- sample.int(edge.n)
        max.fold.num <- ceiling(edge.n/cv)
        fold.index <- rep(1:cv,each=max.fold.num)[edge.n]
        cv.index <- fold.index[sample.index]
        B <- cv
        for(j in 1:B){
            holdout.index.list[[j]] <- which(cv.index==j)
        }
    }
    #print(fast)
    result <- lapply(holdout.index.list,holdout.evaluation.fast.all,A=A,max.K=max.K,soft=soft,tau=tau,dc.est=dc.est,fast=fast,p.sample=1-holdout.p,kappa=kappa)
    dc.block.err.mat <- dc.loglike.mat <- bin.dev.mat <- roc.auc.mat <- impute.err.mat <- block.err.mat <- loglike.mat <- matrix(0,nrow=B,ncol=max.K)
    no.edge.seq <- rep(0,B)
    Omega.list <- A.list <- Imputed.A.list <- list()
    for(b in 1:B){
        impute.err.mat[b,] <- result[[b]]$impute.sq.err
        block.err.mat[b,] <- result[[b]]$block.sq.err
        loglike.mat[b,] <- result[[b]]$loglike
        roc.auc.mat[b,] <- result[[b]]$roc.auc
        bin.dev.mat[b,] <- result[[b]]$bin.dev
        no.edge.seq[b] <- result[[b]]$no.edge
        dc.block.err.mat[b,] <- result[[b]]$dc.block.sq.err
        dc.loglike.mat[b,] <- result[[b]]$dc.loglike

    }

    return(list(impute.err=colMeans(impute.err.mat),block.err=colMeans(block.err.mat),loglike=colSums(loglike.mat),auc=colMeans(roc.auc.mat),no.edge.seq=no.edge.seq,dc.block.err=colMeans(dc.block.err.mat),dc.loglike=colSums(dc.loglike.mat),bin= colMeans(bin.dev.mat),auc.mat=roc.auc.mat,dev.mat=loglike.mat,l2.mat=block.err.mat,SSE.mat=impute.err.mat,dc.dev.mat=dc.loglike.mat,dc.l2.mat=dc.block.err.mat))
}



holdout.evaluation.fast.all <- function(holdout.index,A,max.K,soft=TRUE,tau=0,dc.est=1,fast=FALSE,p.sample=1,kappa=NULL){
    n <- nrow(A)
    edge.index <- which(upper.tri(A))
    edge.n <- length(edge.index)
    A.new <- matrix(0,n,n)
    A.new[upper.tri(A.new)] <- A[edge.index]
    A.new[edge.index[holdout.index]] <- NA
    A.new <- A.new + t(A.new)
    degrees <- colSums(A.new,na.rm=TRUE)
    no.edge <- 0
    no.edge <- sum(degrees==0)

    Omega <- which(is.na(A.new))
    non.miss <- which(!is.na(A.new))
    #A.new[non.miss] <- A.new[non.miss] + 0.5
    SVD.result <- iter.SVD.core.fast.all(A.new,max.K,fast=TRUE,p.sample=p.sample)
    dc.block.sq.err <-  dc.loglike <- roc.auc <- bin.dev <- block.sq.err <- impute.sq.err <- loglike <- rep(0,max.K)

    for(k in 1:max.K){
        print(k)
        #print(fast)
        tmp.est <- SVD.result[[k]]
        A.approx <- tmp.est$A.thr
        impute.sq.err[k] <- sum((A.approx[Omega]-A[Omega])^2)
        response <- A[edge.index[holdout.index]]#A[Omega]
        predictors <- A.approx[edge.index[holdout.index]]#A.approx[Omega]
        print("AUC claculation")
        #print(system.time(tmp.roc <- pROC::roc(response=response,predictor=predictors)))
        #print(length(unique(predictors)))
        aa <- AUC::roc(predictions=predictors,labels=factor(response))
        #tmp.roc.smooth <- smooth(tmp.roc,method="binormal")
        roc.auc[k] <- auc(aa)#as.numeric(tmp.roc$auc)
        #print(tmp.roc$auc)
        #print(auc(aa))
        #roc.auc[k] <- as.numeric(tmp.roc.smooth$auc)
        trunc.predictors <- predictors
        trunc.predictors[predictors>(1-1e-6)] <- 1-1e-6
        trunc.predictors[predictors<1e-6] <- 1e-6
        bin.dev[k] <- sum((response-trunc.predictors)^2)#-sum(response*log(trunc.predictors)) - sum((1-response)*log(1-trunc.predictors))
        if(k==1){
            pb <- (sum(A.new,na.rm=TRUE)+1)/(sum(!is.na(A.new)) -sum(!is.na(diag(A.new)))+1)
            if(pb < 1e-6) pb <- 1e-6
            if(pb > 1-1e-6) pb <- 1-1e-6
            A.Omega <- A[Omega]
            block.sq.err[k] <- sum((pb-A[Omega])^2)
            loglike[k] <- -sum(A.Omega*log(pb)) - sum((1-A.Omega)*log(1-pb))

        }

        #U.approx <- eigen(A.approx)$vectors[,1:k]
        print("SBM calculation")
        #print(k)
        #print(dim(tmp.est$SVD$v))
        ptm <- proc.time()
        if(k==1) {U.approx <- matrix(tmp.est$SVD$v,ncol=k)}else{
            U.approx <- tmp.est$SVD$v[,1:k]
            if(tau>0){
            A.approx <- A.approx + tau*mean(colSums(A.approx))/n
            d.approx <- colSums(A.approx)
            L.approx <- diag(1/sqrt(d.approx))%*%A.approx%*%diag(1/sqrt(d.approx))
            A.approx.svd <- irlba(L.approx,nu=k,nv=k)
            U.approx <- A.approx.svd$v[,1:k]
            }
        }

        km <- kmeans(U.approx,centers=k,nstart=30,iter.max=30)
        B <- matrix(0,k,k)
        Theta <- matrix(0,n,k)
        for(i in 1:k){
            for(j in i:k){
                N.i <- which(km$cluster==i)
                N.j <- which(km$cluster==j)
                if(i!=j){
                    B[i,j] <- B[j,i] <- (sum(A.new[N.i,N.j],na.rm=TRUE)+1)/(sum(!is.na(A.new[N.i,N.j]))+1)
                } else{
                    #print(max(N.i))
                    #print(max(N.j))
                    #print(dim(A.new))
                   B[i,j] <- B[j,i] <- (sum(A.new[N.i,N.j],na.rm=TRUE)+1)/(sum(!is.na(A.new[N.i,N.j])) -sum(!is.na(diag(A.new[N.i,N.j])))+1)
                }

            }
            Theta[N.i,i] <- 1
        }
        P.hat <- Theta%*%B%*%t(Theta)
        diag(P.hat) <- 0
        block.sq.err[k] <- sum((P.hat[Omega]-A[Omega])^2)
        P.hat.Omega <- P.hat[Omega]
        A.Omega <- A[Omega]
        P.hat.Omega[P.hat.Omega < 1e-6] <- 1e-6
        P.hat.Omega[P.hat.Omega > (1-1e-6)] <- 1-1e-6
        loglike[k] <- -sum(A.Omega*log(P.hat.Omega)) - sum((1-A.Omega)*log(1-P.hat.Omega))
        print(proc.time() - ptm)
#### Degree correct model
        V <- U.approx
        print("DCSBM calculation")
        ptm <- proc.time()
        #V.norms <- apply(V,1,function(x) sqrt(sum(x^2)))
        if(k==1) {V.norms <- as.numeric(abs(V))}else{
            V.norms <- apply(V,1,function(x) sqrt(sum(x^2)))
        }

        iso.index <- which(V.norms==0)
        Psi <- V.norms
        Psi <- Psi / max(V.norms)
        inv.V.norms <- 1/V.norms
        inv.V.norms[iso.index] <- 1

        V.normalized <- diag(as.numeric(inv.V.norms))%*%V

        #V.norms <- apply(V,1,function(x) sqrt(sum(x^2)))
        #Psi <- V.norms
        #Psi <- Psi / max(V.norms)
        #V.normalized <- diag(1/V.norms)%*%V
        #Psi.outer <- outer(Psi,Psi)
        if(k==1){
        if(dc.est>1){
            B <- sum(A.new,na.rm=TRUE)+0.01

            partial.d <- colSums(A.new,na.rm=TRUE)
            partial.gd <- B
            phi <- rep(0,n)
            B.g <- partial.gd
            phi <- as.numeric(partial.d/B.g)
            B <- B/p.sample
            P.hat <- t(t(matrix(B,n,n)*phi)*phi)
            #P.hat <- diag(phi)%*%matrix(B,n,n)%*%diag(phi)
            diag(P.hat) <- 0
        }
            dc.block.sq.err[k] <- sum((pb-A[Omega])^2)
            P.hat.Omega <- P.hat[Omega]
            A.Omega <- A[Omega]
            P.hat.Omega[P.hat.Omega < 1e-6] <- 1e-6
            P.hat.Omega[P.hat.Omega > (1-1e-6)] <- 1-1e-6

            dc.loglike[k] <- -sum(A.Omega*log(P.hat.Omega)) - sum((1-A.Omega)*log(1-P.hat.Omega))


        }else{
        km <- kmeans(V.normalized,centers=k,nstart=30,iter.max=30)
        if(dc.est>1){
            B <- matrix(0,k,k)
            Theta <- matrix(0,n,k)
            for(i in 1:k){
                for(j in 1:k){
                N.i <- which(km$cluster==i)
                N.j <- which(km$cluster==j)
                B[i,j] <- sum(A.new[N.i,N.j],na.rm=TRUE)+0.01
                }
                Theta[N.i,i] <- 1
            }
            Theta <- Matrix(Theta,sparse=TRUE)
            partial.d <- colSums(A.new,na.rm=TRUE)
            partial.gd <- colSums(B)
            phi <- rep(0,n)
            B.g <- Theta%*%partial.gd
            phi <- as.numeric(partial.d/B.g)
            B <- B/p.sample
            tmp.int.mat <- Theta*phi
            P.hat <-as.matrix(tmp.int.mat%*%B%*%t(tmp.int.mat))
            #P.hat <- diag(phi)%*%Theta%*%B%*%t(Theta)%*%diag(phi)
            diag(P.hat) <- 0
        }
        dc.block.sq.err[k] <- sum((P.hat[Omega]-A[Omega])^2)
        P.hat.Omega <- P.hat[Omega]
        A.Omega <- A[Omega]
        P.hat.Omega[P.hat.Omega < 1e-6] <- 1e-6
        P.hat.Omega[P.hat.Omega > (1-1e-6)] <- 1-1e-6
        dc.loglike[k] <- -sum(A.Omega*log(P.hat.Omega)) - sum((1-A.Omega)*log(1-P.hat.Omega))
    }
        print(proc.time() - ptm)



    }
    return(list(impute.sq.err=impute.sq.err,block.sq.err=block.sq.err,loglike=loglike,roc.auc=roc.auc,no.edge=no.edge,dc.block.sq.err=dc.block.sq.err,dc.loglike=dc.loglike,bin.dev=bin.dev))
}











#### Using AUC to evaluate intrinsic rank
ECV.undirected.Rank <- function(A,max.K,B=3,holdout.p=0.1,soft=FALSE,weighted=FALSE){
    n <- nrow(A)
    edge.index <- which(upper.tri(A))
    edge.n <- length(edge.index)
    #edge.index <- which(upper.tri(A))
    #edge.n <- length(edge.index)

    holdout.index.list <- list()

    holdout.n <- floor(holdout.p*edge.n)

    for(j in 1:B){
        holdout.index.list[[j]] <- sample(x=edge.n,size=holdout.n)
    }
    result <- lapply(holdout.index.list,missing.Rank.fast.all.both,A=A,max.K=max.K,p.sample=1-holdout.p,weighted)
    sse.mat <- roc.auc.mat <- matrix(0,nrow=B,ncol=max.K)

    for(b in 1:B){
        roc.auc.mat[b,] <- result[[b]]$roc.auc
        sse.mat[b,] <- result[[b]]$sse
    }

    auc.seq <- colMeans(roc.auc.mat)
    auc.sd <- apply(roc.auc.mat,2,sd)/sqrt(B)
    sse.seq <- colMeans(sse.mat)
    sse.sd <- apply(sse.mat,2,sd)/sqrt(B)
    return(list(rank=which.max(auc.seq),auc=auc.seq,auc.sd=auc.sd,sse=sse.seq,sse.sd=sse.sd))
    #return(list(rank=which.max(auc.seq),auc=auc.seq))
}







missing.Rank.fast.all.both <- function(holdout.index,A,max.K,p.sample=NULL,weighted=FALSE){
    n <- nrow(A)
    edge.index <- which(upper.tri(A))
    edge.n <- length(edge.index)
    A.new <- matrix(0,n,n)
    A.new[upper.tri(A.new)] <- A[edge.index]
    A.new[edge.index[holdout.index]] <- NA
    A.new <- A.new + t(A.new)
    if(is.null(p.sample)){
        p.sample <- 1-length(holdout.index)/n^2
    }
    degrees <- colSums(A.new,na.rm=TRUE)
    no.edge <- 0
    no.edge <- sum(degrees==0)

    Omega <- which(is.na(A.new))
    imputed.A <- list()
    sse <- roc.auc <- rep(0,max.K)
    SVD.result <- iter.SVD.core.fast.all(A.new,max.K,fast=TRUE,p.sample=p.sample)

    for(k in 1:max.K){
        print(k)
        tmp.est <- SVD.result[[k]]
        A.approx <- tmp.est$A.thr
        response <- A[Omega]
        predictors <- A.approx[Omega]
        if(!weighted){
        aa <- AUC::roc(predictions=predictors,labels=factor(response))
        roc.auc[k] <- auc(aa)
        }
        sse[k] <- mean((response-predictors)^2)
        imputed.A[[k]] <- A.approx
    }
    return(list(roc.auc=roc.auc,imputed.A=imputed.A,Omega=Omega,sse=sse))
}
