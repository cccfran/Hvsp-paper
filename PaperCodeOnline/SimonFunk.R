##### A batch gradient descent implementation of Simon Funk's low-rank SVD

library('Matrix')
library('irlba')



SF.SVD.core <- function(A,K,learning.rate=0.1,tol=1e-6,max.iter=100,sparse=TRUE,init=NULL,verbose=FALSE){
    if(sparse) A <- Matrix(A,sparse=TRUE)
    Omega <- which(is.na(A))
    A.col.means <- colMeans(A,na.rm=TRUE)
    A.impute <- A
    n <- nrow(A)
    p <- ncol(A)
    if(is.null(init)){
        A.impute[Omega] <- runif(n=length(Omega))
        init.SVD <- svd(A.impute,nu=K,nv=K)
    }else{
        init.SVD <- init
    }
    ## init.SVD <- irlba(A,nu=K,nv=K) ## if you are working on large problems
    U.old <- matrix(init.SVD$u[,1:K],ncol=K)%*%diag(sqrt(init.SVD$d[1:K]))
    V.old <- matrix(init.SVD$v[,1:K],ncol=K)%*%diag(sqrt(init.SVD$d[1:K]))
    A.old <- U.old %*% t(V.old)
    R <- A - A.old

    R[Omega] <- 0
    if(verbose) print(norm(R))
    ### begin iteration
    flag <- 0
    iter <- 0
    err.seq <- Inf
    shrink <- 0
    while((iter < max.iter) && (flag != 1)){
        iter <- iter + 1
        #print(learning.rate)
        R <- Matrix(R,sparse=TRUE)
        U.new <- U.old + 2*learning.rate*R%*%V.old
        t.V.new <- t(V.old) + 2*learning.rate*t(U.old)%*% R

        A.new <- U.new%*%t.V.new
        if(norm(A.new-A.old,"f")< tol*norm(A.old,"f")) flag <- 1
        A.new <- U.new %*%t.V.new
        tmp.R <- A - A.new
        tmp.R[Omega] <- 0
        err <- norm(tmp.R,"f")
        if(err>err.seq[length(err.seq)]){
            learning.rate <- learning.rate/2
            shrink <- shrink+1
            iter <- iter - 1
            if(shrink>30) {
                print("trapped in local area!")
                break
            }
            if(verbose) print("Stepsize shrunken!")
            next
        }
        if(verbose) print(err)
        shrink <- 0
        err.seq <- c(err.seq,err)
        R <- tmp.R
        U.old <- U.new
        V.old <- t(t.V.new)
        A.old <- A.new
    }
    trapped <- ifelse(shrink>30,1,0)
    return(list(U=U.old,V=V.old,A=A.new,iter=iter,SVD=svd(A.new,nu=K,nv=K),err.seq=err.seq,step=learning.rate,trapped=trapped))
}


SF.SVD <- function(A,K,learning.rate=0.1,tol=1e-7,max.iter=100,sparse=TRUE,init=NULL,verbose=FALSE,nstart=1){
	if(is.null(init)){
    approx.list <- list()
    err.seq <- rep(0,nstart)
    for(D in 1:nstart){
        approx.list[[D]] <- SF.SVD.core(A=A,K=K,learning.rate=learning.rate,tol=tol,max.iter=max.iter,sparse=sparse,init=init,verbose=verbose)
        err.seq[D] <- min(approx.list[[D]]$err.seq)
    }
    select.index <- which.min(err.seq)
    return(approx.list[[select.index]])
    }
    return(SF.SVD.core(A=A,K=K,learning.rate=learning.rate,tol=tol,max.iter=max.iter,sparse=sparse,init=init,verbose=verbose))
}



#### An alternative: using coordinate descent idea, which turns out to be a iterative SVD



iter.SVD.core <- function(A,K,tol=1e-5,max.iter=100,sparse=TRUE,init=NULL,verbose=FALSE,tau=0,fast=FALSE,p.sample=1,kappa=NULL){
    if(sparse) A <- Matrix(A,sparse=TRUE)
     avg.p <- mean(as.numeric(A),na.rm=TRUE)
   if(is.null(kappa))kappa <- 2/avg.p
    cap <- 1#kappa*avg.p
    if(cap>1-tau) cap <- 1-tau
    if(fast){
        if(verbose) print("Matrix completion with fast approximation!")
        A[which(is.na(A))] <- 0
        A <- A/p.sample
        #svd.new <- svd(A,nu=K,nv=K)
        svd.new <- irlba(A,nu=K,nv=K)
        if(K==1){ A.new <- svd.new$d[1]*matrix(svd.new$u,ncol=1)%*%t(matrix(svd.new$v,ncol=1))}else{
        A.new <- svd.new$u%*%(t(svd.new$v)*svd.new$d[1:K])}
        A.new[A.new < 0+tau] <- 0+tau
        A.new[A.new >cap] <- cap
        return(list(iter=NA,SVD=svd.new,A=A.new,err.seq=NA))
    }
    #if(sparse) A <- Matrix(A,sparse=TRUE)
    Omega <- which(is.na(A))
    A.col.means <- colMeans(A,na.rm=TRUE)
    #avg.p <- mean(as.numeric(A),na.rm=TRUE)
    A.impute <- A
    n <- nrow(A)
    p <- ncol(A)
    if(is.null(init)){
        A.impute[Omega] <- runif(n=length(Omega))
        init.SVD <- irlba(A.impute,nu=K,nv=K)
    }else{
        init.SVD <- init
    }
    #print(init.SVD$u)
    ## init.SVD <- irlba(A,nu=K,nv=K) ## if you are working on large problems
    if(K==1){U.old <- V.old <- matrix(0,n,1)}else{
    if(K==2){
    U.old <- matrix(init.SVD$u[,1:(K-1)],ncol=K-1)*(sqrt(init.SVD$d[1:(K-1)]))
    V.old <- matrix(init.SVD$v[,1:(K-1)],ncol=K-1)*(sqrt(init.SVD$d[1:(K-1)]))
    }else{
    #print(init.SVD$u)
    U.old <- matrix(init.SVD$u[,1:(K-1)],ncol=K-1)%*%diag(sqrt(init.SVD$d[1:(K-1)]))
    V.old <- matrix(init.SVD$v[,1:(K-1)],ncol=K-1)%*%diag(sqrt(init.SVD$d[1:(K-1)]))
    }
    }
    A.old <- U.old %*% t(V.old)
    R <- A - A.old

    R[Omega] <- 0
    if(verbose) print(norm(R))
    A.impute <- A
    A.impute[Omega] <- A.old[Omega]
    ### begin iteration
    flag <- 0
    iter <- 0
    err.seq <- norm(R,"F")
    shrink <- 0
    while((iter < max.iter) && (flag != 1)){
        #print(iter)
        iter <- iter + 1
        svd.new <- irlba(A.impute,nu=K,nv=K)
        if(K==1){ A.new <- svd.new$d[1]*matrix(svd.new$u,ncol=1)%*%t(matrix(svd.new$v,ncol=1))}else{
        A.new <- svd.new$u%*%diag(svd.new$d)%*%t(svd.new$v)}
        A.new[A.new < 0+tau] <- 0+tau
        A.new[A.new >cap] <- cap
        A.impute[Omega] <- A.new[Omega]
        A.old <- A.new
        R <- A.impute - A.new
        err <- norm(R,"F")
        if(verbose) print(err)
        err.seq <- c(err.seq,err)
        if(abs(err.seq[iter+1]-err.seq[iter])<tol*err.seq[iter]) flag <- 1
     }
    #print(iter)
    return(list(iter=iter,SVD=svd.new,A=A.new,err.seq=err.seq))
}
# source("LeiFunc.R")
# A <- LeiSim1(0.1,K=8,n=1000,n1=300)
# system.time(tt2<-SF.SVD.core(A,K=5,verbose=FALSE))
# system.time(tt3<-iter.SVD.core(A,K=5,verbose=FALSE))


ExactRankApprox <- function(X,K,maxit){
    #### we want to use the hardImpute for approximation, with initialization by softImpute.
    #### initialization
    lam0 <- lambda0(X)
    print("get lambda0")
    lam <- lam0/2
    tmp.init <- softImpute(x=X,rank.max=K+1,type="svd",lambda=lam)
    lam <- lam/2
    solved <- FALSE
    while(!solved){
        if(sum(tmp.init$d>1e-6)>=K){
            solved <- TRUE
        }
        tmp.init <- softImpute(x=X,rank.max=K+1,type="svd",lambda=lam,warm.start=tmp.init)
        lam <- lam/2
        print("Shrinking!")
    }
    print("Start hard impute!")
    SVD <- softImpute(X,rank.max=K,maxit=maxit,type="svd",warm.start=tmp.init)
    return(SVD)
}







iter.SVD.core.L <- function(A,K,tol=1e-5,max.iter=100,sparse=TRUE,init=NULL,verbose=FALSE,tau=0){
    if(sparse) A <- Matrix(A,sparse=TRUE)
    Omega <- which(is.na(A))
    A.col.means <- colMeans(A,na.rm=TRUE)
    A.impute <- A
    n <- nrow(A)
    p <- ncol(A)
    if(is.null(init)){
        A.impute[Omega] <- runif(n=length(Omega))
        init.SVD <- svd(A.impute,nu=K,nv=K)
    }else{
        init.SVD <- init
    }
    #print(init.SVD$u)
    ## init.SVD <- irlba(A,nu=K,nv=K) ## if you are working on large problems
    if(K==1){U.old <- V.old <- matrix(0,n,1)}else{
    if(K==2){
    U.old <- matrix(init.SVD$u[,1:(K-1)],ncol=K-1)*(sqrt(init.SVD$d[1:(K-1)]))
    V.old <- matrix(init.SVD$v[,1:(K-1)],ncol=K-1)*(sqrt(init.SVD$d[1:(K-1)]))
    }else{
    #print(init.SVD$u)
    U.old <- matrix(init.SVD$u[,1:(K-1)],ncol=K-1)%*%diag(sqrt(init.SVD$d[1:(K-1)]))
    V.old <- matrix(init.SVD$v[,1:(K-1)],ncol=K-1)%*%diag(sqrt(init.SVD$d[1:(K-1)]))
    }
    }
    A.old <- U.old %*% t(V.old)
    R <- A - A.old

    R[Omega] <- 0
    if(verbose) print(norm(R))
    A.impute <- A
    A.impute[Omega] <- A.old[Omega]
    ### begin iteration
    flag <- 0
    iter <- 0
    err.seq <- norm(R,"F")
    shrink <- 0
    while((iter < max.iter) && (flag != 1)){
        #print(iter)
        iter <- iter + 1
        svd.new <- irlba(A.impute,nu=K,nv=K)
        if(K==1){ A.new <- svd.new$d[1]*matrix(svd.new$u,ncol=1)%*%t(matrix(svd.new$v,ncol=1))}else{
        A.new <- svd.new$u%*%diag(svd.new$d)%*%t(svd.new$v)}
        A.new[A.new < -1+tau] <- -1+tau
        A.new[A.new >1+tau] <- 1+tau
        A.impute[Omega] <- A.new[Omega]
        A.old <- A.new
        R <- A.impute - A.new
        err <- norm(R,"F")
        if(verbose) print(err)
        err.seq <- c(err.seq,err)
        if(abs(err.seq[iter+1]-err.seq[iter])<tol*err.seq[iter]) flag <- 1
     }
    #print(iter)
    return(list(iter=iter,SVD=svd.new,A=A.new,err.seq=err.seq))
}












iter.SVD.core.fast.all <- function(A,Kmax,tol=1e-5,max.iter=100,sparse=TRUE,init=NULL,verbose=FALSE,tau=0,fast=FALSE,p.sample=1){
    if(sparse) A <- Matrix(A,sparse=TRUE)
     avg.p <- mean(as.numeric(A),na.rm=TRUE)
    cap <- 1#kappa*avg.p
        A[which(is.na(A))] <- 0
        A <- A/p.sample
        #svd.new <- svd(A,nu=K,nv=K)
     #print("begin SVD")
        svd.new <- irlba(A,nu=Kmax,nv=Kmax)
    #print("end SVD")
        result <- list()
        for(K in 1:Kmax){
            print(K)
        if(K==1){
            A.new <- svd.new$d[1]*matrix(svd.new$u[,1],ncol=1)%*%t(matrix(svd.new$v[,1],ncol=1))
              }else{
        A.new <- A.new + svd.new$d[K]*matrix(svd.new$u[,K],ncol=1)%*%t(matrix(svd.new$v[,K],ncol=1))
        }
        A.new.thr <- A.new
        A.new.thr[A.new < 0+tau] <- 0+tau
        A.new.thr[A.new >cap] <- cap
        
        tmp.SVD <- list(u=svd.new$u[,1:K],v=svd.new$v[,1:K],d=svd.new$d[1:K])
        result[[K]] <- list(iter=NA,SVD=tmp.SVD,A=A.new,err.seq=NA,A.thr=A.new.thr)
        }
        return(result)

}
