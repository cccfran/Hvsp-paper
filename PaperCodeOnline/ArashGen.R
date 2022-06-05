## lambda: average degree
## n: network size
## beta: in-out ratio
## w:block specific degree
## Pi: block size proportion
## rho: 1-rho is hub-proportion. If rho >0, it indicates that we use degree corrected model.
## simple: if rho >0, then fix the degree of hubs and normal nodes by 1 and 0.2. If simple == FALSE and rho >0, go to more general degree distribution
## power: if TRUE, generate degrees according to exponential with (3*log(10)) as the rate, truncated between [0.2,1]. This is approximatly power law. If FALSE, generate uniform degrees over [0.2,1].

ArashSBM <- function(lambda,n,beta=0,K=3,w=rep(1,K),Pi=rep(1,K)/K,rho=0,simple=TRUE,power=TRUE){
    P0 <- diag(w)
    if(beta > 0){
        P0 <- matrix(1,K,K)
        diag(P0) <- w/beta
    }
    Pi.vec <- matrix(Pi,ncol=1)
    P <- lambda*P0/((n-1)*as.numeric(t(Pi.vec)%*%P0%*%Pi.vec)*(rho*0.2+(1-rho))^2)
    if((rho >0) && (!simple) && (!power)){ P <- lambda*P0/((n-1)*as.numeric(t(Pi.vec)%*%P0%*%Pi.vec)*(0.6)^2)   }
    if((rho >0) && (!simple) && (power)){ P <- lambda*P0/((n-1)*as.numeric(t(Pi.vec)%*%P0%*%Pi.vec)*((0.2364)^2))   }

    M <- matrix(0,n,K)
    membership <- sample(x=K,size=n,replace=TRUE,prob=Pi)

        M[cbind(1:n,membership)] <- 1

    A.bar <- M%*%P%*%t(M)
    node.degree <- rep(1,n)
    if(rho>0){
    if(simple){
        node.degree[runif(n)<rho] <- 0.2
    }else{
        if(power==FALSE){
            node.degree <- runif(n)*0.8 + 0.2
        }else{
            node.degree <- rexp(n,rate=3*log(10))
            node.degree[node.degree<0.2] <- 0.2
            node.degree[node.degree>1] <- 1
        }
    }}
    DD <- diag(node.degree)
    A.bar <- DD%*%A.bar%*%DD
    diag(A.bar) <- 0
    upper.index <- which(upper.tri(A.bar))
    upper.p <- A.bar[upper.index]
    upper.u <- runif(n=length(upper.p))
    upper.A <- rep(0,length(upper.p))
    upper.A[upper.u < upper.p] <- 1
    A <- matrix(0,n,n)
    A[upper.index] <- upper.A
    A <- A + t(A)
    return(A)
}


ArashSBM.full <- function(lambda,n,beta=0,K=3,w=rep(1,K),Pi=rep(1,K)/K,rho=0,simple=TRUE,power=TRUE,B=NULL){
    if(is.null(B)){
    P0 <- diag(w)
    if(beta > 0){
        P0 <- matrix(1,K,K)
        diag(P0) <- w/beta
    }
    Pi.vec <- matrix(Pi,ncol=1)
    P <- lambda*P0/((n-1)*as.numeric(t(Pi.vec)%*%P0%*%Pi.vec)*(rho*0.2+(1-rho))^2)
    if((rho >0) && (!simple) && (!power)){ P <- lambda*P0/((n-1)*as.numeric(t(Pi.vec)%*%P0%*%Pi.vec)*(0.6)^2)   }
    if((rho >0) && (!simple) && (power)){ P <- lambda*P0/((n-1)*as.numeric(t(Pi.vec)%*%P0%*%Pi.vec)*((0.2364)^2))   }
    }else{
        P <- B
    }
    M <- matrix(0,n,K)
    membership <- sample(x=K,size=n,replace=TRUE,prob=Pi)
        M[cbind(1:n,membership)] <- 1
    A.bar <- M%*%P%*%t(M)
    node.degree <- rep(1,n)
    if(rho>0){
    if(simple){
        node.degree[runif(n)<rho] <- 0.2
    }else{
        if(power==FALSE){
            node.degree <- runif(n)*0.8 + 0.2
        }else{
            node.degree <- rexp(n,rate=3*log(10))
            node.degree[node.degree<0.2] <- 0.2
            node.degree[node.degree>1] <- 1
        }
    }}
    DD <- diag(node.degree)
    A.bar <- DD%*%A.bar%*%DD
    PP <- A.bar
    diag(A.bar) <- 0
    upper.index <- which(upper.tri(A.bar))
    upper.p <- A.bar[upper.index]
    upper.u <- runif(n=length(upper.p))
    upper.A <- rep(0,length(upper.p))
    upper.A[upper.u < upper.p] <- 1
    A <- matrix(0,n,n)
    A[upper.index] <- upper.A
    A <- A + t(A)
    return(list(A=A,g=membership,B=P))
}






ArashSBM.pure <- function(lambda,n,beta=0,K=3,w=rep(1,K),Pi=rep(1,K)/K,B=NULL,membership=NULL){
    if(is.null(B)){
    P0 <- diag(w)
    if(beta > 0){
        P0 <- matrix(1,K,K)
        diag(P0) <- w/beta
    }
    Pi.vec <- matrix(Pi,ncol=1)
    P <- lambda*P0/((n-1)*as.numeric(t(Pi.vec)%*%P0%*%Pi.vec))
    }else{
        P <- B
        if(!is.na(lambda)){
        expect.d <- n*as.numeric(matrix(Pi,nrow=1)%*%B%*%matrix(Pi,ncol=1))
        P <- P*lambda/expect.d
        }
    }
    M <- matrix(0,n,K)
    if(is.null(membership)){
    membership <- sample(x=K,size=n,replace=TRUE,prob=Pi)
    }
        M[cbind(1:n,membership)] <- 1
    print(dim(M))
    A.bar <- M%*%P%*%t(M)
    PP <- A.bar
    diag(A.bar) <- 0
    upper.index <- which(upper.tri(A.bar))
    upper.p <- A.bar[upper.index]
    upper.u <- runif(n=length(upper.p))
    upper.A <- rep(0,length(upper.p))
    upper.A[upper.u < upper.p] <- 1
    A <- matrix(0,n,n)
    A[upper.index] <- upper.A
    A <- A + t(A)
    return(list(A=A,g=membership,B=P))
}





BlurredSBM <- function(lambda,n,beta=0,K=3,w=rep(1,K),Pi=rep(1,K)/K,rho=0.2){
    P0 <- diag(w)
    if(beta > 0){
        P0 <- matrix(1,K,K)
        diag(P0) <- w/beta
    }
    theta <- diag(runif(K)*(1-rho)+rho)
    for(i in 1:(K-1)){
        theta[i,i+1] <- runif(1)*rho
        #theta[i,i+1] <- runif(1)*rho/2

    }
    P0 <- theta%*%P0%*%t(theta)
    Pi.vec <- matrix(Pi,ncol=1)
    P <- lambda*P0/((n-1)*as.numeric(t(Pi.vec)%*%P0%*%Pi.vec))
    M <- matrix(0,n,K)
    membership <- sample(x=K,size=n,replace=TRUE,prob=Pi)
    for(i in 1:n){
        M[i,membership[i]] <- 1
    }
    A.bar <- M%*%P%*%t(M)
    A.bar.2 <- A.bar
    diag(A.bar) <- 0
    upper.index <- which(upper.tri(A.bar))
    upper.p <- A.bar[upper.index]
    upper.u <- runif(n=length(upper.p))
    upper.A <- rep(0,length(upper.p))
    upper.A[upper.u < upper.p] <- 1
    A <- matrix(0,n,n)
    A[upper.index] <- upper.A
    A <- A + t(A)
    return(list(A=A,P=P,A.bar=A.bar.2))
}


GenericLowRank <- function(n,K){
    Z <- matrix(abs(runif(n*K)),nrow=n,ncol=K)
    S <- Z%*%t(Z)
    P <- S/max(S)

    A <- matrix(0,n,n)
    upper.tri.index <- which(upper.tri(P))
    prob.seq <- P[upper.tri.index]
    T <- length(prob.seq)
    binary.seq <- rep(0,T)
    for(t in 1:T){
        binary.seq[t] <- rbinom(1,1,prob=prob.seq[t])
    }
    A[upper.tri.index] <- binary.seq
    A <- A + t(A)
    return(A)
}





ArashSBMX <- function(lambda,n,beta=0,K=3,w=rep(1,K),Pi=rep(1,K)/K,rho=0,simple=TRUE,power=TRUE,p=1,sigma=1){
    P0 <- diag(w)
    if(beta > 0){
        P0 <- matrix(1,K,K)
        diag(P0) <- w/beta
    }
    Pi.vec <- matrix(Pi,ncol=1)
    P <- lambda*P0/((n-1)*as.numeric(t(Pi.vec)%*%P0%*%Pi.vec)*(rho*0.2+(1-rho))^2)
    if((rho >0) && (!simple) && (!power)){ P <- lambda*P0/((n-1)*as.numeric(t(Pi.vec)%*%P0%*%Pi.vec)*(0.6)^2)   }
    if((rho >0) && (!simple) && (power)){ P <- lambda*P0/((n-1)*as.numeric(t(Pi.vec)%*%P0%*%Pi.vec)*((0.2364)^2))   }

    M <- matrix(0,n,K)
    membership <- sample(x=K,size=n,replace=TRUE,prob=Pi)
    for(i in 1:n){
        M[i,membership[i]] <- 1
    }
    A.bar <- M%*%P%*%t(M)
    node.degree <- rep(1,n)
    if(rho>0){
    if(simple){
        node.degree[runif(n)<rho] <- 0.2
    }else{
        if(power==FALSE){
            node.degree <- runif(n)*0.8 + 0.2
        }else{
            node.degree <- rexp(n,rate=3*log(10))
            node.degree[node.degree<0.2] <- 0.2
            node.degree[node.degree>1] <- 1
        }
    }}
    DD <- diag(node.degree)
    A.bar <- DD%*%A.bar%*%DD
    diag(A.bar) <- 0
    upper.index <- which(upper.tri(A.bar))
    upper.p <- A.bar[upper.index]
    upper.u <- runif(n=length(upper.p))
    upper.A <- rep(0,length(upper.p))
    upper.A[upper.u < upper.p] <- 1
    A <- matrix(0,n,n)
    A[upper.index] <- upper.A
    A <- A + t(A)


    X <- matrix(rnorm(n*p),nrow=n,ncol=p)*sigma
    MU <- seq(-1,1,length.out=K)
    for(i in 1:n){
        if(runif(1)>=0.2){
        X[i,] <- X[i,] + MU[membership[i]]
        }
    }

    return(list(A=A,X=X,g=membership))
}
