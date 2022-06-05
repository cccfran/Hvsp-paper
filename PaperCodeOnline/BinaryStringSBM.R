BinaryStringSBM <- function(n,d,beta=0.3,lambda){
    K <- 2^d
    #outin <- beta
    #beta <- (2*K*outin/(K-1))^((1/((d-1))))
    ## generate binary strings
    b.list <- list()
    for(k in 1:K){
        b.list[[k]] <- as.character(intToBits(k-1))[d:1]
    }
    ## construct B
    comm.sim.mat <- B <- matrix(0,K,K)
    for(i in 1:(K-1)){
        for(j in (i+1):K){
            s <- Binary.Similarity(b.list[[i]],b.list[[j]])-1
            comm.sim.mat[i,j] <- s+1
            
            dd <- d-s
            B[i,j] <- exp(dd*log(beta))
        }
    }
    B <- B+t(B)
    comm.sim.mat <- comm.sim.mat+t(comm.sim.mat)
    diag(comm.sim.mat) <- d+1
    diag(B) <- 1
    w <- floor(n/K)
    g <- c(rep(seq(1,K),each=w),rep(K,n-w*K))
    Z <- matrix(0,n,K)
    Z[cbind(1:n,g)] <- 1
    P <- Z%*%B%*%t(Z)
    P <- P*lambda/mean(colSums(P))
    node.sim.mat <- Z%*%comm.sim.mat%*%t(Z)
    diag(node.sim.mat) <- 0
    print(paste("Within community expected edges: ",P[1,1]*w,sep=""))
    return(list(b.list=b.list,B=B,label=g,P=P,comm.sim.mat=comm.sim.mat,node.sim.mat=node.sim.mat))
}



BinaryStringSBM.general.seq <- function(n,d,a.seq,lambda,alpha=NULL){
    K <- 2^d
    #outin <- beta
    #beta <- (2*K*outin/(K-1))^((1/((d-1))))
    ## generate binary strings
    b.list <- list()
    for(k in 1:K){
        b.list[[k]] <- as.character(intToBits(k-1))[d:1]
    }
    ## construct B
    comm.sim.mat <- B <- matrix(0,K,K)
    for(i in 1:(K-1)){
        for(j in (i+1):K){
            s <- Binary.Similarity(b.list[[i]],b.list[[j]])-1
            comm.sim.mat[i,j] <- s+1            
        }
    }
    comm.sim.mat <- comm.sim.mat+t(comm.sim.mat)
    diag(comm.sim.mat) <- d+1
    B[1:(K^2)] <- a.seq[d+2-as.numeric(comm.sim.mat)]
    w <- floor(n/K)
    g <- c(rep(seq(1,K),each=w),rep(K,n-w*K))
    Z <- matrix(0,n,K)
    Z[cbind(1:n,g)] <- 1
    P <- Z%*%B%*%t(Z)
    if(is.null(alpha)){
    P <- P*lambda/mean(colSums(P))
    }else{
        P <- P*alpha
    }
    node.sim.mat <- Z%*%comm.sim.mat%*%t(Z)
    diag(node.sim.mat) <- 0
    print(paste("Within community expected edges: ",P[1,1]*w,sep=""))
    return(list(b.list=b.list,B=B,label=g,P=P,comm.sim.mat=comm.sim.mat,node.sim.mat=node.sim.mat))
}




BinaryStringSBM.simple.imbalance <- function(n,d,beta=0.3,lambda,rho=0){
    K <- 2^d
    outin <- beta
    beta <- (2*K*outin/(K-1))^((1/((d-1))))
    ## generate binary strings
    b.list <- list()
    for(k in 1:K){
        b.list[[k]] <- as.character(intToBits(k-1))[d:1]
    }
    ## construct B
    comm.sim.mat <- B <- matrix(0,K,K)
    for(i in 1:(K-1)){
        for(j in (i+1):K){
            s <- Binary.Similarity(b.list[[i]],b.list[[j]])-1
            comm.sim.mat[i,j] <- s+1
            
            dd <- d-s
            B[i,j] <- exp(dd*log(beta))
        }
    }
    B <- B+t(B)
    comm.sim.mat <- comm.sim.mat+t(comm.sim.mat)
    diag(comm.sim.mat) <- d+1
    diag(B) <- 1
    size.weight <- (1:K)^rho
    size.weight <- size.weight/sum(size.weight)
    w <- floor(n/K)
    size.n <- unlist(lapply(n*size.weight,floor))
    size.index <- c(0,cumsum(size.n))
    g <- rep(K,n)
    for(k in 1:(K-1)){
        g[(size.index[k]+1):size.index[k+1]] <- k
    }
    #g <- c(rep(seq(1,K),each=w),rep(K,n-w*K))
    Z <- matrix(0,n,K)
    Z[cbind(1:n,g)] <- 1
    P <- Z%*%B%*%t(Z)
    P <- P*lambda/mean(colSums(P))
    node.sim.mat <- Z%*%comm.sim.mat%*%t(Z)
    diag(node.sim.mat) <- 0
    print(paste("Within community expected edges: ",P[1,1]*w,sep=""))
    return(list(b.list=b.list,B=B,label=g,P=P,comm.sim.mat=comm.sim.mat,node.sim.mat=node.sim.mat))
}



#dt <- BinaryStringSBM(2000,5,beta=0.3,lambda=10)


BinStringShoot <- function(x,dmax,rho=0){
    result <- list()
    if(length(x)==dmax){
        result[[1]] <- x
        return(result)
    }
    shoot <- as.numeric(runif(1)<rho)
    if(shoot==1){
        result[[1]] <- x
        return(result)
    }else{
        
    }
}



BinString.RS.shoot <- function(V,rho,split.flag=TRUE){
    prob.s <- rho
    n <- length(V)
    randUnif <- runif(1)
    next.split.flag <- FALSE
    if(randUnif<prob.s){
        next.split.flag <- TRUE
   }
    if(length(V)==1){
        split.flag=FALSE
    }
    if(split.flag){
        n1 <- n/2
        sub.BS1 <- BinString.RS.shoot(V[1:n1],rho,next.split.flag)
        sub.BS2 <- BinString.RS.shoot(V[(n1+1):n],rho,next.split.flag)
        V <- list(sub.BS1$V,sub.BS2$V)
        left.path <- paste("L",sub.BS1$tree.path,sep="/")
        right.path <- paste("R",sub.BS2$tree.path,sep="/")
        tree.path <- c(left.path,right.path)
    }else{
        print("Stop splitting")
        tree.path <- as.character(V)
        tree.path2 <- c("",tree.path)
        V <- list(V)
    }
    return(list(tree.path=tree.path))
}



BinaryStringSBM.RandomMerging <- function(n,d,beta=0.3,lambda,rho = 0){
    V <- 1:(2^d)
    tt <- BinString.RS.shoot(V,rho=rho,split.flag=TRUE)
    tt2 <- BinString.RS.shoot(V,rho=1,split.flag=TRUE)
    comm.tree.path <- strsplit(tt$tree.path,"/")
    comm.tree.path <- lapply(comm.tree.path,function(x) return(x[1:(length(x)-1)]))
    comm.tree.path2 <- strsplit(tt2$tree.path,"/")
    comm.tree.path2 <- lapply(comm.tree.path2,function(x) return(x[1:(length(x)-1)]))
    b.list <- comm.tree.path2
    K <- length(b.list)
    outin <- beta
    beta <- (2*K*outin/(K-1))^((1/((d-1))))
    #beta <- (outin)^((1/(sqrt(d-1))))
    ## generate binary strings
    ## construct B
    comm.sim.mat <- B <- matrix(0,K,K)
    for(i in 1:(K-1)){
        for(j in (i+1):K){
            s <- Binary.Similarity(b.list[[i]],b.list[[j]])-1
            comm.sim.mat[i,j] <- s+1
            
            dd <- d-s
            B[i,j] <- exp(dd*log(beta))
        }
    }
    B <- B+t(B)
    comm.sim.mat <- comm.sim.mat+t(comm.sim.mat)
    diag(comm.sim.mat) <- d+1
    diag(B) <- 1

#### Now begin to merge
####
    merge.list <- list()
    b.list.unique <- unique(comm.tree.path)
    merge.list <- lapply(b.list.unique,function(x){
               md <- length(x)
               result <- NULL
               for(ii in 1:length(b.list)){
                   if(sum(b.list[[ii]][1:md]!=x)==0){
                       result <- c(result,ii)
                   }
               }
               return(result)
           })
    ## construct true comm.sim.mat and B
    K <- length(b.list.unique)
    true.C <- true.B <- matrix(0,K,K)
    for(i in 1:K){
        for(j in 1:K){
            true.C[i,j] <- mean(comm.sim.mat[merge.list[[i]],merge.list[[j]]])
            true.B[i,j] <- mean(B[merge.list[[i]],merge.list[[j]]])
        }
    }
    w <- floor(n/length(b.list))
    size.ratio <- unlist(lapply(merge.list,length))
    size.index <- cumsum(size.ratio)
    size.index <- c(0,size.index)
    g <- rep(K,n)
    for(k in 1:K){
        g[(size.index[k]*w+1):(size.index[k+1]*w)] <- k
    }
    Z <- matrix(0,n,K)
    Z[cbind(1:n,g)] <- 1
    P <- Z%*%true.B%*%t(Z)
    P <- P*lambda/mean(colSums(P))
    node.sim.mat <- Z%*%true.C%*%t(Z)
    diag(node.sim.mat) <- 0
    exp.out.in <- rep(0,K)
    for(k in 1:K){
        exp.out.in[k] <- (sum(true.B[k,-k]*size.ratio[-k])/sum(size.ratio[-k]))/((true.B[k,k]*size.ratio[k])/(size.ratio[k]))
    }
    print(paste("Expected out-in ratio:  ",mean(exp.out.in),sep=""))
    return(list(b.list.unique=b.list.unique,B=true.B,label=g,P=P,comm.sim.mat=true.C,node.sim.mat=node.sim.mat))
}

