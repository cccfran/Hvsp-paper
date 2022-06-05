library(Matrix)

BSHSBM <- function(n,a,b,alpha,beta,p0,p1,n.min=10,plot.it=FALSE){
    P <- Matrix(0,n,n)
    V <- 1:n
    tree.SBM <- BS(V,a,b,alpha,beta,p0,p1,n.min)
    P <- tree.SBM$P
    V <- tree.SBM$V
    tree.path <- tree.SBM$tree.path
    tree.path <- paste("All",tree.path,sep="/")
    tree.path2 <- tree.SBM$tree.path2
    tree.path2 <- paste("All",tree.path2,sep="/")
    category <- unlist(lapply(tree.path,function(x) gsub("[[:digit:]]","",x)))
    ind <- as.numeric(unlist(lapply(tree.path,function(x) grep("[[:digit:]]",unlist(strsplit(x,"/")),value=TRUE))))
    label <- rep(0,n)
    label[ind] <- as.numeric(factor(category))
    tree.df <- data.frame(V1=1:length(tree.path2),pathString=tree.path2)
    #print(tree.path)
    tree.df$pathString <- as.character(tree.df$pathString)
    cluster.tree <- as.Node(tree.df)
    if(plot.it) {
        #par(mfrow=c(1,2))
        filled.contour(z=t(P[,n:1]))
        #plot(as.dendrogram(cluster.tree),center=TRUE,leaflab="none")
        #plot(cluster.tree)
        #plot(as.dendrogram(cluster.tree))
    }
    ncl <- length(unique(label))
    block.distance.mat <- matrix(0,ncl,ncl)
    node.distance.mat <- matrix(0,nrow=n,ncol=n)
    for(k1 in 1:(ncl-1)){
        c1.index <- which(label==k1)
        for(k2 in (k1+1):ncl){
            c2.index <- which(label==k2)
            tmp.dist <- Distance(node1=FindNode(cluster.tree,as.character(c1.index[1])),node2=FindNode(cluster.tree,as.character(c2.index[1])))-2
            block.distance.mat[k1,k2] <- block.distance.mat[k2,k1] <- tmp.dist
            node.distance.mat[c1.index,c2.index] <- tmp.dist
        }
    }
    node.distance.mat <- node.distance.mat+t(node.distance.mat)
    block.distance.mat <- as.dist(block.distance.mat)
    node.distance.mat <- as.dist(node.distance.mat)
    return(list(P=P,V=V,tree.path=tree.path2,label=label,cluster.tree=cluster.tree,block.dmat=block.distance.mat,node.dmat=node.distance.mat))
}


BS <- function(V,a,b,alpha,beta,p0,p1,n.min=10){
    n <- length(V)
    P <- matrix(p1,n,n)
    prob.s <- runif(1)*(b-a)+a
    if(n<=n.min){
        print("too small!")
        prob.s <- 0
    }
    randUnif <- runif(1)
    split.flag <- FALSE
    if(randUnif<prob.s){
        split.flag <- TRUE
        if((beta*p1)<=(alpha*p0)){
            split.flag <- FALSE
            print("gap zero!")
        }
    }
    if(split.flag){
        n1 <- floor(n/2)
        p.between <- runif(1)*(beta*p1-alpha*p0)+alpha*p0
        P[1:n1,(n1+1):n] <- P[(n1+1):n,1:n1] <- p.between
        sub.BS1 <- BS(V[1:n1],a/1.5,b/1.5,alpha,beta,p.between,p1)
        sub.BS2 <- BS(V[(n1+1):n],a/1.5,b/1.5,alpha,beta,p.between,p1)
        P[1:n1,1:n1] <- sub.BS1$P
        P[(n1+1):n,(n1+1):n] <-  sub.BS2$P
        V <- list(sub.BS1$V,sub.BS2$V)
        left.path <- paste("L",sub.BS1$tree.path,sep="/")
        right.path <- paste("R",sub.BS2$tree.path,sep="/")
        tree.path <- c(left.path,right.path)
        left.path2 <- paste("L",sub.BS1$tree.path2,sep="/")
        right.path2 <- paste("R",sub.BS2$tree.path2,sep="/")
        tree.path2 <- c("",left.path2,right.path2)
    }else{
        print("Stop splitting")
        tree.path <- as.character(V)
        tree.path2 <- c("",tree.path)
        V <- list(V)

    }
    return(list(V=V,P=P,tree.path=tree.path,tree.path2=tree.path2))
}


#tt <- BSHSBM(400,1,1,1.2,0.6,0,0.5,50,TRUE)
#table(tt$label)
#filled.contour(tt$P)
#plot(as.dendrogram(tt$tree),center=TRUE,leaflab="none")







BSHSBM.RS<- function(n,a,b,beta,p1,n.min=10,plot.it=FALSE){
    P <- Matrix(0,n,n)
    V <- 1:n
    tree.SBM <- BS.RS(V,a,b,beta,p1,n.min,TRUE)
    P <- tree.SBM$P
    V <- tree.SBM$V
    tree.path <- tree.SBM$tree.path
    tree.path0 <- tree.path
    tree.path <- paste("All",tree.path,sep="/")
    tree.path2 <- tree.SBM$tree.path2
    tree.path2 <- paste("All",tree.path2,sep="/")
    category <- unlist(lapply(tree.path,function(x) gsub("[[:digit:]]","",x)))
    ind <- as.numeric(unlist(lapply(tree.path,function(x) grep("[[:digit:]]",unlist(strsplit(x,"/")),value=TRUE))))
    label <- rep(0,n)
    label[ind] <- as.numeric(factor(category))
   

    node.tree.path <- strsplit(tree.path0,"/")
    node.index <- which(unlist(lapply(node.tree.path,function(x) sum(!is.na(as.numeric(x)))>0)))
    node.number <- unlist(lapply(node.tree.path,function(x) as.numeric(x)[which(!is.na(as.numeric(x)))]))
    node.index <- node.index[sort(node.number,index.return=TRUE)$ix]
    #Binary.Similarity(node.tree.path[[node.index[1]]],node.tree.path[[node.index[1000]]])
    ncl <- length(unique(label))
    representers <- rep(0,ncl)
    for(i in 1:ncl){
        representers[i] <- which(label==i)[1]
    }
    community.bin.sim.mat <- matrix(0,ncl,ncl)
    diag(community.bin.sim.mat) <- 0
    for(i in 1:ncl){
        for(j in 1:ncl){
            if(i==j){
                community.bin.sim.mat[i,i] <- length(node.tree.path[[node.index[representers[i]]]])
            }else{
                community.bin.sim.mat[i,j] <- Binary.Similarity(node.tree.path[[node.index[representers[i]]]],node.tree.path[[node.index[representers[j]]]])
                community.bin.sim.mat[j,i] <- community.bin.sim.mat[i,j]
            }
        }
    }
    Z.mat <- matrix(0,n,ncl)
    Z.mat[cbind(1:n,label)] <- 1
    node.bin.sim.mat <- Z.mat%*%community.bin.sim.mat%*%t(Z.mat)
    diag(node.bin.sim.mat) <- 0

    
    tree.df <- data.frame(V1=1:length(tree.path2),pathString=tree.path2)
    #print(tree.path)
    tree.df$pathString <- as.character(tree.df$pathString)
    cluster.tree <- as.Node(tree.df)
    if(plot.it) {
        #par(mfrow=c(1,2))
        filled.contour(z=t(P[,n:1]))
        #plot(as.dendrogram(cluster.tree),center=TRUE,leaflab="none")
        #plot(cluster.tree)
        #plot(as.dendrogram(cluster.tree))
    }
    ncl <- length(unique(label))
    block.distance.mat <- matrix(0,ncl,ncl)
    node.distance.mat <- matrix(0,nrow=n,ncol=n)
    for(k1 in 1:(ncl-1)){
        c1.index <- which(label==k1)
        for(k2 in (k1+1):ncl){
            c2.index <- which(label==k2)
            tmp.dist <- Distance(node1=FindNode(cluster.tree,as.character(c1.index[1])),node2=FindNode(cluster.tree,as.character(c2.index[1])))-2
            block.distance.mat[k1,k2] <- block.distance.mat[k2,k1] <- tmp.dist
            node.distance.mat[c1.index,c2.index] <- tmp.dist
        }
    }
    node.distance.mat <- node.distance.mat+t(node.distance.mat)
    block.distance.mat <- as.dist(block.distance.mat)
    node.distance.mat <- as.dist(node.distance.mat)
    return(list(P=P,V=V,tree.path=tree.path2,label=label,cluster.tree=cluster.tree,block.dmat=block.distance.mat,node.dmat=node.distance.mat,node.bin.sim.mat=node.bin.sim.mat))
}


BS.RS <- function(V,a,b,beta,p1,n.min=10,split.flag=TRUE){
    n <- length(V)
    P <- matrix(p1,n,n)
    prob.s <- (a+b)/2
    if(n<=n.min){
        print("too small!")
        prob.s <- 0
        split.flag <- FALSE
    }
    randUnif <- runif(1)
    next.split.flag <- FALSE
    if(randUnif<prob.s){
        next.split.flag <- TRUE
   }
    if(split.flag){
        n1 <- floor(n/2)
        p.between <- beta*p1*0.8
        P[1:n1,(n1+1):n] <- P[(n1+1):n,1:n1] <- p.between
        sub.BS1 <- BS.RS(V[1:n1],a/1.5,b/1.5,sqrt(beta),p1,n.min,next.split.flag)
        sub.BS2 <- BS.RS(V[(n1+1):n],a/1.5,b/1.5,sqrt(beta),p1,n.min,next.split.flag)
        P[1:n1,1:n1] <- sub.BS1$P
        P[(n1+1):n,(n1+1):n] <-  sub.BS2$P
        V <- list(sub.BS1$V,sub.BS2$V)
        left.path <- paste("L",sub.BS1$tree.path,sep="/")
        right.path <- paste("R",sub.BS2$tree.path,sep="/")
        tree.path <- c(left.path,right.path)
        left.path2 <- paste("L",sub.BS1$tree.path2,sep="/")
        right.path2 <- paste("R",sub.BS2$tree.path2,sep="/")
        #tree.path <- c("",left.path,left.path2)
        tree.path2 <- c("",left.path2,right.path2)
    }else{
        print("Stop splitting")
        tree.path <- as.character(V)
        tree.path2 <- c("",tree.path)
        V <- list(V)

    }
    return(list(V=V,P=P,tree.path=tree.path,tree.path2=tree.path2))
}







BSHSBM.Balance<- function(n,a,b,beta,p1,n.min=10,plot.it=FALSE,layer.n=3){
    P <- Matrix(0,n,n)
    V <- 1:n
    tree.SBM <- BS.Balance(V,a,b,beta,p1,n.min,TRUE,layer.n)
    P <- tree.SBM$P
    V <- tree.SBM$V
    tree.path <- tree.SBM$tree.path
    tree.path <- paste("All",tree.path,sep="/")
    tree.path2 <- tree.SBM$tree.path2
    tree.path2 <- paste("All",tree.path2,sep="/")
    category <- unlist(lapply(tree.path,function(x) gsub("[[:digit:]]","",x)))
    ind <- as.numeric(unlist(lapply(tree.path,function(x) grep("[[:digit:]]",unlist(strsplit(x,"/")),value=TRUE))))
    label <- rep(0,n)
    label[ind] <- as.numeric(factor(category))
    tree.df <- data.frame(V1=1:length(tree.path2),pathString=tree.path2)
    #print(tree.path)
    tree.df$pathString <- as.character(tree.df$pathString)
    cluster.tree <- as.Node(tree.df)
    if(plot.it) {
        #par(mfrow=c(1,2))
        filled.contour(z=t(P[,n:1]))
        #plot(as.dendrogram(cluster.tree),center=TRUE,leaflab="none")
        #plot(cluster.tree)
        #plot(as.dendrogram(cluster.tree))
    }
    ncl <- length(unique(label))
    block.distance.mat <- matrix(0,ncl,ncl)
    node.distance.mat <- matrix(0,nrow=n,ncol=n)
    for(k1 in 1:(ncl-1)){
        c1.index <- which(label==k1)
        for(k2 in (k1+1):ncl){
            c2.index <- which(label==k2)
            tmp.dist <- Distance(node1=FindNode(cluster.tree,as.character(c1.index[1])),node2=FindNode(cluster.tree,as.character(c2.index[1])))-2
            block.distance.mat[k1,k2] <- block.distance.mat[k2,k1] <- tmp.dist
            node.distance.mat[c1.index,c2.index] <- tmp.dist
        }
    }
    node.distance.mat <- node.distance.mat+t(node.distance.mat)
    block.distance.mat <- as.dist(block.distance.mat)
    node.distance.mat <- as.dist(node.distance.mat)
    return(list(P=P,V=V,tree.path=tree.path2,label=label,cluster.tree=cluster.tree,block.dmat=block.distance.mat,node.dmat=node.distance.mat))
}


BS.Balance <- function(V,a,b,beta,p1,n.min=10,split.flag=TRUE,layer.n){
    n <- length(V)
    P <- matrix(p1,n,n)
    prob.s <- (a+b)/2
    if(n<=n.min){
        print("too small!")
        prob.s <- 0
        split.flag <- FALSE
    }
    randUnif <- runif(1)
    next.split.flag <- FALSE
    if(layer.n>1){
        next.split.flag <- TRUE
   }
    if(split.flag){
        n1 <- floor(n/2)
        p.between <- beta*p1*0.8
        P[1:n1,(n1+1):n] <- P[(n1+1):n,1:n1] <- p.between
        sub.BS1 <- BS.Balance(V[1:n1],a/1.5,b/1.5,sqrt(beta),p1,n.min,next.split.flag,layer.n-1)
        sub.BS2 <- BS.Balance(V[(n1+1):n],a/1.5,b/1.5,sqrt(beta),p1,n.min,next.split.flag,layer.n-1)
        P[1:n1,1:n1] <- sub.BS1$P
        P[(n1+1):n,(n1+1):n] <-  sub.BS2$P
        V <- list(sub.BS1$V,sub.BS2$V)
        left.path <- paste("L",sub.BS1$tree.path,sep="/")
        right.path <- paste("R",sub.BS2$tree.path,sep="/")
        tree.path <- c(left.path,right.path)
        left.path2 <- paste("L",sub.BS1$tree.path2,sep="/")
        right.path2 <- paste("R",sub.BS2$tree.path2,sep="/")
        #tree.path <- c("",left.path,left.path2)
        tree.path2 <- c("",left.path2,right.path2)
    }else{
        print("Stop splitting")
        tree.path <- as.character(V)
        tree.path2 <- c("",tree.path)
        V <- list(V)

    }
    return(list(V=V,P=P,tree.path=tree.path,tree.path2=tree.path2))
}
