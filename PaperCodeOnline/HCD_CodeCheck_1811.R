## Input
## A: adjacency matrix
## a: true block label
## K: number of true blocks

## Return: node cluster labels, true labels, cross-table between node cluster and true labels, NMI
library(data.table)

library(data.tree)
source("LRT.R")
source("DistanceMetric.R")
HCD <- function(A,a,K,BH=2,adj=TRUE,plot.it=FALSE,lap=FALSE,reg=FALSE,n.min=25,D=NULL){
  f <- A
  #graph.f <- graph.adjacency(f,mode="undirected")
  #f <- get.adjacency(graph.f, sparse = T)
  n <- nrow(f)
  ncl <- 0
  xi.loc.labels <- list()
  clusters <- break.cl(f, a, xi.loc.labels, ncl, 1:n,BH,lap=lap,reg=reg,n.min,D=D)
  ncl <- clusters$ncl
  xi.loc.labels <- clusters$xi.loc.labels
  labels = rep(0, n)
  for(i in 1:ncl) {
    labels[xi.loc.labels[[i]]] <- i
  }
  cross.table <- table(labels, a)
  tree.path <- clusters$tree.path
  node.tree.path <- strsplit(tree.path,"/")
  node.index <- which(unlist(lapply(node.tree.path,function(x) sum(!is.na(as.numeric(x)))>0)))
  node.number <- unlist(lapply(node.tree.path,function(x) as.numeric(x)[which(!is.na(as.numeric(x)))]))
  node.dt <- data.table(node.number=node.number,node.index=node.index)
  node.dt2 <- unique(node.dt, by = "node.number")
  node.index <- node.dt2$node.index[sort(node.dt2$node.number,index.return=TRUE)$ix]
  #Binary.Similarity(node.tree.path[[node.index[1]]],node.tree.path[[node.index[1000]]])
  representers <- rep(0,ncl)
  for(i in 1:ncl){
    representers[i] <- which(labels==i)[1]
  }
  community.bin.sim.mat <- matrix(0,ncl,ncl)
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
  Z.mat[cbind(1:n,labels)] <- 1
  node.bin.sim.mat <- Z.mat%*%community.bin.sim.mat%*%t(Z.mat)
  diag(node.bin.sim.mat) <- 0
  
  tree.path <- paste("All",tree.path,sep="/")
  tree.df <- data.frame(V1=1:length(tree.path),pathString=tree.path)
  #print(tree.path)
  tree.df$pathString <- as.character(tree.df$pathString)
  cluster.tree <- as.Node(tree.df)
  block.distance.mat <- matrix(0,ncl,ncl)
  node.distance.mat <- matrix(0,nrow=n,ncol=n)
  for(k1 in 1:(ncl-1)){
    c1.index <- which(labels==k1)
    if(ncl>1){
      for(k2 in (k1+1):ncl){
        c2.index <- which(labels==k2)
        #print(c(k1,k2))
        tmp.dist <- Distance(node1=FindNode(cluster.tree,as.character(c1.index[1])),node2=FindNode(cluster.tree,as.character(c2.index[1])))-2
        block.distance.mat[k1,k2] <- block.distance.mat[k2,k1] <- tmp.dist
        node.distance.mat[c1.index,c2.index] <- tmp.dist
      }
    }
  }
  node.distance.mat <- node.distance.mat+t(node.distance.mat)
  block.distance.mat <- as.dist(block.distance.mat)
  node.distance.mat <- as.dist(node.distance.mat)#as(node.distance.mat,"dgCMatrix")
  NMI <- mi.empirical(as.matrix(cross.table))/entropy.empirical(as.vector(cross.table))
  P.est.node.dmat <- P.est <- NULL
  if(adj){
    P.est <- SBM.estimate(A,labels)
    HCD.P.est <- HCD.pop(A=P.est,a=labels,K=ncl)
    P.est.node.dmat <- HCD.P.est$node.dmat
    P.node.bin.sim.mat <- HCD.P.est$node.bin.sim.mat
  }
  result <- list(labels=labels,true.label=a,cross.table=cross.table,NMI=NMI,ncl=ncl,cluster.tree=cluster.tree,mod.path=clusters$mod.path,block.dmat=block.distance.mat,node.dmat=node.distance.mat,P=P.est,P.node.dmat=P.est.node.dmat,ncl.P=HCD.P.est$ncl,node.bin.sim.mat=node.bin.sim.mat,P.node.bin.sim.mat=P.node.bin.sim.mat,comm.bin.sim.mat=community.bin.sim.mat,tree.path=tree.path)
  
  if(plot.it) plot(as.dendrogram(cluster.tree,heightAttribute="V1"),center=T,leaflab="none")
  return(result)
}


## cl.labels is the current node index involved in the function call
break.cl = function(f, per.a, xi.loc.labels, ncl, cl.labels,BH=2,lap=FALSE,reg=FALSE,n.min=25,D=NULL) {
  nisol = which(rowSums(f) > 0)
  isol = which(rowSums(f) == 0)
  if((length(nisol)<=8)||(length(isol)>=5*length(nisol))||(length(nisol)<2*n.min)){
    ncl = ncl + 1
    xi.loc.labels[[ncl]] = cl.labels
    #if(length(isol)>0) xi.loc.labels[[ncl]] = c(cl.labels,cl.labels) ## ????? need to check
    tree.path <- c("",as.character(cl.labels))
    mod.path <- c(0,rep(0,length(cl.labels)))
    if(length(isol)>0) tree.path <- c(tree.path,as.character(cl.labels[isol]))
    print('Too few connected nodes, not even started!')
    return(list(xi.loc.labels = xi.loc.labels, ncl = ncl,tree.path=tree.path,mod.path=mod.path))
  }
  print(paste(length(isol),"isolated nodes",cl.labels[isol]))
  all.deg = rowSums(f)
  f = f[nisol, nisol] ### only focus on non-isolated nodes - how to label isolated nodes?
  cl.labels.full <- cl.labels
  per.a = per.a[nisol] ### true label
  all.deg = all.deg[nisol]
  cl.labels = cl.labels[nisol]
  n1 = dim(f)[1]
  phat = sum(f)/(n1*(n1-1))
  pnull <- phat
  sum.f <- sum(f)
  norm.f = (f - phat)/(sqrt((n1-1)*phat*(1-phat)))
  
  K = 2
  if(lap){
    if(reg){
      f.reg <- f + 0.1*mean(all.deg)/nrow(f) ## use regularized laplacian for spectral clustering
      all.deg.reg <- colSums(f.reg)
      L <- t(t(f.reg/sqrt(all.deg.reg))/sqrt(all.deg.reg))
    }else{
      L <- t(t(f/sqrt(all.deg))/sqrt(all.deg))
      
    }
    svd.L <- irlba(L,nu=K,nv=K)
  }
  split.flag <- FALSE
  print(dim(norm.f))
  #eig.nf = irlba(norm.f, nu = K, nv = K)   ## old code from Sharmo, keep it for stability, should not be useful if we use laplacian. will remove
  if(BH<=2){
    if(BH==2){
      split.flag <- NB.check(f)
    }
    if(BH==1){
      split.flag <- ECV.check(f)
    }
  }else{
    if(BH==3){
      lrt <- WangLRT(f,3)
      if(lrt$SBM.K>2) split.flag <- TRUE
    }else{
      if(is.null(D)){
        eig.nf1 <- irlba(f, nu = K, nv = K)
        #    plot(eig.nf$u[,1], ylab = 'Eigenvector',col=per.a)
        bnd <- 2.05*sqrt(n1*(phat*0.25 + phat*(1-phat)*0.75))
        n.imp.eig <- length(which(eig.nf1$d >= bnd))
        split.flag <- (n.imp.eig>1)}else{
          split.flag <- (D>0)
        }
    }
  }
  if(split.flag) {
    if(lap){
      clus = pam(svd.L$v[,1:2], 2, metric = 'euclidean')
    }else{
      eig.nf = irlba(norm.f, nu = K, nv = K)
      clus = pam(eig.nf$v[,1], 2, metric = 'euclidean')
    }
    xi.f = clus$clustering
    xi.labels = lapply(1:2, function(x){which(xi.f == x)})
    #print(length(xi.f))
    #print(length(per.a))
    print(as.matrix(table(xi.f, per.a)))
    smaller.cluster <- xi.labels[[which.min(sapply(xi.labels,length))]]
    f1 <- f[smaller.cluster,smaller.cluster]
    sum.f1 <- sum(f1-pnull)
    num.f1 <- nrow(f1)
    #pf1 <- sum(f1)/(num.f1^2-num.f1)
    #print(paste("pf1 =",pf1))
    #print(dim(f1))
    a1 <- per.a[smaller.cluster]
    a1.labels <- cl.labels[smaller.cluster]
    if(length(dim(f1)) > 0) {
      #nisol = which(rowSums(f1) > 0)
      f1n = f1
      #print(dim(f1n))
      if(length(dim(f1n)) > 0) {
        n1 = dim(f1n)[1]
        phat = sum(f1n)/(n1*(n1-1))
      } else {
        n1 = 0
      }
    } else if(length(f1) > 0) { ### case when only 1 node is in this cluster, make it 2, so the later on rank check code still works
      f1n = diag(rep(f1, 2))
      n1 = dim(f1n)[1]
      phat = sum(f1n)/(n1*(n1-1))
    } else {
      n1 = 0
    }
    if(n1 > 2*n.min) { ## only do further clustering on cluster larger than 50, should change to 2*n.min
      #print(min(nrow(f1n),ncol(f1n)))
      if((phat != 1) && (phat != 0) && (sum(irlba(f1n,nu=2,nv=2)$d>1e-6)==2) && (nrow(f1n) > 4)) {
        #print("keep going!")
        res = break.cl(f1, a1, xi.loc.labels, ncl, a1.labels,BH,lap=lap,reg=reg,n.min,D=D-1)
        xi.loc.labels = res$xi.loc.labels
        ncl = res$ncl
        if(length(isol)>0) xi.loc.labels[[ncl]] = c(xi.loc.labels[[ncl]],cl.labels.full[isol]) ### attached the isolated nodes in this level with the clusters under the smaller split
        L.tree.path <- res$tree.path
        if(length(isol)>0){
          path.head <- L.tree.path[length(L.tree.path)]
          path.head <- gsub('[[:digit:]]+', '', path.head)
          iso.path <- paste(path.head,cl.labels.full[isol],sep="")
          L.tree.path <- c(L.tree.path,iso.path)
        }
        L.mod.path <- res$mod.path
      } else {
        ncl = ncl + 1
        xi.loc.labels[[ncl]] = a1.labels
        if(length(isol)>0) xi.loc.labels[[ncl]] = c(cl.labels.full[isol],xi.loc.labels[[ncl]])
        L.tree.path <- as.character(xi.loc.labels[[ncl]])
        L.mod.path <- rep(0,length(xi.loc.labels[[ncl]]))
        print('Homogeneous, Branch End')
      }
    } else {
      ncl = ncl + 1
      xi.loc.labels[[ncl]] = a1.labels
      if(length(isol)>0) xi.loc.labels[[ncl]] = c(cl.labels.full[isol],xi.loc.labels[[ncl]])
      L.tree.path <- as.character(xi.loc.labels[[ncl]])
      L.mod.path <- rep(0,length(xi.loc.labels[[ncl]]))
      print('Too small cluster, Branch End')
    }
    f2 = f[-smaller.cluster, -smaller.cluster]
    sum.f2 <- sum(f2-pnull)
    num.f2 <- nrow(f2)
    #pf2 <- sum.f2/(num.f2^2-num.f2)
    #print(paste("pf2 =",pf2))
    
    a2 = per.a[-smaller.cluster]
    a2.labels = cl.labels[-smaller.cluster]
    if(length(dim(f2)) > 0) {
      #nisol = which(rowSums(f2) > 0)
      f2n = f2#[nisol, nisol]
      #print(dim(f2n))
      if(length(dim(f2n)) > 0) {
        n1 <- nrow(f2n)
        phat = sum(f2n)/(n1*(n1-1))
      } else {
        n1 <- 0
      }
    } else if(length(f2) > 0) {
      f2n = diag(rep(f2, 2))
      n1 <- nrow(f2n)
      phat = sum(f2n)/(n1*(n1-1))
    } else {
      n1 <- 0
    }
    if(n1 > 2*n.min) {  ## change to 2*n.min
      #print(dim(f2n))
      if((phat != 1) && (phat != 0) && (sum(irlba(f2n,nu=2,nv=2)$d>1e-6)==2) && (nrow(f2n) > 4)) {
        res = break.cl(f2, a2, xi.loc.labels, ncl, a2.labels,BH,lap=lap,reg=reg,n.min,D=D-1)
        xi.loc.labels = res$xi.loc.labels
        R.tree.path <- res$tree.path
        R.mod.path <- res$mod.path
        ncl = res$ncl
      } else {
        ncl = ncl + 1
        xi.loc.labels[[ncl]] = a2.labels
        R.tree.path <- as.character(a2.labels)
        R.mod.path <- rep(0,length(a2.labels))
        print('Homogenous, Branch End')
      }
    } else {
      ncl = ncl + 1
      xi.loc.labels[[ncl]] = a2.labels
      R.tree.path <- as.character(a2.labels)
      R.mod.path <- rep(0,length(a2.labels))
      print('Too small cluster, Branch End')
    }
    L.tree.path <- paste("L",L.tree.path,sep="/")
    #if(length(isol)>0) R.tree.path <- c(R.tree.path,as.character(cl.labels.full[isol]))
    R.tree.path <- paste("R",R.tree.path,sep="/")
    tree.path <- c("",L.tree.path,R.tree.path)
    
    modularity <- (sum.f1+sum.f2)/sum(f)  ## check modularity, not useful, to be removed
    mod.path <- c(modularity,L.mod.path,R.mod.path)
  } else {
    ncl = ncl + 1
    xi.loc.labels[[ncl]] = cl.labels
    if(length(isol)>0) xi.loc.labels[[ncl]] = c(cl.labels.full[isol],cl.labels)
    tree.path <- c("",as.character(cl.labels))
    mod.path <- c(0,rep(0,length(cl.labels)))
    if(length(isol)>0) tree.path <- c(tree.path,as.character(cl.labels.full[isol]))
    print('One cluster, Branch End, not even started!')
  }
  return(list(xi.loc.labels = xi.loc.labels, ncl = ncl,tree.path=tree.path,mod.path=mod.path))
}





SpectralClust <- function(A,a,K,lap=FALSE,reg=FALSE){
  graph.f <- graph.adjacency(A,mode="undirected")
  A <- get.adjacency(graph.f, sparse = T)
  f <- A
  all.deg <- colSums(f)
  nisol <- which(all.deg>0)
  isol <- which(all.deg==0)
  f <- f[nisol,nisol]
  n1 = dim(f)[1]
  all.deg <- all.deg[nisol]
  phat = sum(f)/(n1*(n1-1))
  norm.f = (f - phat)/(sqrt((n1-1)*phat*(1-phat)))
  #d <- colSums(A)
  #D.inv <- diag(1/sqrt(d))
  #L <- D.inv%*%A%*%D.inv
  #SVD <- irlba(A,nu=K,nv=K)
  if(lap){
    if(reg){
      f.reg <- f + 0.1*mean(all.deg)/n1
      all.deg.reg <- colSums(f.reg)
      L <- t(t(f.reg/sqrt(all.deg.reg))/sqrt(all.deg.reg))
    }else{
      L <- t(t(f/sqrt(all.deg))/sqrt(all.deg))
    }
    SVD <- irlba(L,nu=K,nv=K)
  }else{
    SVD <- irlba(norm.f,nu=K,nv=K)
  }
  if(lap){
    clus <- pam(SVD$v,k=K)
  }else{
    if(K>1){
      clus <- pam(SVD$v[,1:(K-1)],k=K)
    }else{
      clus <- pam(SVD$v,k=K)
    }
  }
  labels <- rep(0,nrow(A))
  labels[nisol] <- clus$clustering
  min.label <- which.min(table(clus$clustering))
  labels[isol] <- min.label
  cross.table <- table(labels,a)
  NMI <- mi.empirical(as.matrix(cross.table))/entropy.empirical(as.vector(cross.table))
  P.est.node.dmat <- P.est <- NULL
  P.est <- SBM.estimate(A,labels)
  HCD.P.est <- HCD.pop(A=P.est,a=labels,K=K)
  P.est.node.dmat <- HCD.P.est$node.dmat
  P.node.bin.sim.mat <- HCD.P.est$node.bin.sim.mat
  result <- list(labels=labels,true.label=a,cross.table=cross.table,NMI=NMI,P=P.est,P.node.dmat=P.est.node.dmat,ncl.P=HCD.P.est$ncl,P.node.bin.sim.mat=P.node.bin.sim.mat)
  
}




### HCD procedure on population probability matrix: the only difference from HCD is that the stopping rule is to check whether all elements of the population matrix are the same.
HCD.pop <- function(A,a,K,BH=0,D=NULL){
  print("call population version")
  f <- A
  #graph.f <- graph.adjacency(f,mode="undirected")
  #f <- get.adjacency(graph.f, sparse = T)
  n <- nrow(f)
  
  ncl <- 0
  xi.loc.labels <- list()
  clusters <- break.cl.pop(f, a, xi.loc.labels, ncl, 1:n,BH,D=D)
  ncl <- clusters$ncl
  xi.loc.labels <- clusters$xi.loc.labels
  labels = rep(0, n)
  for(i in 1:ncl) {
    labels[xi.loc.labels[[i]]] <- i
  }
  cross.table <- table(labels, a)
  tree.path <- clusters$tree.path
  node.tree.path <- strsplit(tree.path,"/")
  node.index <- which(unlist(lapply(node.tree.path,function(x) sum(!is.na(as.numeric(x)))>0)))
  node.number <- unlist(lapply(node.tree.path,function(x) as.numeric(x)[which(!is.na(as.numeric(x)))]))
  node.dt <- data.table(node.number=node.number,node.index=node.index)
  node.dt2 <- unique(node.dt, by = "node.number")
  node.index <- node.dt2$node.index[sort(node.dt2$node.number,index.return=TRUE)$ix]
  #node.index <- node.index[sort(node.number,index.return=TRUE)$ix]
  #Binary.Similarity(node.tree.path[[node.index[1]]],node.tree.path[[node.index[1000]]])
  representers <- rep(0,ncl)
  for(i in 1:ncl){
    representers[i] <- which(labels==i)[1]
  }
  community.bin.sim.mat <- matrix(0,ncl,ncl)
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
  Z.mat[cbind(1:n,labels)] <- 1
  node.bin.sim.mat <- Z.mat%*%community.bin.sim.mat%*%t(Z.mat)
  diag(node.bin.sim.mat) <- 0
  tree.path <- paste("All",tree.path,sep="/")
  tree.df <- data.frame(V1=clusters$mod.path,pathString=tree.path)
  #print(tree.path)
  tree.df$pathString <- as.character(tree.df$pathString)
  cluster.tree <- as.Node(tree.df)
  block.distance.mat <- matrix(0,ncl,ncl)
  node.distance.mat <- matrix(0,nrow=n,ncol=n)
  for(k1 in 1:(ncl-1)){
    c1.index <- which(labels==k1)
    if(ncl>1){
      for(k2 in (k1+1):ncl){
        c2.index <- which(labels==k2)
        tmp.dist <- Distance(node1=FindNode(cluster.tree,as.character(c1.index[1])),node2=FindNode(cluster.tree,as.character(c2.index[1])))-2
        block.distance.mat[k1,k2] <- block.distance.mat[k2,k1] <- tmp.dist
        node.distance.mat[c1.index,c2.index] <- tmp.dist
      }
    }
  }
  node.distance.mat <- node.distance.mat+t(node.distance.mat)
  block.distance.mat <- as.dist(block.distance.mat)
  node.distance.mat <- as.dist(node.distance.mat)#as(node.distance.mat,"dgCMatrix")
  NMI <- mi.empirical(as.matrix(cross.table))/entropy.empirical(as.vector(cross.table))
  
  result <- list(labels=labels,true.label=a,cross.table=cross.table,NMI=NMI,ncl=ncl,cluster.tree=cluster.tree,mod.path=clusters$mod.path,block.dmat=block.distance.mat,node.dmat=node.distance.mat,node.bin.sim.mat=node.bin.sim.mat)
  
  return(result)
}





break.cl.pop = function(f, per.a, xi.loc.labels, ncl, cl.labels,BH=0,D=NULL) {
  nisol = which(rowSums(f) > 0)
  isol = which(rowSums(f) == 0)
  all.deg = rowSums(f)
  f = f[nisol, nisol]
  cl.labels.full <- cl.labels
  per.a = per.a[nisol]
  all.deg = all.deg[nisol]
  cl.labels = cl.labels[nisol]
  n1 = dim(f)[1]
  phat = sum(f)/(n1*(n1-1))
  pnull <- phat
  sum.f <- sum(f)
  norm.f = (f - phat)/(sqrt((n1-1)*phat*(1-phat)))
  
  
  K = 2
  split.flag <- FALSE
  #L <- t(t(f*sqrt(1/all.deg))*(sqrt(1/all.deg)))
  #svd.L = irlba(L, nu = K, nv = K)
  #svd.L = irlba(norm.f, nu = K, nv = K)
  eig.nf = irlba(norm.f, nu = K, nv = K)
  f.reg <- f + 0.1*mean(all.deg)/n1
  all.deg.reg <- colSums(f.reg)
  L <- t(t(f.reg/sqrt(all.deg.reg))/sqrt(all.deg.reg))
  svd.L <- irlba(L,nu=K,nv=K)
  if(BH<=2){
    if(BH==0){
      split.flag <- (sum(unique(f))>1)
      #print(paste("BH=0 check",sum(unique(f))))
    }
    if(BH==2){
      split.flag <- NB.check(f)
    }
    if(BH==1){
      split.flag <- ECV.check(f)
    }
  }else{
    if(BH==3){
      lrt <- WangLRT(f,3)
      if(lrt$SBM.K>2) split.flag <- TRUE
    }else{
      if(is.null(D)){
        eig.nf1 <- irlba(f, nu = K, nv = K)
        #    plot(eig.nf$u[,1], ylab = 'Eigenvector',col=per.a)
        bnd <- 2.05*sqrt(n1*(phat*0.25 + phat*(1-phat)*0.75))
        n.imp.eig <- length(which(eig.nf1$d >= bnd))
        split.flag <- (n.imp.eig>1)}else{
          split.flag <- (D>0)
        }
    }
  }
  #print(paste("Split flag is ",split.flag))
  if(split.flag) {
    #clus = pam(eig.nf$v[,1], 2, metric = 'euclidean')
    clus = pam(svd.L$v[,1:2], 2, metric = 'euclidean')
    xi.f = clus$clustering
    xi.labels = lapply(1:2, function(x){which(xi.f == x)})
    #print(length(xi.f))
    #print(length(per.a))
    #print(as.matrix(table(xi.f, per.a)))
    smaller.cluster <- xi.labels[[which.min(sapply(xi.labels,length))]]
    f1 <- f[smaller.cluster,smaller.cluster]
    sum.f1 <- sum(f1-pnull)
    num.f1 <- nrow(f1)
    #pf1 <- sum(f1)/(num.f1^2-num.f1)
    #print(paste("pf1 =",pf1))
    #print(dim(f1))
    a1 <- per.a[smaller.cluster]
    a1.labels <- cl.labels[smaller.cluster]
    if(length(dim(f1)) > 0) {
      nisol = which(rowSums(f1) > 0)
      f1n = f1[nisol, nisol]
      #print(dim(f1n))
      if(length(dim(f1n)) > 0) {
        n1 = dim(f1n)[1]
        phat = sum(f1n)/(n1*(n1-1))
      } else {
        n1 = 0
      }
    } else if(length(f1) > 0) {
      f1n = diag(rep(f1, 2))
      n1 = dim(f1n)[1]
      phat = sum(f1n)/(n1*(n1-1))
    } else {
      n1 = 0
    }
    if(n1 > 8) {
      #print(min(nrow(f1n),ncol(f1n)))
      if((phat != 1) && (phat != 0) && (sum(irlba(f1n,nu=2,nv=2)$d>1e-6)==2) && (nrow(f1n) > 4)) {
        #print("keep going!")
        res = break.cl.pop(f1, a1, xi.loc.labels, ncl, a1.labels,BH,D=D-1)
        xi.loc.labels = res$xi.loc.labels
        ncl = res$ncl
        if(length(isol)>0) xi.loc.labels[[ncl]] = c(cl.labels.full[isol],xi.loc.labels[[ncl]])
        L.tree.path <- res$tree.path
        L.mod.path <- res$mod.path
      } else {
        ncl = ncl + 1
        xi.loc.labels[[ncl]] = a1.labels
        if(length(isol)>0) xi.loc.labels[[ncl]] = c(cl.labels.full[isol],xi.loc.labels[[ncl]])
        L.tree.path <- as.character(a1.labels)
        L.mod.path <- rep(0,length(a1.labels))
        print('Left Split: Too small cluster, Branch End, n>8')
      }
    } else {
      ncl = ncl + 1
      xi.loc.labels[[ncl]] = a1.labels
      if(length(isol)>0) xi.loc.labels[[ncl]] = c(cl.labels.full[isol],xi.loc.labels[[ncl]])
      L.tree.path <- as.character(a1.labels)
      L.mod.path <- rep(0,length(a1.labels))
      print('Left Split: Too small cluster, Branch End, n<=8')
    }
    f2 = f[-smaller.cluster, -smaller.cluster]
    sum.f2 <- sum(f2-pnull)
    num.f2 <- nrow(f2)
    #pf2 <- sum.f2/(num.f2^2-num.f2)
    #print(paste("pf2 =",pf2))
    
    a2 = per.a[-smaller.cluster]
    a2.labels = cl.labels[-smaller.cluster]
    if(length(dim(f2)) > 0) {
      nisol = which(rowSums(f2) > 0)
      f2n = f2[nisol, nisol]
      #print(dim(f2n))
      if(length(dim(f2n)) > 0) {
        n1 <- nrow(f2n)
        phat = sum(f2n)/(n1*(n1-1))
      } else {
        n1 <- 0
      }
    } else if(length(f2) > 0) {
      f2n = diag(rep(f2, 2))
      n1 <- nrow(f2n)
      phat = sum(f2n)/(n1*(n1-1))
    } else {
      n1 <- 0
    }
    if(n1 > 8) {
      #print(dim(f2n))
      if((phat != 1) && (phat != 0) && (sum(irlba(f2n,nu=2,nv=2)$d>1e-6)==2) && (nrow(f2n) > 4)) {
        res = break.cl.pop(f2, a2, xi.loc.labels, ncl, a2.labels,BH,D=D-1)
        xi.loc.labels = res$xi.loc.labels
        R.tree.path <- res$tree.path
        R.mod.path <- res$mod.path
        ncl = res$ncl
      } else {
        ncl = ncl + 1
        xi.loc.labels[[ncl]] = a2.labels
        R.tree.path <- as.character(a2.labels)
        R.mod.path <- rep(0,length(a2.labels))
        print('Right split: Too small cluster, Branch End, n>8')
      }
    } else {
      ncl = ncl + 1
      xi.loc.labels[[ncl]] = a2.labels
      R.tree.path <- as.character(a2.labels)
      R.mod.path <- rep(0,length(a2.labels))
      print('Right split: Too small cluster, Branch End, n<=8')
    }
    L.tree.path <- paste("L",L.tree.path,sep="/")
    if(length(isol)>0) R.tree.path <- c(R.tree.path,as.character(cl.labels.full[isol]))
    R.tree.path <- paste("R",R.tree.path,sep="/")
    tree.path <- c("",L.tree.path,R.tree.path)
    
    modularity <- (sum.f1+sum.f2)/sum(f)
    mod.path <- c(modularity,L.mod.path,R.mod.path)
  } else {
    ncl = ncl + 1
    xi.loc.labels[[ncl]] = cl.labels
    if(length(isol)>0) xi.loc.labels[[ncl]] = c(cl.labels.full[isol],cl.labels)
    tree.path <- c("",as.character(cl.labels))
    mod.path <- c(0,rep(0,length(cl.labels)))
    if(length(isol)>0) tree.path <- c(tree.path,as.character(cl.labels.full[isol]))
    print('One cluster, Branch End, not even started!')
    #print(split.flag)
  }
  return(list(xi.loc.labels = xi.loc.labels, ncl = ncl,tree.path=tree.path,mod.path=mod.path))
}




## Use the non-back tracking method to check stoppping rule
library(RSpectra)
BHMC.check <- function(A){
  d <- colSums(A)
  n <- nrow(A)
  I <- as(diag(rep(1,n)),"dgCMatrix")
  D <- as(diag(d),"dgCMatrix")
  r <- sqrt(sum(d^2)/sum(d)-1)
  B <- as(matrix(0,nrow=2*n,ncol=2*n),"dgCMatrix")
  B[(n+1):(2*n),1:n] <- -1*I
  B[1:n,(n+1):(2*n)] <- D-I
  B[(n+1):(2*n),(n+1):(2*n)] <- A
  #r <- sqrt(sum(d)/n)
  #r <- sqrt(irlba(B,nu=1,nv=1)$d[1])
  BH <- (r^2-1)*I-r*A+D
  rho <- eigs_sym(BH,15,which="LM",sigma=0)$values
  diff <- rho[1:13]-5*rho[2:14]
  #if(rho[1]>5*rho[2]) return(TRUE)
  return(sum(diff>0)>0)
}


BHMC.estimate <- function(A){
  d <- colSums(A)
  n <- nrow(A)
  I <- as(diag(rep(1,n)),"dgCMatrix")
  D <- as(diag(d),"dgCMatrix")
  r <- sqrt(sum(d^2)/sum(d)-1)
  B <- as(matrix(0,nrow=2*n,ncol=2*n),"dgCMatrix")
  B[(n+1):(2*n),1:n] <- -1*I
  B[1:n,(n+1):(2*n)] <- D-I
  B[(n+1):(2*n),(n+1):(2*n)] <- A
  #r <- sqrt(sum(d)/n)
  #r <- sqrt(irlba(B,nu=1,nv=1)$d[1])
  BH <- (r^2-1)*I-r*A+D
  rho <- eigs_sym(BH,15,which="LM",sigma=0)$values
  diff <- rho[1:14]-5*rho[2:15]
  #if(rho[1]>5*rho[2]) return(TRUE)
  return(15-min(which(diff>0)))
}



NB.check <- function(A){
  d <- colSums(A)
  n <- nrow(A)
  I <- as(diag(rep(1,n)),"dgCMatrix")
  D <- as(diag(d),"dgCMatrix")
  r <- sqrt(sum(d^2)/sum(d)-1)
  B <- as(matrix(0,nrow=2*n,ncol=2*n),"dgCMatrix")
  B[(n+1):(2*n),1:n] <- -1*I
  B[1:n,(n+1):(2*n)] <- D-I
  B[(n+1):(2*n),(n+1):(2*n)] <- A
  ss <- Re(eigs(B,k=2,which="LM")$values)
  #r <- sqrt(irlba(B,nu=1,nv=1)$d[1])
  return(sum(abs(ss)>r)>1)
}



NB.estimate <- function(A,K=NULL){
  d <- colSums(A)
  n <- nrow(A)
  I <- as(diag(rep(1,n)),"dgCMatrix")
  D <- as(diag(d),"dgCMatrix")
  r <- sqrt(sum(d^2)/sum(d)-1)
  B <- as(matrix(0,nrow=2*n,ncol=2*n),"dgCMatrix")
  B[(n+1):(2*n),1:n] <- -1*I
  B[1:n,(n+1):(2*n)] <- D-I
  B[(n+1):(2*n),(n+1):(2*n)] <- A
  if(is.null(K)){
    K <- ceiling(sqrt(n))
  }
  ss <- Re(eigs(B,k=K,which="LM")$values)
  #r <- sqrt(irlba(B,nu=1,nv=1)$d[1])
  return(sum(abs(ss)>r))
}

partition_leaves <- function(dend, ...) {
  if (!is.dendrogram(dend)) stop("'dend' is not a dendrogram")
  
  nodes_labels <- vector("list", length = nnodes(dend))
  
  i_counter <- 0
  push_node_labels <- function(dend_node) {
    i_counter <<- i_counter + 1
    
    nodes_labels[[i_counter]] <<- labels(dend_node)
    return(NULL)
  }
  dendrapply(dend, push_node_labels)
  
  return(nodes_labels)
}


count.edge.dend <- function(dend){
  sort_a_character <- function(dend) dend %>% as.character %>% sort
  
  bp1 <- partition_leaves(dend)
  bp1 <- lapply(bp1, sort_a_character)
  return(length(bp1))
  
}




SBM.estimate <- function(A,g){
  K <- length(unique(g))
  n <- nrow(A)
  Z <- matrix(0,n,K)
  Z[cbind(1:n,g)] <- 1
  B <- matrix(0,K,K)
  for(i in 1:K){
    for(j in i:K){
      if(i==j){
        i.index <- which(g==i)
        i.num <- length(i.index)
        B[i,i] <- sum(A[i.index,i.index])/(i.num^2-i.num)
      }else{
        i.index <- which(g==i)
        j.index <- which(g==j)
        B[i,j] <- B[j,i] <- mean(A[i.index,j.index])
      }
    }
  }
  P <- Z%*%B%*%t(Z)
  return(P)
}

source("RandomHoldout.R")
ECV.check <- function(A,B=3){
  ecv <- ECV.undirected.Rank(A,max.K=2,B=B,holdout.p=0.1,weighted=TRUE)
  if(ecv$sse[2] < ecv$sse[1]){
    return(TRUE)
  }else{
    return(FALSE)
  }
}






## cl.labels is the current node index involved in the function call - this is by sign splitting - use D = maximum level of splitting if want to use HCD-SP by fixed D with BH=4
break.cl.sp = function(f, per.a, xi.loc.labels, ncl, cl.labels,BH=2,lap=FALSE,reg=FALSE,n.min=25,D=NULL) {
  nisol = which(rowSums(f) > 0)
  isol = which(rowSums(f) == 0)
  cl.labels.full <- cl.labels
  if((length(nisol)<=8)||(length(isol)>=5*length(nisol))||(length(nisol)<2*n.min)){
    ncl = ncl + 1
    xi.loc.labels[[ncl]] = cl.labels
    tree.path <- c("",as.character(cl.labels))
    mod.path <- c(0,rep(0,length(cl.labels)))
    if(length(isol)>0) tree.path <- c(tree.path,as.character(cl.labels[isol]))
    print('Too few connected nodes, not even started!')
    return(list(xi.loc.labels = xi.loc.labels, ncl = ncl,tree.path=tree.path,mod.path=mod.path))
  }
  #print(paste(length(isol),"isolated nodes",cl.labels[isol]))
  all.deg = rowSums(f)
  f = f[nisol, nisol] ### only focus on non-isolated nodes - how to label isolated nodes?
  cl.labels.full <- cl.labels
  per.a = per.a[nisol] ### true label
  all.deg = all.deg[nisol]
  cl.labels = cl.labels[nisol]
  n1 = dim(f)[1]
  phat = sum(f)/(n1*(n1-1))
  pnull <- phat
  sum.f <- sum(f)
  norm.f = (f - phat)/(sqrt((n1-1)*phat*(1-phat)))
  
  K = 2
  if(lap){
    if(reg){
      f.reg <- f + 0.1*mean(all.deg)/nrow(f) ## use regularized laplacian for spectral clustering
      all.deg.reg <- colSums(f.reg)
      L <- t(t(f.reg/sqrt(all.deg.reg))/sqrt(all.deg.reg))
    }else{
      L <- t(t(f/sqrt(all.deg))/sqrt(all.deg))
      
    }
    svd.L <- irlba(L,nu=K,nv=K)
  }
  split.flag <- FALSE
  #print(dim(norm.f))
  #eig.nf = irlba(norm.f, nu = K, nv = K)   ## old code from Sharmo, keep it for stability, should not be useful if we use laplacian. will remove
  if(BH<=2){
    if(BH==2){
      split.flag <- NB.check(f)
    }
    if(BH==1){
      split.flag <- ECV.check(f)
    }
  }else{
    if(BH==3){
      lrt <- WangLRT(f,3)
      if(lrt$SBM.K>2) split.flag <- TRUE
    }else{
      #eig.nf = irlba(norm.f, nu = K, nv = K)
      if(is.null(D)){
        eig.nf1 <- irlba(f, nu = K, nv = K)
        #    plot(eig.nf$u[,1], ylab = 'Eigenvector',col=per.a)
        bnd <- 2.05*sqrt(n1*(phat*0.25 + phat*(1-phat)*0.75))
        n.imp.eig <- length(which(eig.nf1$d >= bnd))
        split.flag <- (n.imp.eig>1)}else{
          split.flag <- (D>0)
        }
    }
  }
  if(split.flag) {
    if(lap){
      clus = pam(svd.L$v[,1:2], 2, metric = 'euclidean')
    }else{
      if(reg){
        f.reg <- f + 0.1*mean(all.deg)/nrow(f)
        eig.nf = irlba(f.reg, nu = 2, nv = 2)
      }else{
        eig.nf = irlba(f, nu = 2, nv = 2)
      }
      clustering <- rep(0,length(eig.nf$v[,2]))
      clustering[eig.nf$v[,2]<=0] <- 1
      clustering[eig.nf$v[,2]>0] <- 2
      clus = list(clustering=clustering)
    }
    xi.f = clus$clustering
    xi.labels = lapply(1:2, function(x){which(xi.f == x)})
    #print(length(xi.f))
    #print(length(per.a))
    #print(as.matrix(table(xi.f, per.a)))
    smaller.cluster <- xi.labels[[which.min(sapply(xi.labels,length))]]
    f1 <- f[smaller.cluster,smaller.cluster]
    sum.f1 <- sum(f1-pnull)
    num.f1 <- nrow(f1)
    #pf1 <- sum(f1)/(num.f1^2-num.f1)
    #print(paste("pf1 =",pf1))
    #print(dim(f1))
    a1 <- per.a[smaller.cluster]
    a1.labels <- cl.labels[smaller.cluster]
    if(length(dim(f1)) > 0) {
      #nisol = which(rowSums(f1) > 0)
      f1n = f1
      #print(dim(f1n))
      if(length(dim(f1n)) > 0) {
        n1 = dim(f1n)[1]
        phat = sum(f1n)/(n1*(n1-1))
      } else {
        n1 = 0
      }
    } else if(length(f1) > 0) { ### case when only 1 node is in this cluster, make it 2, so the later on rank check code still works
      #print("1 node community produced, this is")
      #print(smaller.cluster)
      f1n = diag(rep(f1, 2))
      n1 = dim(f1n)[1]
      phat = sum(f1n)/(n1*(n1-1))
    } else {
      n1 = 0
    }
    if(n1 > 2*n.min) { ## only do further clustering on cluster larger than 50, should change to 2*n.min
      #print(min(nrow(f1n),ncol(f1n)))
      if((phat != 1) && (phat != 0) && (sum(irlba(f1n,nu=2,nv=2)$d>1e-6)==2) && (nrow(f1n) > 4)) {
        #print("keep going!")
        res = break.cl.sp(f1, a1, xi.loc.labels, ncl, a1.labels,BH,lap=lap,reg=reg,n.min,D=D-1)
        xi.loc.labels = res$xi.loc.labels
        ncl = res$ncl
        L.tree.path <- res$tree.path
        if(length(isol)>0){
          xi.loc.labels[[ncl]] = c(xi.loc.labels[[ncl]],cl.labels.full[isol]) ### attached the isolated nodes in this level with the clusters under the smaller split
          
          path.head <- L.tree.path[length(L.tree.path)]
          path.head <- gsub('[[:digit:]]+', '', path.head)
          iso.path <- paste(path.head,cl.labels.full[isol],sep="")
          L.tree.path <- c(L.tree.path,iso.path)
        }
        
        L.mod.path <- res$mod.path
      } else {
        ncl = ncl + 1
        xi.loc.labels[[ncl]] = a1.labels
        if(length(isol)>0) xi.loc.labels[[ncl]] = c(cl.labels.full[isol],xi.loc.labels[[ncl]])
        L.tree.path <- as.character(xi.loc.labels[[ncl]])
        L.mod.path <- rep(0,length(xi.loc.labels[[ncl]]))
        print('Homogeneous, Branch End')
      }
    } else {
      ncl = ncl + 1
      xi.loc.labels[[ncl]] = a1.labels
      if(length(isol)>0) xi.loc.labels[[ncl]] = c(cl.labels.full[isol],xi.loc.labels[[ncl]])
      L.tree.path <- as.character(xi.loc.labels[[ncl]])
      L.mod.path <- rep(0,length(xi.loc.labels[[ncl]]))
      print('Too small cluster, Branch End')
    }
    f2 = f[-smaller.cluster, -smaller.cluster]
    sum.f2 <- sum(f2-pnull)
    num.f2 <- nrow(f2)
    #pf2 <- sum.f2/(num.f2^2-num.f2)
    #print(paste("pf2 =",pf2))
    
    a2 = per.a[-smaller.cluster]
    a2.labels = cl.labels[-smaller.cluster]
    if(length(dim(f2)) > 0) {
      #nisol = which(rowSums(f2) > 0)
      f2n = f2#[nisol, nisol]
      #print(dim(f2n))
      if(length(dim(f2n)) > 0) {
        n1 <- nrow(f2n)
        phat = sum(f2n)/(n1*(n1-1))
      } else {
        n1 <- 0
      }
    } else if(length(f2) > 0) {
      f2n = diag(rep(f2, 2))
      n1 <- nrow(f2n)
      phat = sum(f2n)/(n1*(n1-1))
    } else {
      n1 <- 0
    }
    if(n1 > 2*n.min) {  ## change to 2*n.min
      #print(dim(f2n))
      if((phat != 1) && (phat != 0) && (sum(irlba(f2n,nu=2,nv=2)$d>1e-6)==2) && (nrow(f2n) > 4)) {
        res = break.cl.sp(f2, a2, xi.loc.labels, ncl, a2.labels,BH,lap=lap,reg=reg,n.min,D=D-1)
        xi.loc.labels = res$xi.loc.labels
        R.tree.path <- res$tree.path
        R.mod.path <- res$mod.path
        ncl = res$ncl
      } else {
        ncl = ncl + 1
        xi.loc.labels[[ncl]] = a2.labels
        R.tree.path <- as.character(a2.labels)
        R.mod.path <- rep(0,length(a2.labels))
        print('Homogenous, Branch End')
      }
    } else {
      ncl = ncl + 1
      xi.loc.labels[[ncl]] = a2.labels
      R.tree.path <- as.character(a2.labels)
      R.mod.path <- rep(0,length(a2.labels))
      print('Too small cluster, Branch End')
    }
    L.tree.path <- paste("L",L.tree.path,sep="/")
    #if(length(isol)>0) R.tree.path <- c(R.tree.path,as.character(cl.labels.full[isol]))
    R.tree.path <- paste("R",R.tree.path,sep="/")
    tree.path <- c("",L.tree.path,R.tree.path)
    
    modularity <- (sum.f1+sum.f2)/sum(f)  ## check modularity, not useful, to be removed
    mod.path <- c(modularity,L.mod.path,R.mod.path)
  } else {
    ncl = ncl + 1
    xi.loc.labels[[ncl]] = cl.labels
    if(length(isol)>0) xi.loc.labels[[ncl]] = c(cl.labels.full[isol],cl.labels)
    tree.path <- c("",as.character(cl.labels))
    mod.path <- c(0,rep(0,length(cl.labels)))
    if(length(isol)>0) tree.path <- c(tree.path,as.character(cl.labels.full[isol]))
    print('One cluster, Branch End, not even started!')
  }
  #print(tree.path)
  return(list(xi.loc.labels = xi.loc.labels, ncl = ncl,tree.path=tree.path,mod.path=mod.path))
}



## HCD-SS algorithm. The population probability matrix is still based on standard HCD
HCD.SS <- function(A,a,K,BH=2,adj=TRUE,plot.it=FALSE,lap=FALSE,reg=FALSE,n.min=25,D=NULL){
  f <- A
  #graph.f <- graph.adjacency(f,mode="undirected")
  #f <- get.adjacency(graph.f, sparse = T)
  n <- nrow(f)
  ncl <- 0
  xi.loc.labels <- list()
  clusters <- break.cl.sp(f, a, xi.loc.labels, ncl, 1:n,BH,lap=lap,reg=reg,n.min,D=D)
  ncl <- clusters$ncl
  xi.loc.labels <- clusters$xi.loc.labels
  labels = rep(0, n)
  for(i in 1:ncl) {
    labels[xi.loc.labels[[i]]] <- i
  }
  cross.table <- table(labels, a)
  tree.path <- clusters$tree.path
  node.tree.path <- strsplit(tree.path,"/")
  node.index <- which(unlist(lapply(node.tree.path,function(x) sum(!is.na(as.numeric(x)))>0))) ### find which path string is for individual nodes
  node.number <- unlist(lapply(node.tree.path,function(x) as.numeric(x)[which(!is.na(as.numeric(x)))])) ### find the specific node name. The non-node string will be removed
  node.dt <- data.table(node.number=node.number,node.index=node.index)
  node.dt2 <- unique(node.dt, by = "node.number")
  node.index <- node.dt2$node.index[sort(node.dt2$node.number,index.return=TRUE)$ix]
  #Binary.Similarity(node.tree.path[[node.index[1]]],node.tree.path[[node.index[1000]]])
  representers <- rep(0,ncl)
  for(i in 1:ncl){
    representers[i] <- which(labels==i)[1]
  }
  community.bin.sim.mat <- matrix(0,ncl,ncl)
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
  Z.mat[cbind(1:n,labels)] <- 1
  node.bin.sim.mat <- Z.mat%*%community.bin.sim.mat%*%t(Z.mat)
  diag(node.bin.sim.mat) <- 0
  
  tree.path <- paste("All",tree.path,sep="/")
  tree.df <- data.frame(V1=1:length(tree.path),pathString=tree.path)
  #print(tree.path)
  tree.df$pathString <- as.character(tree.df$pathString)
  cluster.tree <- as.Node(tree.df)
  block.distance.mat <- matrix(0,ncl,ncl)
  node.distance.mat <- matrix(0,nrow=n,ncol=n)
  for(k1 in 1:(ncl-1)){
    c1.index <- which(labels==k1)
    if(ncl>1){
      for(k2 in (k1+1):ncl){
        c2.index <- which(labels==k2)
        #print(c(k1,k2))
        tmp.dist <- Distance(node1=FindNode(cluster.tree,as.character(c1.index[1])),node2=FindNode(cluster.tree,as.character(c2.index[1])))-2
        block.distance.mat[k1,k2] <- block.distance.mat[k2,k1] <- tmp.dist
        node.distance.mat[c1.index,c2.index] <- tmp.dist
      }
    }
  }
  node.distance.mat <- node.distance.mat+t(node.distance.mat)
  block.distance.mat <- as.dist(block.distance.mat)
  node.distance.mat <- as.dist(node.distance.mat)#as(node.distance.mat,"dgCMatrix")
  NMI <- mi.empirical(as.matrix(cross.table))/entropy.empirical(as.vector(cross.table))
  P.est.node.dmat <- P.est <- NULL
  if(adj){
    P.est <- SBM.estimate(A,labels)
    HCD.P.est <- HCD.pop(A=P.est,a=labels,K=ncl)
    P.est.node.dmat <- HCD.P.est$node.dmat
    P.node.bin.sim.mat <- HCD.P.est$node.bin.sim.mat
  }
  result <- list(labels=labels,true.label=a,cross.table=cross.table,NMI=NMI,ncl=ncl,cluster.tree=cluster.tree,mod.path=clusters$mod.path,block.dmat=block.distance.mat,node.dmat=node.distance.mat,P=P.est,P.node.dmat=P.est.node.dmat,ncl.P=HCD.P.est$ncl,node.bin.sim.mat=node.bin.sim.mat,P.node.bin.sim.mat=P.node.bin.sim.mat,comm.bin.sim.mat=community.bin.sim.mat,tree.path=tree.path)
  
  if(plot.it) plot(as.dendrogram(cluster.tree,heightAttribute="V1"),center=T,leaflab="none")
  return(result)
}




NMI <- function(g1,g2){
  cross.table <- table(g1, g2)
  nmi <- mi.empirical(as.matrix(cross.table))/entropy.empirical(as.vector(cross.table))
  return(nmi)
}