rm(list=ls())
#.libPaths(c(.libPaths(),'/afs/umich.edu/user/t/i/tianxili/R/x86_64-redhat-linux-gnu-library/3.2','/usr/lib64/R/library','/usr/share/R/library'))
#.libPaths(new='/usr/lib64/R/library')
#.libPaths(new='/usr/share/R/library')
if (!require("pacman")) install.packages("pacman")
pacman::p_load(optparse, data.table, vsp)

option_list = list(
  make_option(c("--nsamp"), type="integer", default="256", 
              help="Sample size", metavar="number"),
  make_option(c("--lambda"), type="numeric", default="50", 
              help="Gamma", metavar="number"),
  make_option(c("--k"), type="integer", default="2", 
              help="Cluster", metavar="number"),
  make_option(c("--jobid"), type="character", default="1", 
              help="jobid", metavar="number"),
  make_option(c("--ibatch"), type="numeric", default="1", 
              help="batch_size", metavar="number"),
  make_option(c("--batch_size"), type="numeric", default="2", 
              help="batch_size", metavar="number")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
opt

new_dir <- paste0("../hpcc_output/vsp_", opt$jobid, "/")
dir.create(new_dir)

job_id <- opt$jobid

#library(fpc)
#library(mvtnorm)
#library(ergm)
library(igraph)
library(entropy)
library(e1071)
#library(mclust)
library(Matrix)
library(cluster)
library(irlba)
source("evaluate_estimation.R")
source("FindBeta.R")
source("ArashGen.R")
source("HCD_CodeCheck_1811.R")
source("BSHSBM.R")
source("DistanceMetric.R")
source("BinaryStringSBM.R")
#rm(list=ls())
# require(doParallel)

# registerDoParallel(cores=4)

# DD.seq <- c(2,3,4,5,6)
system.time({
  
  # DD <- DD.seq[KK]
  DD <- opt$k
  II <- opt$batch_size
  K.table <- rep(0,100)
  # layer.n.seq <- rep(c(1,2,3,4),each=50)
  # set.seed(500+KK)
  sp2.acc.seq <- sp4.acc.seq <- mega4.acc.sp.seq <- mega4.acc.HCD.seq <- mega4.acc.HCD2.seq <- mega4.coclust.sp.seq <- mega4.coclust.HCD.seq <- mega4.coclust.HCD2.seq <-
    mega2.acc.sp.seq <- mega2.acc.HCD.seq <- mega2.acc.HCD2.seq <- mega2.coclust.sp.seq <- mega2.coclust.HCD.seq <- mega2.coclust.HCD2.seq <- mega2.nmi.sp.seq <- mega4.nmi.sp.seq <- mega2.nmi.HCD.seq <- mega2.nmi.HCD2.seq <- mega4.nmi.HCD.seq <- mega4.nmi.HCD2.seq <- bin.err.HCD.seq <-bin.err.HCD.P.seq <-bin.err.HCD2.seq <-bin.err.HCD2.P.seq <- bin.err.sp.P.seq <- Standard.FP.HCD.seq <- Standard.FP.HCD2.seq <- Standard.FP.sp.seq <- Scaled.FN.HCD.seq <- Scaled.FP.HCD.seq <- Scaled.FN.HCD2.seq <- Scaled.FP.HCD2.seq <- Scaled.FN.HCD.P.seq <- Scaled.FP.HCD.P.seq <- Scaled.FN.HCD2.P.seq <- Scaled.FP.HCD2.P.seq <- Scaled.FN.sp.seq <- Scaled.FP.sp.seq <-    P.err.HCD.seq <- P.err.HCD2.seq <- P.err.sp.seq <- ncl.P.HCD.seq <- ncl.P.HCD2.seq <- ncl.P.sp.seq <- size.var.seq <- size.ratio.seq <- coclust.HCD.seq <- coclust.HCD2.seq <- coclust.sp.seq <- coclust.true.HCD.seq <- coclust.true.HCD2.seq <- coclust.true.sp.seq <- coclust.false.HCD.seq <- coclust.false.HCD2.seq <- coclust.false.sp.seq <-  nmi.HCD.seq <- nmi.HCD2.seq <- nmi.sp.seq <- ncl.NB.seq <- ncl.HCD.seq <- ncl.HCD2.seq <- K.seq  <- rep(0,II)
  
  mega2.nmi.HCD.varimax.seq <- mega4.nmi.HCD.varimax.seq <- 
    mega2.nmi.sp.varimax.seq <- mega4.nmi.sp.varimax.seq <-
    ncl.P.HCD.varimax.seq <- ncl.P.sp.varimax.seq <-
    P.err.HCD.varimax.seq <- P.err.sp.varimax.seq <-
    bin.err.HCD.varimax.seq <- bin.err.HCD.varimax.P.seq <- bin.err.sp.varimax.P.seq <-
    mega2.coclust.HCD.varimax.seq <- mega2.coclust.sp.varimax.seq <- mega2.acc.HCD.varimax.seq <-
    mega2.acc.sp.varimax.seq <- sp2.varimax.acc.seq   <- 
    mega4.coclust.HCD.varimax.seq <- mega4.coclust.sp.varimax.seq <- mega4.acc.HCD.varimax.seq <-
    mega4.acc.sp.varimax.seq <- sp4.varimax.acc.seq <- 
    nmi.HCD.varimax.seq <- nmi.sp.varimax.seq <- ncl.HCD.varimax.seq <- rep(0,II)
  
  for(I in 1:II){
    
    set.seed((opt$ibatch-1)*opt$batch_size + I)
    
    n <- opt$nsamp
    
    beta <- FindBeta(DD,0.15)
    # beta <- FindBeta(DD,0.15)
    
    tt <- BinaryStringSBM(n,DD,beta=beta,lambda=opt$lambda)
    B <- tt$B
    
    
    mean(tt$B[1,-1])/tt$B[1,1]
    
    
    
    K <- length(unique(tt$label))
    
    K
    
    comm.size <- as.numeric(table(tt$label))
    size.var <- var(comm.size)
    size.ratio <- max(comm.size)/min(comm.size)
    size.var.seq[I] <- size.var
    size.ratio.seq[I] <- size.ratio
    
    K.table[K] <- K.table[K]+1
    
    P <- tt$P
    upper.tri.index <- which(upper.tri(P))
    tmp.rand <- runif(n=length(upper.tri.index))
    A <- matrix(0,n,n)
    A[upper.tri.index[tmp.rand<P[upper.tri.index]]] <- 1
    A <- A+t(A)
    diag(A) <- 0
    mean(colSums(A))
    
    true.hc <- hclust(as.dist(50-tt$comm.sim.mat))
    #png("tmp.png")
    #plot(HCD.hc)
    #dev.off()
    
    true.Mega2 <- cutree(true.hc,k=2)
    true.Mega2.label <- rep(0,n)
    true.Mega2.label <- true.Mega2[tt$label]
    true.Mega4 <- cutree(true.hc,k=4)
    true.Mega4.label <- rep(0,n)
    true.Mega4.label <- true.Mega4[tt$label]
    
    
    # HCD.varimax.result <- HCD.SS(A,tt$label,K=K,BH=2,lap=FALSE,reg=FALSE,varimax = T)
    HCD.varimax.result <- HCD(A,tt$label,K=K,BH=2,lap=FALSE,reg=FALSE,varimax = T)
    sp.varimax.result <- SpectralClust.varimax(A,tt$label,K=HCD.varimax.result$ncl)
    sp.varimax.result2 <- SpectralClust.varimax(A,tt$label,K=2)
    sp.varimax.result4 <- SpectralClust.varimax(A,tt$label,K=4)
    
    HCD.result <- HCD.SS(A,tt$label,K=K,BH=2,lap=FALSE,reg=FALSE)
    HCD.result2 <- HCD(A,tt$label,K=K,2,lap=TRUE,reg=TRUE)
    sp.result <- SpectralClust(A,tt$label,K=HCD.result2$ncl,lap=TRUE,reg=TRUE)
    sp.result2 <- SpectralClust(A,tt$label,K=2,lap=TRUE,reg=TRUE)
    sp.result4 <- SpectralClust(A,tt$label,K=4,lap=TRUE,reg=TRUE)
    #cluster.tree <- HCD.result$cluster.tree
    #labels <- HCD.result$labels
    #ncl <- HCD.result$ncl
    #names(HCD.result)
    
    HCD.hc <- hclust(as.dist(50-HCD.result$comm.bin.sim.mat),method="single")
    HCD.Mega2 <- cutree(HCD.hc,k=2)
    HCD.Mega2.label <- rep(0,n)
    HCD.Mega2.label <- HCD.Mega2[HCD.result$label]
    if(HCD.result$ncl>=4){
      HCD.Mega4 <- cutree(HCD.hc,k=4)
      HCD.Mega4.label <- rep(0,n)
      HCD.Mega4.label <- HCD.Mega4[HCD.result$label]
    }else{
      HCD.Mega4.label <- HCD.Mega2.label
      sub.index <- sample(1:n,n/2)
      HCD.Mega4.label[sub.index] <- HCD.Mega4.label[sub.index]+2
    }
    HCD2.hc <- hclust(as.dist(50-HCD.result2$comm.bin.sim.mat),method="single")
    HCD2.Mega2 <- cutree(HCD2.hc,k=2)
    HCD2.Mega2.label <- rep(0,n)
    HCD2.Mega2.label <- HCD2.Mega2[HCD.result2$label]
    if(HCD.result2$ncl>=4){
      HCD2.Mega4 <- cutree(HCD2.hc,k=4)
      HCD2.Mega4.label <- rep(0,n)
      HCD2.Mega4.label <- HCD2.Mega4[HCD.result2$label]
    }else{
      HCD2.Mega4.label <- HCD2.Mega2.label
      sub.index <- sample(1:n,n/2)
      HCD2.Mega4.label[sub.index] <- HCD2.Mega4.label[sub.index]+2
    }
    
    ### VSP
    HCD.varimax.hc <- hclust(as.dist(50-HCD.varimax.result$comm.bin.sim.mat),method="single")
    HCD.varimax.Mega2 <- cutree(HCD.varimax.hc,k=2)
    HCD.varimax.Mega2.label <- rep(0,n)
    HCD.varimax.Mega2.label <- HCD.varimax.Mega2[HCD.varimax.result$label]
    if(HCD.varimax.result$ncl>=4){
      HCD.varimax.Mega4 <- cutree(HCD.varimax.hc,k=4)
      HCD.varimax.Mega4.label <- rep(0,n)
      HCD.varimax.Mega4.label <- HCD.varimax.Mega4[HCD.varimax.result$label]
    }else{
      HCD.varimax.Mega4.label <- HCD.varimax.Mega2.label
      sub.index <- sample(1:n,n/2)
      HCD.varimax.Mega4.label[sub.index] <- HCD.varimax.Mega4.label[sub.index]+2
    }
    
    sp.hc <- hclust(sp.result$P.node.dmat,method="complete")
    sp.Mega2.label <- cutree(sp.hc,k=2)
    sp2.label <- sp.result2$labels
    if(HCD.result$ncl>=4){
      sp.Mega4.label <- cutree(sp.hc,k=4)
      sp4.label <- sp.result4$labels
    }else{
      sp.Mega4.label <- sp.Mega2.label
      sub.index <- sample(1:n,n/2)
      sp.Mega4.label[sub.index] <- sp.Mega4.label[sub.index]+2
    }
    sp4.label <- sp.result4$labels
    
    ### VSP
    sp.varimax.hc <- hclust(sp.varimax.result$P.node.dmat,method="complete")
    sp.varimax.Mega2.label <- cutree(sp.varimax.hc,k=2)
    sp.varimax2.label <- sp.varimax.result2$labels
    if(HCD.varimax.result$ncl>=4){
      sp.varimax.Mega4.label <- cutree(sp.varimax.hc,k=4)
      sp.varimax4.label <- sp.varimax.result4$labels
    }else{
      sp.varimax.Mega4.label <- sp.varimax.Mega2.label
      sub.index <- sample(1:n,n/2)
      sp.varimax.Mega4.label[sub.index] <- sp.varimax.Mega4.label[sub.index]+2
    }
    sp.varimax4.label <- sp.varimax.result4$labels
    
    
    mega2.nmi.HCD.seq[I] <- NMI(true.Mega2.label,HCD.Mega2.label)
    mega2.nmi.HCD2.seq[I] <- NMI(true.Mega2.label,HCD2.Mega2.label)
    mega2.nmi.HCD.varimax.seq[I] <- NMI(true.Mega2.label,HCD.varimax.Mega2.label)
    mega4.nmi.HCD.seq[I] <- NMI(true.Mega4.label,HCD.Mega4.label)
    mega4.nmi.HCD2.seq[I] <- NMI(true.Mega4.label,HCD2.Mega4.label)
    mega4.nmi.HCD.varimax.seq[I] <- NMI(true.Mega4.label,HCD.varimax.Mega4.label)
    mega2.nmi.sp.seq[I] <- NMI(true.Mega2.label,sp.Mega2.label)
    mega4.nmi.sp.seq[I] <- NMI(true.Mega4.label,sp.Mega4.label)
    mega2.nmi.sp.varimax.seq[I] <- NMI(true.Mega2.label,sp.varimax.Mega2.label)
    mega4.nmi.sp.varimax.seq[I] <- NMI(true.Mega4.label,sp.varimax.Mega4.label)
    
    ncl.P.HCD.seq[I] <- HCD.result$ncl.P
    ncl.P.HCD2.seq[I] <- HCD.result2$ncl.P
    ncl.P.HCD.varimax.seq[I] <- HCD.varimax.result$ncl.P
    ncl.P.sp.seq[I] <- sp.result$ncl.P
    ncl.P.sp.varimax.seq[I] <- sp.varimax.result$ncl.P
    
    P.err.HCD.seq[I] <- norm(P-HCD.result$P,"F")^2/norm(P,"F")^2
    P.err.HCD.varimax.seq[I] <- norm(P-HCD.varimax.result$P,"F")^2/norm(P,"F")^2
    P.err.HCD2.seq[I] <- norm(P-HCD.result2$P,"F")^2/norm(P,"F")^2
    P.err.sp.seq[I] <- norm(P-sp.result$P,"F")^2/norm(P,"F")^2
    P.err.sp.varimax.seq[I] <- norm(P-sp.varimax.result$P,"F")^2/norm(P,"F")^2
    
    
    bin.err.HCD.seq[I] <- norm(tt$node.sim.mat-HCD.result$node.bin.sim.mat,"F")^2/norm(tt$node.sim.mat,"F")^2
    bin.err.HCD.P.seq[I] <- norm(tt$node.sim.mat-HCD.result$P.node.bin.sim.mat,"F")^2/norm(tt$node.sim.mat,"F")^2
    bin.err.HCD.varimax.seq[I] <- norm(tt$node.sim.mat-HCD.varimax.result$node.bin.sim.mat,"F")^2/norm(tt$node.sim.mat,"F")^2
    bin.err.HCD.varimax.P.seq[I] <- norm(tt$node.sim.mat-HCD.varimax.result$P.node.bin.sim.mat,"F")^2/norm(tt$node.sim.mat,"F")^2
    bin.err.HCD2.seq[I] <- norm(tt$node.sim.mat-HCD.result2$node.bin.sim.mat,"F")^2/norm(tt$node.sim.mat,"F")^2
    bin.err.HCD2.P.seq[I] <- norm(tt$node.sim.mat-HCD.result2$P.node.bin.sim.mat,"F")^2/norm(tt$node.sim.mat,"F")^2
    bin.err.sp.P.seq[I] <- norm(tt$node.sim.mat-sp.result$P.node.bin.sim.mat,"F")^2/norm(tt$node.sim.mat,"F")^2
    bin.err.sp.varimax.P.seq[I] <- norm(tt$node.sim.mat-sp.varimax.result$P.node.bin.sim.mat,"F")^2/norm(tt$node.sim.mat,"F")^2
    
    
    
    
    Z0 <- matrix(0,nrow=n,ncol=K)
    Z0[cbind(1:n,tt$label)] <- 1
    G0 <- Z0%*%t(Z0)
    Z11 <- matrix(0,nrow=n,ncol=HCD.result$ncl)
    Z2 <- Z12 <- matrix(0,nrow=n,ncol=HCD.result2$ncl)
    Z2 <- matrix(0,nrow=n,ncol=length(unique(sp.result$labels)))
    
    Z11[cbind(1:n,HCD.result$labels)] <- 1
    Z12[cbind(1:n,HCD.result2$labels)] <- 1
    Z2[cbind(1:n,sp.result$labels)] <- 1
    G12 <- Z12%*%t(Z12)
    G2 <- Z2%*%t(Z2)
    G11 <- Z11%*%t(Z11)
    
    
    coclust.HCD.seq[I] <- mean(G0==G11)
    coclust.HCD2.seq[I] <- mean(G0==G12)
    coclust.sp.seq[I] <- mean(G0==G2)
    which.true.coclust <- which(G0==1)
    coclust.true.HCD.seq[I] <- mean(G11[which.true.coclust])
    coclust.true.HCD2.seq[I] <- mean(G12[which.true.coclust])
    coclust.true.sp.seq[I] <- mean(G2[which.true.coclust])
    which.false.coclust <- which(G0==0)
    coclust.false.HCD.seq[I] <- mean(G11[which.false.coclust]==1)
    coclust.false.HCD2.seq[I] <- mean(G12[which.false.coclust]==1)
    coclust.false.sp.seq[I] <- mean(G2[which.false.coclust]==1)
    
    
    mega2.Z0 <- matrix(0,nrow=n,ncol=2)
    mega2.Z0[cbind(1:n,true.Mega2.label)] <- 1
    mega2.G0 <- mega2.Z0%*%t(mega2.Z0)
    mega2.HCD.Z <- matrix(0,nrow=n,ncol=2)
    sp2.Z <- mega2.sp.Z <- mega2.HCD2.Z <- matrix(0,nrow=n,ncol=2)
    sp2.varimax.Z <- mega2.sp.varimax.Z <- mega2.HCD.varimax.Z <- matrix(0,nrow=n,ncol=2)
    
    mega2.HCD.Z[cbind(1:n,HCD.Mega2.label)] <- 1
    mega2.HCD.varimax.Z[cbind(1:n,HCD.varimax.Mega2.label)] <- 1
    mega2.HCD2.Z[cbind(1:n,HCD2.Mega2.label)] <- 1
    mega2.sp.Z[cbind(1:n,sp.Mega2.label)] <- 1
    sp2.Z[cbind(1:n,sp2.label)] <- 1
    mega2.sp.varimax.Z[cbind(1:n,sp.varimax.Mega2.label)] <- 1
    sp2.varimax.Z[cbind(1:n,sp.varimax2.label)] <- 1
    
    mega2.HCD.G <- mega2.HCD.Z%*%t(mega2.HCD.Z)
    mega2.HCD.varimax.G <- mega2.HCD.varimax.Z%*%t(mega2.HCD.varimax.Z)
    mega2.HCD2.G <- mega2.HCD2.Z%*%t(mega2.HCD2.Z)
    mega2.sp.G <- mega2.sp.Z%*%t(mega2.sp.Z)
    mega2.sp.varimax.G <- mega2.sp.varimax.Z%*%t(mega2.sp.varimax.Z)
    
    
    
    mega2.coclust.HCD.seq[I] <- mean(mega2.G0==mega2.HCD.G)
    mega2.coclust.HCD.varimax.seq[I] <- mean(mega2.G0==mega2.HCD.varimax.G)
    mega2.coclust.HCD2.seq[I] <- mean(mega2.G0==mega2.HCD2.G)
    mega2.coclust.sp.seq[I] <- mean(mega2.G0==mega2.sp.G)
    mega2.coclust.sp.varimax.seq[I] <- mean(mega2.G0==mega2.sp.varimax.G)
    mega2.acc.HCD.seq[I] <- accuracy(mega2.HCD.Z,mega2.Z0)
    mega2.acc.HCD.varimax.seq[I] <- accuracy(mega2.HCD.varimax.Z,mega2.Z0)
    mega2.acc.HCD2.seq[I] <- accuracy(mega2.HCD2.Z,mega2.Z0)
    mega2.acc.sp.varimax.seq[I] <- accuracy(mega2.sp.varimax.Z,mega2.Z0)
    sp2.acc.seq[I] <- accuracy(sp2.Z,mega2.Z0)
    sp2.varimax.acc.seq[I] <- accuracy(sp2.varimax.Z,mega2.Z0)
    
    
    
    mega4.Z0 <- matrix(0,nrow=n,ncol=4)
    mega4.Z0[cbind(1:n,true.Mega4.label)] <- 1
    mega4.G0 <- mega4.Z0%*%t(mega4.Z0)
    mega4.HCD.Z <- matrix(0,nrow=n,ncol=4)
    sp4.Z <- mega4.sp.Z <- mega4.HCD2.Z <- matrix(0,nrow=n,ncol=4)
    sp4.varimax.Z <- mega4.sp.varimax.Z <- mega4.HCD.varimax.Z <- matrix(0,nrow=n,ncol=4)
    
    mega4.HCD.Z[cbind(1:n,HCD.Mega4.label)] <- 1
    mega4.HCD.varimax.Z[cbind(1:n,HCD.varimax.Mega4.label)] <- 1
    mega4.HCD2.Z[cbind(1:n,HCD2.Mega4.label)] <- 1
    mega4.sp.Z[cbind(1:n,sp.Mega4.label)] <- 1
    sp4.Z[cbind(1:n,sp4.label)] <- 1
    mega4.sp.varimax.Z[cbind(1:n,sp.varimax.Mega4.label)] <- 1
    sp4.varimax.Z[cbind(1:n,sp.varimax4.label)] <- 1
    
    mega4.HCD.G <- mega4.HCD.Z%*%t(mega4.HCD.Z)
    mega4.HCD.varimax.G <- mega4.HCD.varimax.Z%*%t(mega4.HCD.varimax.Z)
    mega4.HCD2.G <- mega4.HCD2.Z%*%t(mega4.HCD2.Z)
    mega4.sp.G <- mega4.sp.Z%*%t(mega4.sp.Z)
    mega4.sp.varimax.G <- mega4.sp.varimax.Z%*%t(mega4.sp.varimax.Z)
    
    
    mega4.coclust.HCD.seq[I] <- mean(mega4.G0==mega4.HCD.G)
    mega4.coclust.HCD.varimax.seq[I] <- mean(mega4.G0==mega4.HCD.varimax.G)
    mega4.coclust.HCD2.seq[I] <- mean(mega4.G0==mega4.HCD2.G)
    mega4.coclust.sp.seq[I] <- mean(mega4.G0==mega4.sp.G)
    mega4.coclust.sp.varimax.seq[I] <- mean(mega4.G0==mega4.sp.varimax.G)
    mega4.acc.HCD.seq[I] <- accuracy(mega4.HCD.Z,mega4.Z0)
    mega4.acc.HCD.varimax.seq[I] <- accuracy(mega4.HCD.varimax.Z,mega4.Z0)
    mega4.acc.HCD2.seq[I] <- accuracy(mega4.HCD2.Z,mega4.Z0)
    mega4.acc.sp.varimax.seq[I] <- accuracy(mega4.sp.varimax.Z,mega4.Z0)
    sp4.acc.seq[I] <- accuracy(sp4.Z,mega4.Z0)
    sp4.varimax.acc.seq[I] <- accuracy(sp4.varimax.Z,mega4.Z0)
    
    
    
    
    
    
    #Bk(tree1=t1, tree2=t2, k = 2)
    nmi.HCD.seq[I] <- HCD.result$NMI
    nmi.HCD.varimax.seq[I] <- HCD.varimax.result$NMI
    nmi.HCD2.seq[I] <- HCD.result2$NMI
    nmi.sp.seq[I] <- sp.result$NMI
    nmi.sp.varimax.seq[I] <- sp.varimax.result$NMI
    
    K.seq[I] <- K
    ncl.HCD.seq[I] <- HCD.result$ncl
    ncl.HCD.varimax.seq[I] <- HCD.varimax.result$ncl
    ncl.HCD2.seq[I] <- HCD.result2$ncl
    ncl.NB.seq[I] <- NB.estimate(A)
    #num.edge.HCD.seq[I] <- count.edge.dend(t2)
    
    tmp <- list(nmi.HCD.seq=nmi.HCD.seq,nmi.HCD2.seq=nmi.HCD2.seq, 
                nmi.sp.seq=nmi.sp.seq,  ncl.NB.seq=ncl.NB.seq,  
                ncl.HCD.seq = ncl.HCD.seq, ncl.HCD2.seq = ncl.HCD2.seq,  
                K.seq=K.seq,coclust.HCD.seq=coclust.HCD.seq, 
                coclust.HCD2.seq=coclust.HCD2.seq,
                coclust.sp.seq=coclust.sp.seq,coclust.true.HCD.seq=coclust.true.HCD.seq,
                coclust.true.HCD2.seq=coclust.true.HCD2.seq,
                coclust.true.sp.seq=coclust.true.sp.seq,
                coclust.false.HCD.seq=coclust.false.HCD.seq,
                coclust.false.HCD2.seq=coclust.false.HCD2.seq,coclust.false.sp.seq=coclust.false.sp.seq,size.var=size.var.seq,size.ratio=size.ratio.seq,ncl.P.HCD=ncl.P.HCD.seq,ncl.P.HCD2=ncl.P.HCD2.seq,ncl.P.sp=ncl.P.sp.seq,P.err.HCD=P.err.HCD.seq,P.err.HCD2=P.err.HCD2.seq,P.err.sp=P.err.sp.seq,Scaled.FP.HCD=Scaled.FP.HCD.seq,Scaled.FP.HCD2=Scaled.FP.HCD2.seq,Scaled.FP.sp = Scaled.FP.sp.seq,Standard.FP.HCD=Standard.FP.HCD.seq,Standard.FP.HCD2=Standard.FP.HCD2.seq,Standard.FP.sp = Standard.FP.sp.seq,bin.err.HCD.seq=bin.err.HCD.seq,bin.err.HCD2.seq=bin.err.HCD2.seq,bin.err.HCD.P.seq=bin.err.HCD.P.seq,bin.err.HCD2.P.seq=bin.err.HCD2.P.seq,bin.err.sp.P.seq=bin.err.sp.P.seq,mega2.nmi.sp.seq=mega2.nmi.sp.seq, mega4.nmi.sp.seq=mega4.nmi.sp.seq , mega2.nmi.HCD.seq=mega2.nmi.HCD.seq , mega2.nmi.HCD2.seq=mega2.nmi.HCD2.seq , mega4.nmi.HCD.seq=mega4.nmi.HCD.seq , mega4.nmi.HCD2.seq=mega4.nmi.HCD2.seq, mega2.coclust.HCD.seq=mega2.coclust.HCD.seq, mega2.coclust.HCD2.seq=mega2.coclust.HCD2.seq,mega2.coclust.sp.seq=mega2.coclust.sp.seq,mega2.acc.HCD.seq=mega2.acc.HCD.seq, mega2.acc.HCD2.seq=mega2.acc.HCD2.seq,mega2.acc.sp.seq=mega2.acc.sp.seq, mega4.coclust.HCD.seq=mega4.coclust.HCD.seq, mega4.coclust.HCD2.seq=mega4.coclust.HCD2.seq,mega4.coclust.sp.seq=mega4.coclust.sp.seq,mega4.acc.HCD.seq=mega4.acc.HCD.seq, mega4.acc.HCD2.seq=mega4.acc.HCD2.seq,mega4.acc.sp.seq=mega4.acc.sp.seq,sp2.acc.seq=sp2.acc.seq,sp4.acc.seq=sp4.acc.seq,
                # varimax
                mega2.nmi.HCD.varimax.seq = mega2.nmi.HCD.varimax.seq,
                mega4.nmi.HCD.varimax.seq = mega4.nmi.HCD.varimax.seq,
                mega2.nmi.sp.varimax.seq = mega2.nmi.sp.varimax.seq,
                mega4.nmi.sp.varimax.seq = mega4.nmi.sp.varimax.seq,
                ncl.P.HCD.varimax.seq = ncl.P.HCD.varimax.seq,
                ncl.P.sp.varimax.seq = ncl.P.sp.varimax.seq,
                P.err.HCD.varimax.seq = P.err.HCD.varimax.seq,
                P.err.sp.varimax.seq = P.err.sp.varimax.seq,
                bin.err.HCD.varimax.seq = bin.err.HCD.varimax.seq,
                bin.err.HCD.varimax.P.seq = bin.err.HCD.varimax.P.seq,
                bin.err.sp.varimax.P.seq = bin.err.sp.varimax.P.seq,
                mega2.coclust.HCD.varimax.seq = mega2.coclust.HCD.varimax.seq,
                mega2.coclust.sp.varimax.seq = mega2.coclust.sp.varimax.seq,
                mega2.acc.HCD.varimax.seq = mega2.acc.HCD.varimax.seq,
                mega2.acc.sp.varimax.seq = mega2.acc.sp.varimax.seq,
                sp2.varimax.acc.seq   = sp2.varimax.acc.seq  ,
                mega4.coclust.HCD.varimax.seq = mega4.coclust.HCD.varimax.seq,
                mega4.coclust.sp.varimax.seq = mega4.coclust.sp.varimax.seq,
                mega4.acc.HCD.varimax.seq = mega4.acc.HCD.varimax.seq,
                mega4.acc.sp.varimax.seq = mega4.acc.sp.varimax.seq,
                sp4.varimax.acc.seq = sp4.varimax.acc.seq,
                nmi.HCD.varimax.seq = nmi.HCD.varimax.seq,
                nmi.sp.varimax.seq = nmi.sp.varimax.seq,
                ncl.HCD.varimax.seq = ncl.HCD.varimax.seq)
    
    ret_list <- as.data.table(tmp)
    ret_list[, id := (opt$ibatch-1)*opt$batch_size + (1:opt$batch_size)]
    
    print(ret_list)
    
    fwrite(ret_list, paste0(new_dir,
                            "btsbm_vsp",
                            "_n", opt$nsamp,
                            "_lambda", opt$lambda,
                            "_k", opt$k,
                            "_s", (opt$ibatch-1)*opt$batch_size + 1,
                            ".csv"))
    
  }
  
  tmp <- list(nmi.HCD.seq=nmi.HCD.seq,nmi.HCD2.seq=nmi.HCD2.seq, 
              nmi.sp.seq=nmi.sp.seq,  ncl.NB.seq=ncl.NB.seq,  
              ncl.HCD.seq = ncl.HCD.seq, ncl.HCD2.seq = ncl.HCD2.seq,  
              K.seq=K.seq,coclust.HCD.seq=coclust.HCD.seq, 
              coclust.HCD2.seq=coclust.HCD2.seq,
              coclust.sp.seq=coclust.sp.seq,coclust.true.HCD.seq=coclust.true.HCD.seq,
              coclust.true.HCD2.seq=coclust.true.HCD2.seq,
              coclust.true.sp.seq=coclust.true.sp.seq,
              coclust.false.HCD.seq=coclust.false.HCD.seq,
              coclust.false.HCD2.seq=coclust.false.HCD2.seq,coclust.false.sp.seq=coclust.false.sp.seq,size.var=size.var.seq,size.ratio=size.ratio.seq,ncl.P.HCD=ncl.P.HCD.seq,ncl.P.HCD2=ncl.P.HCD2.seq,ncl.P.sp=ncl.P.sp.seq,P.err.HCD=P.err.HCD.seq,P.err.HCD2=P.err.HCD2.seq,P.err.sp=P.err.sp.seq,Scaled.FP.HCD=Scaled.FP.HCD.seq,Scaled.FP.HCD2=Scaled.FP.HCD2.seq,Scaled.FP.sp = Scaled.FP.sp.seq,Standard.FP.HCD=Standard.FP.HCD.seq,Standard.FP.HCD2=Standard.FP.HCD2.seq,Standard.FP.sp = Standard.FP.sp.seq,bin.err.HCD.seq=bin.err.HCD.seq,bin.err.HCD2.seq=bin.err.HCD2.seq,bin.err.HCD.P.seq=bin.err.HCD.P.seq,bin.err.HCD2.P.seq=bin.err.HCD2.P.seq,bin.err.sp.P.seq=bin.err.sp.P.seq,mega2.nmi.sp.seq=mega2.nmi.sp.seq, mega4.nmi.sp.seq=mega4.nmi.sp.seq , mega2.nmi.HCD.seq=mega2.nmi.HCD.seq , mega2.nmi.HCD2.seq=mega2.nmi.HCD2.seq , mega4.nmi.HCD.seq=mega4.nmi.HCD.seq , mega4.nmi.HCD2.seq=mega4.nmi.HCD2.seq, mega2.coclust.HCD.seq=mega2.coclust.HCD.seq, mega2.coclust.HCD2.seq=mega2.coclust.HCD2.seq,mega2.coclust.sp.seq=mega2.coclust.sp.seq,mega2.acc.HCD.seq=mega2.acc.HCD.seq, mega2.acc.HCD2.seq=mega2.acc.HCD2.seq,mega2.acc.sp.seq=mega2.acc.sp.seq, mega4.coclust.HCD.seq=mega4.coclust.HCD.seq, mega4.coclust.HCD2.seq=mega4.coclust.HCD2.seq,mega4.coclust.sp.seq=mega4.coclust.sp.seq,mega4.acc.HCD.seq=mega4.acc.HCD.seq, mega4.acc.HCD2.seq=mega4.acc.HCD2.seq,mega4.acc.sp.seq=mega4.acc.sp.seq,sp2.acc.seq=sp2.acc.seq,sp4.acc.seq=sp4.acc.seq,
              # varimax
              mega2.nmi.HCD.varimax.seq = mega2.nmi.HCD.varimax.seq,
              mega4.nmi.HCD.varimax.seq = mega4.nmi.HCD.varimax.seq,
              mega2.nmi.sp.varimax.seq = mega2.nmi.sp.varimax.seq,
              mega4.nmi.sp.varimax.seq = mega4.nmi.sp.varimax.seq,
              ncl.P.HCD.varimax.seq = ncl.P.HCD.varimax.seq,
              ncl.P.sp.varimax.seq = ncl.P.sp.varimax.seq,
              P.err.HCD.varimax.seq = P.err.HCD.varimax.seq,
              P.err.sp.varimax.seq = P.err.sp.varimax.seq,
              bin.err.HCD.varimax.seq = bin.err.HCD.varimax.seq,
              bin.err.HCD.varimax.P.seq = bin.err.HCD.varimax.P.seq,
              bin.err.sp.varimax.P.seq = bin.err.sp.varimax.P.seq,
              mega2.coclust.HCD.varimax.seq = mega2.coclust.HCD.varimax.seq,
              mega2.coclust.sp.varimax.seq = mega2.coclust.sp.varimax.seq,
              mega2.acc.HCD.varimax.seq = mega2.acc.HCD.varimax.seq,
              mega2.acc.sp.varimax.seq = mega2.acc.sp.varimax.seq,
              sp2.varimax.acc.seq   = sp2.varimax.acc.seq  ,
              mega4.coclust.HCD.varimax.seq = mega4.coclust.HCD.varimax.seq,
              mega4.coclust.sp.varimax.seq = mega4.coclust.sp.varimax.seq,
              mega4.acc.HCD.varimax.seq = mega4.acc.HCD.varimax.seq,
              mega4.acc.sp.varimax.seq = mega4.acc.sp.varimax.seq,
              sp4.varimax.acc.seq = sp4.varimax.acc.seq,
              nmi.HCD.varimax.seq = nmi.HCD.varimax.seq,
              nmi.sp.varimax.seq = nmi.sp.varimax.seq,
              ncl.HCD.varimax.seq = ncl.HCD.varimax.seq)
  
  ret_list <- as.data.table(tmp)
  ret_list[, id := (opt$ibatch-1)*opt$batch_size + (1:opt$batch_size)]
  
  ret_list
  
  fwrite(ret_list, paste0(new_dir,
                           "btsbm_vsp",
                           "_n", opt$nsamp,
                           "_lambda", opt$lambda,
                           "_k", opt$k,
                           "_s", (opt$ibatch-1)*opt$batch_size + 1,
                           ".csv"))
  
  # save(result,file="BTSBM-Geometric-VaryingK-FixOutIn-withSS-MetaCommunity.Rda")
  
})
