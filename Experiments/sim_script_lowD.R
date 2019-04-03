###############
### Copyright
##############
# EqVarDAG methods and simulations:
# Copyright (c) 2018 - 2019  Wenyu Chen [wenyuc@uw.edu]
# All rights reserved.  See the file COPYING for license terms. 
# GDS methods:
# Copyright (c) 2010 - 2012  Jonas Peters  [peters@stat.math.ethz.ch]
# All rights reserved.  See the file GDS/COPYING for license terms. 

###############
### Simulation
##############

# Low dimensional DAGs with m variables and n samples
# generated as described in the main paper, 
# simulation with N replications, run the code 
#         simulation_lowD(m,n,N)

###############
# load packages
###############
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
cat("Initializing...")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(foreach,doSNOW,doParallel,
               lattice,bnlearn,matrixcalc,Kendall,glmnet)

wd<-getwd()
source("../EqVarDAG_BU.R")
source("../EqVarDAG_TD.R")
setwd("../GDS/startups/")
source("startupGDS.R")
setwd(wd)
source("../GDS/util_DAGs/hammingDistance.R")

###############
### helper functions
###############

# Visualize the adjacency matrix
plot_matrix<-function(M){
  pl<-levelplot(M[ncol(M):1,1:ncol(M)],
                xaxt='n',yaxt='n', 
                col.regions=c(0,1),colorkey=FALSE,xlab=F,ylab=F, 
                scales=list(alternating=0))
  return(pl)
}

# Convert adjacency matrix to rank of ordering
Adj_to_TO_rank<-function(A){
  p=dim(A)[2]
  TR<-rep(0,p)
  CS<-colSums(A)
  r<-0
  while (sum(CS)){
    s<-which(TR+CS==0)
    if (!length(s)){break}
    if (length(s)==1){
      TR[s]<-r+1
      CS<-CS-A[s,]
      r<-r+1
    } else {
      TR[s]<-r+length(s)/3
      CS<-CS-colSums(A[s,])
      r<-r+length(s)
    }
  }
  TR[which(TR==0)]<-p
  return(TR)
}

###############
### Generate data
###############
# Generate random DAG
randomDAG2 <- function(p,probConnect)
  # This function is modified from randomDAG2 function by Jonas Peters
{
  DAG <- diag(rep(0,p))
  causalOrder <- sample(p)
  for(i in 1:(p-2))
  {
    node <- causalOrder[i]
    possibleParents <- causalOrder[(i+1):p]
    numberParents <- rbinom(n=1,size=(p-i),prob=probConnect)
    Parents <- sample(x = possibleParents, size = numberParents, replace = FALSE)
    Parents<-c(causalOrder[i+1],Parents)
    DAG[Parents,node] <- rep(1,numberParents+1)
  }
  # Sample does not work properly when choosing from sets with one element. We thus consider the last case separately.  
  node <- causalOrder[p-1]
  ParentYesNo <- rbinom(n=1,size=1,prob=probConnect)
  DAG[causalOrder[p],node] <- 1
  return(list(DAG=DAG,TO=rev(causalOrder)))
}

# Generate data from SEM
Bmin<-0.3
get_DAGdata<-function(n,p,pc){
  D<-randomDAG2(p,pc)
  truth<-D$DAG
  TO<-D$TO
  errs <- matrix(rnorm(p * n), nrow = p, ncol = n)
  B<-t(truth)
  B[B==1]<-runif(sum(truth),Bmin,1)*(2*rbinom(sum(truth),1,0.5)-1)
  X <- solve(diag(rep(1, p)) - B, errs)
  X <- t(X)
  return(list(truth=truth,B=B,X=X,TO=TO))
}

###############
### Experiment
##############
do_one<-function(p,n,pc){
  result<-matrix(0,6,4)
  # Data generation
  dat<-get_DAGdata(n,p,pc)
  truth<-dat$truth
  X<-dat$X
  TO<-dat$TO
  
  # TD
  ptm<-proc.time()
  sresult<-EqVarDAG_TD(X)
  sx<-sresult$adj
  result[2,1]<- (proc.time()-ptm)[1]
  result[1,1]<- hammingDistance(sx,truth)
  result[3,1]<- Kendall(
    sapply(1:p,function(i){
      which(TO==i)}),
    sapply(1:p,function(i){
      which(sresult$TO==i)}))$tau[1]
  result[4,1]<- ifelse(sum(truth),sum(truth*sx)/sum(truth),0) # power
  result[5,1]<- ifelse(sum(truth),sum(truth*t(sx))/sum(truth),0) # flipped (power+)
  result[6,1]<- ifelse(sum(sx),1-sum(truth*sx)/sum(sx),0) # FDR
  sxtd <- sx
  # BU
  ptm<-proc.time()
  sresult<-EqVarDAG_BU(X)
  sx<-sresult$adj
  result[2,2]<- (proc.time()-ptm)[1]
  result[1,2]<- hammingDistance(sx,truth)
  result[3,2]<- Kendall(
    sapply(1:p,function(i){
      which(TO==i)}),
    sapply(1:p,function(i){
      which(sresult$TO==i)}))$tau[1]
  result[4,2]<- ifelse(sum(truth),sum(truth*sx)/sum(truth),0) # power
  result[5,2]<- ifelse(sum(truth),sum(truth*t(sx))/sum(truth),0) # flipped (power+)
  result[6,2]<- ifelse(sum(sx),1-sum(truth*sx)/sum(sx),0) # FDR
  # GDS
  ptm<-proc.time()
  gx<-GDS(X)$Adj
  result[2,3]<- (proc.time()-ptm)[1]
  result[1,3]<- hammingDistance(gx,truth)
  result[3,3]<- Kendall(
    sapply(1:p,function(i){
      which(TO==i)}),
    Adj_to_TO_rank(gx))$tau[1]
  result[4,3]<- ifelse(sum(truth),sum(truth*gx)/sum(truth),0) # power
  result[5,3]<- ifelse(sum(truth),sum(truth*t(gx))/sum(truth),0) # flipped (power+)
  result[6,3]<- ifelse(sum(gx),1-sum(truth*gx)/sum(gx),0) # FDR
  # GDS-warm
  ptm<-proc.time()
  gx<-GDS(X,startAt = 'Warm',initG = sxtd,smallK=T)$Adj
  result[2,4]<- (proc.time()-ptm)[1]
  result[1,4]<- hammingDistance(gx,truth)
  result[3,4]<- Kendall(
    sapply(1:p,function(i){
      which(TO==i)}),
    Adj_to_TO_rank(gx))$tau[1]
  result[4,4]<- ifelse(sum(truth),sum(truth*gx)/sum(truth),0) # power
  result[5,4]<- ifelse(sum(truth),sum(truth*t(gx))/sum(truth),0) # flipped (power+)
  result[6,4]<- ifelse(sum(gx),1-sum(truth*gx)/sum(gx),0) # FDR
  return(result)
}

sim_run<-function(p,pc,n,N,par){
  # p:  number of nodes
  # pc: probability of connection
  # n:  number of samples
  # N:  number of simulation runs
  if (par==T){
    cl <- makeCluster(detectCores()-1)
    registerDoSNOW(cl)
    print(detectCores())
    results <- 
      foreach(i = 1:N, 
              .export = ls(globalenv()),
              .packages = c("glmnet","matrixcalc","Kendall")) %dopar% {
                do_one(p,n,pc)
              }
    stopCluster(cl) 
  } else {
    pb <- txtProgressBar(min = 0, max = N, style = 3)
    result<-matrix(0,6,4)
    for(i in 1:N) {
      results <- 
        foreach(i = 1:N, 
                .export = ls(globalenv()),
                .packages = c("glmnet","matrixcalc","Kendall")) %do% {
                  do_one(p,n,pc)
                }
      setTxtProgressBar(pb, i)
    }
    close(pb)
  }
  result<-Reduce('+',results)
  result<-rbind(result,
                apply(FUN = sd,MARGIN = 2,
                      X = Reduce('rbind',
                                 lapply(results,
                                        function(x)x[1,]))))
  result<-rbind(result,
                apply(FUN = sd,MARGIN = 2,
                      X = Reduce('rbind',
                                 lapply(results,
                                        function(x)x[3,]))))
  return(result/N)
}

simulation_lowD<-function(m,n,N){
  set.seed(1)
  results<-array(0,c(2,4));
  # return sparse result
  cat('Setting Sparse: m=', m, " n=",n, "N=",N,"Bmin",Bmin, "\n")
  result<-sim_run(m,3/(2*m-2),n,N,T);
  sink(paste("WarmGDS_ChainDAG_Simresults_p=",m,"n=",n,"Bmin",Bmin,"_sparse.txt"))
  cat('Setting', 'sparse',': p=', m, " n=",n,"N=",N, "\n")
  cat('\tTD\tBU\tGDS\tGDS_warm\n')
  cat('Avg. Dist\t', round(result[1,],4),'\n')
  cat('SD . Dist\t', round(result[7,],4),'\n')
  cat('Avg. tau \t', round(result[3,],4),'\n')
  cat('SD . tau \t', round(result[8,],4),'\n')
  cat('Avg. time\t', round(result[2,],4),'\n')
  cat('Avg.Power\t', round(result[4,],4),'\n')
  cat('Avg. Flip\t', round(result[5,],4),'\n')
  cat('Avg. FDR \t', round(result[6,],4),'\n')
  sink()
  # return dense result
  cat('Setting Dense: m=', m, " n=",n, "N=",N,"Bmin",Bmin, "\n")
  result<-sim_run(m,0.3,n,N,T);
  sink(paste("WarmGDS_ChainDAG_Simresults_p=",m,"n=",n,"Bmin",Bmin,"_dense.txt"))
  cat('Setting', 'dense',': p=', m, " n=",n,"N=",N, "\n")
  cat('\tTD\tBU\tGDS\tGDS_warm\n')
  cat('Avg. Dist\t', round(result[1,],4),'\n')
  cat('SD . Dist\t', round(result[7,],4),'\n')
  cat('Avg. tau \t', round(result[3,],4),'\n')
  cat('SD . tau \t', round(result[8,],4),'\n')
  cat('Avg. time\t', round(result[2,],4),'\n')
  cat('Avg.Power\t', round(result[4,],4),'\n')
  cat('Avg. Flip\t', round(result[5,],4),'\n')
  cat('Avg. FDR \t', round(result[6,],4),'\n')
  sink()
}
cat("done.\n")