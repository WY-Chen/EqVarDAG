library(EqVarDAG)
library(Kendall)
hammingDistance <- function(G1,G2)
  # hammingDistance(G1,G2)
  #
  # Computes Hamming Distance between DAGs G1 and G2 with SHD(->,<-) = 1!!!!
  #
  # INPUT:  G1, G2     adjacency graph containing only zeros and ones: (i,j)=1 means edge from X_i to X_j.
  #
  # OUTPUT: hammingDis Hamming Distance between G1 and G2
  #
  # Copyright (c) 2012-2013  Jonas Peters [peters@stat.math.ethz.ch]
  # All rights reserved.  See the file COPYING for license terms.
{
  allMistakesOne <- FALSE
  if(allMistakesOne)
  {
    Gtmp <- (G1+G2)%%2
    Gtmp <- Gtmp + t(Gtmp)
    nrReversals <- sum(Gtmp == 2)/2
    nrInclDel <- sum(Gtmp == 1)/2
    hammingDis <- nrReversals + nrInclDel
  } else
  {
    hammingDis <- sum(abs(G1 - G2))
    # correction: dist(-,.) = 1, not 2
    hammingDis <- hammingDis - 0.5*sum(G1 * t(G1) * (1-G2) * t(1-G2) + G2 * t(G2) * (1-G1) * t(1-G1))
  }
  return(hammingDis)
}


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

randomDAG2_er <- function(p,probConnect)
{
  # This function is modified from randomDAG2 function by Jonas Peters
  DAG <- diag(rep(0,p))
  causalOrder <- sample(p)
  for(i in 3:(p))
  {
    node <- causalOrder[i]
    possibleParents <- causalOrder[1:(i-1)]
    Parents <- possibleParents[rbinom(length(possibleParents),1,probConnect)==1]
    DAG[Parents,node] <- rep(1,length(Parents))
  }
  node <- causalOrder[p-1]
  ParentYesNo <- rbinom(n=1,size=1,prob=probConnect)
  DAG[causalOrder[1],causalOrder[2]] <- 1

  return(list(DAG=DAG,TO=causalOrder))
}
# Chain graph
randomDAG2_chain <- function(p,probConnect)
{
  # This function is modified from randomDAG2 function by Jonas Peters
  DAG <- diag(rep(0,p))
  causalOrder <- sample(p)
  DAG[causalOrder[1],causalOrder[2]] <- 1
  for(i in 3:(p))
  {
    node <- causalOrder[i]
    possibleParents <- causalOrder[1:(i-1)]
    possibleParents <- possibleParents[which(rowSums(DAG[possibleParents,])<4)]
    if (length(possibleParents)>0){
      Parents <- sample(possibleParents,min(length(possibleParents),2))
      DAG[Parents,node] <- rep(1,length(Parents))
    }
    DAG[causalOrder[i-1],node]<-1
  }
  return(list(DAG=DAG,TO=causalOrder))
}
# Hub-and-chain graph
randomDAG2_hub <- function(p,probConnect)
{
  # This function is modified from randomDAG2 function by Jonas Peters
  DAG <- diag(rep(0,p))
  causalOrder <- sample(p)
  Z<-10
  for(i in 1:(p))
  {
    node <- causalOrder[i]
    DAG[causalOrder[i-1],node]<-1
    if (i>2){DAG[causalOrder[sample(seq(min(i-1,Z)),2)],node]<-1}
  }
  DAG[causalOrder[1],causalOrder[2]] <- 1
  return(list(DAG=DAG,TO=causalOrder))
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

do_one<-function(p,n,pc){
  result<-matrix(0,6,2)
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
  return(result)
}

sim_run<-function(p,pc,n,N){
  # p:  number of nodes
  # pc: probability of connection
  # n:  number of samples
  # N:  number of simulation runs
  pb <- txtProgressBar(min = 0, max = N, style = 3)
  result<-matrix(0,6,2)
  results <-
    lapply(1:N, function(i){
      setTxtProgressBar(pb, i)
      do_one(p,n,pc)
    })
  close(pb)

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
  result<-sim_run(m,3/(2*m-2),n,N);
  #sink(paste("WarmGDS_ChainDAG_Simresults_p=",m,"n=",n,"Bmin",Bmin,"_sparse.txt"))
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
  # return dense result
  cat('Setting Dense: m=', m, " n=",n, "N=",N,"Bmin",Bmin, "\n")
  result<-sim_run(m,0.3,n,N);
  #sink(paste("WarmGDS_ChainDAG_Simresults_p=",m,"n=",n,"Bmin",Bmin,"_dense.txt"))
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
}
cat("done.\n")
