# Copyright (c) 2018 - 2020  Wenyu Chen [wenyuc@uw.edu]
# All rights reserved.  See the file COPYING for license terms.

###############
### Main method with high-dimensional bottom-up CLIME approach
###############
#' Estimate topological ordering and DAG using high dimensional bottom-up CLIME approach
#' Estimate  DAG using topological ordering
#' @param X,Y: n x p and 1 x p matrix
#' @param alpha: desired selection significance level
#' @param mtd: methods for learning DAG from topological orderings.
#'  "ztest": (p<n) [Multiple Testing and Error Control in Gaussian Graphical Model Selection. Drton and Perlman.2007]
#'  "rls": (p<n) fit recursive least squares using ggm package and threshold the regression coefs
#'  "chol": (p<n) perform cholesky decomposition and threshold the regression coefs
#'  "dlasso": debiased lasso (default with FCD=True and precmtd="sqrtlasso");
#'   "lasso": lasso with fixed lambda from [Penalized likelihood methods for estimation of sparse high-dimensional directed acyclic graphs. Shojaie and Michailidis. 2010];
#'   "adalasso": adaptive lasso with fixed lambda from [Shojaie and Michailidis. 2010];
#'   "cvlasso": cross-validated lasso from glmnet;
#'    "scallasso": scaled lasso.
#' @param threshold: for rls and chol, the threshold level.
#' @param FCD: for debiased lasso, use the FCD procedure [False Discovery Rate Control via Debiased Lasso. Javanmard and Montanari. 2018]
#' or use individual tests to select support.
#' @param precmtd: for debiased lasso, how to compute debiasing matrix
#'               "cv": node-wise lasso w/ joint 10 fold cv
#'               "sqrtlasso": square-root lasso (no tune, default)
#' @return Adjacency matrix with ADJ[i,j]!=0 iff i->j, and topological ordering
#' @examples
#' X1<-rnorm(100)
#' X2<-X1+rnorm(100)
#' EqVarDAG_HD_TD(cbind(X1,X2),2)
#'
#' #$adj
#' #[,1] [,2]
#' #[1,]    0    1
#' #[2,]    0    0
#' #
#' #$TO
#' #[1] 1 2
EqVarDAG_HD_CLIME<-function(X,mtd="dlasso",alpha=0.05,
                            threshold=1e-1,FCD=TRUE,precmtd="sqrtlasso"){
  # Input
  # X : n by p matrix of data
  # cv: if true, use cv-ed lambda, else use lambdafix,default True
  # lambdafix: customized lambda, default 0.1
  # Output
  # adj: estimated adjacency matrix
  # TO : estimated topological ordering
  n<-dim(X)[1]
  p<-dim(X)[2]
  TO=EqVarDAG_HD_CLIME_internal(X,NULL)
  adj=DAG_from_Ordering(X,TO,mtd,alpha,threshold,FCD,precmtd)
  return(list(adj=adj,TO=TO))
}


###############
### helper functions
###############
EqVarDAG_HD_CLIME_internal<-function(X,lam=NULL){
  # (i,j)=1 in fixedorder means i is ancestral to j
  n=dim(X)[1]
  p=dim(X)[2]
  Sigma=cov(X)
  if (is.null(lam)){lam = 4/sqrt(n)*sqrt(log(p/sqrt(0.05)))}
  Theta=clime_theta(X)
  TO=NULL
  while (length(TO)<p-1){
    sink=which.min(diag(Theta))
    TO=c(TO,sink)
    s = setdiff(which(Theta[,sink]!=0),sink)
    for (j in s){
      sj = unique(c(j,setdiff(which(Theta[,j]!=0),sink),s))
      if (length(sj)==1){
        Theta[j,sj]=Theta[sj,j]=1/Sigma[sj,sj]
      } else {
        Theta[j,sj]=Theta[sj,j]=clime_lp(Sigma[sj,sj],lam,1)
      }
    }
    Theta[,sink]=Theta[sink,]=rep(0,p)
    Theta[sink,sink]=Inf
  }
  return(rev(unname(c(TO,setdiff(seq(p),TO)))))
}

# clime utilities
unit_vec<-function(p,i){v=rep(0,p);v[i]=1;return(v)}
clime_lp<-function(Sigma,lam,j){
  # Sigma: cov(X)
  # lam: tuning parameter
  # j: the j-th problem
  p = dim(Sigma)[2]
  f.obj = rep(1,p*2) # sum (u+v), u=max(x,0), v=max(-x,0), x=u-v
  const.mat = rbind(
    cbind(Sigma,-Sigma), # Sigma*(u-v) >= lam +ej
    cbind(-Sigma,Sigma), # -Sigma*(u-v) >= lam-ej
    cbind(diag(p),matrix(0,p,p)), # u>0
    cbind(matrix(0,p,p),diag(p))  # v>0
  )
  const.dir = c(
    rep("<=",2*p),rep(">=",2*p)
  )
  const.rhs = c(
    rep(lam,p)+unit_vec(p,j),
    rep(lam,p)-unit_vec(p,j),
    rep(0,2*p)
  )
  lpout=lpSolve::lp(direction = "min",objective.in = f.obj,
                    const.mat = const.mat,const.dir = const.dir,
                    const.rhs = const.rhs)
  return(lpout$solution[1:p]-lpout$solution[(p+1):(2*p)])
}

clime_theta<-function(X,lam=NULL){
  p=dim(X)[2]
  n=dim(X)[1]
  Sigma = cov(X)
  if (is.null(lam)){lam = 4/sqrt(n)*sqrt(log(p/sqrt(0.05)))}
  Omega = sapply(1:p, function(i)clime_lp(Sigma,lam,i))
  Omega = (abs(Omega)<=abs(t(Omega)))*Omega+
    (abs(Omega)>abs(t(Omega)))*t(Omega)
  return(Omega)
}

cv_clime<-function(X){
  p=dim(X)[2]
  n=dim(X)[1]
  Sigma = cov(X)
  ws = sqrt(diag(Sigma))
  lams=exp(seq(log(1e-4),log(0.8),length.out = 100))
  ebics=sapply(1:100, function(j){
    Theta = sapply(1:p, function(i)clime_lp(Sigma,lams[j],i))
    Theta = (abs(Theta)<=abs(t(Theta)))*Theta+
      (abs(Theta)>abs(t(Theta)))*t(Theta)
    loglikGGM(S=Sigma,Theta=Theta)-sum(Theta!=0)/2*(log(p)+0.5*log(n))/n
  })
  return(clime_theta(X,lam = lams[which.max(ebics)]))
}



