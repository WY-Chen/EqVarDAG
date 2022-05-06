# Copyright (c) 2018 - 2020  Wenyu Chen [wenyuc@uw.edu]
# All rights reserved.  See the file COPYING for license terms.

###############
### Main method with bottom-up approach
###############
#' Estimate topological ordering and DAG using bottom-up approach
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
EqVarDAG_BU<-function(X,mtd = "dlasso",
  alpha = 0.05,
  threshold = NULL,
  FCD = TRUE,
  precmtd = "sqrtlasso"){
  # Input
  # X : n by p matrix of data
  # Output
  # adj: estimated adjacency matrix
  # TO : estimated topological ordering
  n<-dim(X)[1]
  p<-dim(X)[2]
  TO=EqVarDAG_BU_internal(X)
  adj=DAG_from_Ordering(X,TO,mtd,alpha,threshold,FCD,precmtd)
  return(list(adj=adj,TO=TO))
}

# compute precision matrix
get_precision<-function(X){
  n<-dim(X)[1]
  p<-dim(X)[2]
  S<-cov(X)
  if (n>p){
    return(solve(S))
  } else {
    message("Warning: EqVarDAG:BU n<p")
  }
}

###############
### helper functions
###############

# compute topological ordering
EqVarDAG_BU_internal<-function(X){
  n<-dim(X)[1]
  p<-dim(X)[2]
  # get precision matrix
  S<-cov(X)
  done<-p+1
  for(i in 1:p){
    varmap<-seq(p)[-done]
    v<-which.min(diag(solve(S[-done,-done])))[1]
    done<-c(done,varmap[v])
  }
  return(rev(done[-1]))
}
