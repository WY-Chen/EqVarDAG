# Copyright (c) 2018 - 2020  Wenyu Chen [wenyuc@uw.edu]
# All rights reserved.  See the file COPYING for license terms.

###############
### Main method with high-dimensional top-down approach (best subset reg)
###############
#' Estimate topological ordering and DAG using high dimensional top-down approach (best subset reg)
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
EqVarDAG_HD_TD<-function(X,J=3,mtd="dlasso",alpha=0.05,
                            threshold=1e-1,FCD=TRUE,precmtd="sqrtlasso"){
  # Input
  # X : n by p matrix of data
  # Output
  # adj: estimated adjacency matrix
  # TO : estimated topological ordering
  n<-dim(X)[1]
  p<-dim(X)[2]
  TO=getOrdering(X,J)
  adj=DAG_from_Ordering(X,TO,mtd,alpha,threshold,FCD,precmtd)
  return(list(adj=adj,TO=TO))
}

###############
### helper functions
###############
# compute best subset search
helper.func <- function(z, Y, Theta, J,mtd="exhaustive"){
  leaps::regsubsets(x = Y[, Theta, drop = F],
                    y = Y[, z, drop = F],
                    method=mtd,
                    nbest = 1,
                    nvmax = min(J, sum(Theta > 0 )), really.big = T)
}
# compute topological ordering
getOrdering <- function(Y, J){
  p <- dim(Y)[2]
  variances <- apply(Y, MAR = 2, sd)
  Theta <- rep(0, p)
  Theta[1] <- which.min(variances)
  out <- sapply(setdiff(1:p, Theta[1]),
                  function(z){
                    sum(resid(RcppEigen
                              ::fastLm(Y[, z] ~ Y[, Theta[1], drop = F]) )^2)})
  Theta[2] <- setdiff(1:p, Theta[1])[which.min(out)]
  for(i in 3:p){
    out <- lapply(setdiff(1:p, Theta),
                  function(jj)helper.func(jj, Y, Theta[seq(i-1)], J))
    nextRoot <- which.min(sapply(out,function(x){min(x$rss)}))
    Theta[i] <- setdiff(1:p, Theta)[nextRoot]
  }
  return(Theta)
}
