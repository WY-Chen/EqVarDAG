# Copyright (c) 2018 - 2019  Wenyu Chen [wenyuc@uw.edu]
# All rights reserved.  See the file COPYING for license terms.



###############
### Main method with top-down approach
###############

#' Estimate topological ordering and DAG using Top-down approach
#'
#' @param X An n-by-p data matrix.
#' @return Estimated Adjacency matrix and topological ordering.
#' @examples
#' X1<-rnorm(100)
#' X2<-X1+rnorm(100)
#' EqVarDAG_TD(cbind(X1,X2))
#'
#' #$adj
#' #[,1] [,2]
#' #[1,]    0    1
#' #[2,]    0    0
#' #
#' #$TO
#' #[1] 1 2
EqVarDAG_TD<-function(X){
  # Input
  # X: n by p matrix of data
  # Output
  # adj: estimated adjacency matrix
  # TO : estimated topological ordering
  n<-dim(X)[1]
  p<-dim(X)[2]
  rr<-rev(EqVarDAG_TD_internal(X)$TO)
  result<-matrix(0,p,p)
  for (ii in 1:(p-1)){
    now<-rr[ii]
    this<-sort(rr[(ii+1):p])
    if (length(this)>1){
      # variable selection
      lassom<-glmnet::cv.glmnet(X[,this],X[,now] )
      bfit<-coefficients(lassom)[-1]
      for (jj in 1:length(this)){
        if(bfit[jj]!=0)
          result[this[jj],now]<-1
      }
    } else {
      # deal with the last two nodes
      lmod<-summary(RcppEigen::fastLm(X[,now]~X[,this]))
      if (lmod$coef[2,4]<0.05) {
        result[this,now]<-1
      }
    }
  }
  return(list(adj=result,TO=rev(rr)))
}

###############
### helper functions
###############

# Compute topological ordering by conditioning
EqVarDAG_TD_internal<-function(X){
  n<-dim(X)[1]
  p<-dim(X)[2]
  done<-NULL
  done<-p+1
  S<-cov(X)
  Sinv<-solve(S)
  for(i in 1:p){
    varmap<-seq(p)[-done]
    v<-which.min(diag(solve(Sinv[-done,-done])))[1]
    done<-c(done,varmap[v])
  }
  return(list(TO=done[-1],support=NULL))
}
