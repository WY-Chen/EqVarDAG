# Copyright (c) 2018 - 2019  Wenyu Chen [wenyuc@uw.edu]
# All rights reserved.  See the file COPYING for license terms.



###############
### Main method with top-down approach
###############
#' Estimate topological ordering and DAG using high dimensional top-down approach
#'
#' @param X An n-by-p data matrix.
#' @param J maximum number of parents to condition on.
#' @return Estimated Adjacency matrix and topological ordering.
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
EqVarDAG_HD_TD<-function(X,J){
  # Input
  # X : n by p matrix of data
  # J : maximum number of parents to condition on
  # Output
  # adj: estimated adjacency matrix
  # TO : estimated topological ordering
  n<-dim(X)[1]
  p<-dim(X)[2]
  rr<-rev(getOrdering(X,J)) # use J=5
  result<-matrix(0,p,p)
  for (ii in 1:(p-1)){
    now<-rr[ii]
    this<-sort(rr[(ii+1):p])
    if (length(this)>1){
      # variable selection
      if (n>100){
        lassom<-glmnet::cv.glmnet(X[,this],X[,now]  )
        bfit<-coefficients(lassom)[-1]
      } else {
        lassom<-glmnet::glmnet(X[,this],X[,now] )
        bic<-n*log(colSums((predict(lassom,
                                    X[,this])-X[,now])^2)/n)+lassom$df*log(n)+
          2*lassom$df*log(p-ii)
        bfit<-coefficients(lassom)[,which(bic==min(bic))[1]][-1]
      }
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
