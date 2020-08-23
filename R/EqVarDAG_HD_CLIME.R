# Copyright (c) 2018 - 2019  Wenyu Chen [wenyuc@uw.edu]
# All rights reserved.  See the file COPYING for license terms.

###############
### Main method with high-dimensional bottom-up CLIME approach
###############
#' Estimate topological ordering and DAG using high dimensional bottom-up CLIME approach
#'
#' @param X An n-by-p data matrix.
#' @param cv Obtain regularization parameter by cross-validation or use speficied value. Default True.
#' @param lambdafix Specified lambda. (If cv=F).
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
EqVarDAG_HD_CLIME<-function(X,cv=T,lambdafix=0.1){
  # Input
  # X : n by p matrix of data
  # cv: if true, use cv-ed lambda, else use lambdafix,default True
  # lambdafix: customized lambda, default 0.1
  # Output
  # adj: estimated adjacency matrix
  # TO : estimated topological ordering
  n<-dim(X)[1]
  p<-dim(X)[2]
  rr<-rev(EqVarDAG_HD_CLIME_internal(X,cv,lambdafix))
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
        bic<-n*log(colSums((predict(lassom,X[,this])-X[,now])^2)/n)+lassom$df*log(n)+
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
# estimate topological ordering
EqVarDAG_HD_CLIME_internal<-function(X,cv=T,lambdafix=0.1){
  p<-dim(X)[2]
  n<-dim(X)[1]
  S<-cov(X)
  Theta<-rep(0,p)
  lambda<- exp(-seq(0.2,7,length.out = 10))
  # 4-fold cv
  if (cv){
    cv.folds<-cvTools::cvFolds(n,4)
    lik<-rep(0,length(lambda))
    for (k in 1:4){
      x.test<-X[cv.folds$which==k,]
      x.fit <-X[cv.folds$which!=k,]
      S.fit <- cov(x.fit)
      S.test<- cov(x.test)
      out<-lapply(lambda,function(l){clime(S.fit,l)})
      lik<-lik+sapply(out, function(t)sum(diag(t%*%S.test))-
                        log(det(t))+
                        (sum(t!=0)-p)/2*(log(p)+log(n)))
    }
    Omega<-clime(S,lambda[which.min(lik)])
  } else {
    Omega<-clime(S,lambdafix)
  }
  # first variable
  Theta[1]<-which.min(diag(Omega))
  # other variables
  for (i in 2:(p-2)){
    supp_changed<-which(Omega[Theta[i-1],]!=0)
    supp_changed<-supp_changed[!supp_changed%in% Theta]
    if (cv){
      # using cv-ed CLIME
      lik<-NULL
      out<-lapply(lambda, function(l){
        Omega_temp<-Omega
        for (j in supp_changed){
          supp_j<-union(supp_changed,which(Omega_temp[j,]!=0))
          supp_j<-sort(supp_j[!supp_j%in% Theta])
          if (!gtools::invalid(supp_j)&length(supp_j)>1){
            w_j<-clime_lp(S[supp_j,supp_j],which(supp_j==j),l)
            Omega_temp[j,supp_j]<-w_j
            Omega_temp[supp_j,j]<-w_j
          }
        }
        Omega_temp[Theta[i-1],]<-0
        Omega_temp[,Theta[i-1]]<-0
        lik<-sum(diag(Omega_temp[-Theta[seq(i)],-Theta[seq(i)]]%*%
                        S[-Theta[seq(i)],-Theta[seq(i)]]-diag(p-i+1))^2)
        return(list(Omega_temp=Omega_temp,lik=lik))
      })
      Omega<-out[[which.min(sapply(out, function(x)x$lik))]]$Omega_temp
    } else {
      # Using fix-lambda CLIME
      for (j in supp_changed){
        supp_j<-union(supp_changed,which(Omega[j,]!=0))
        supp_j<-sort(supp_j[!supp_j%in% Theta])
        if (!gtools::invalid(supp_j)&length(supp_j)>1){
          w_j<-clime_lp(S[supp_j,supp_j],which(supp_j==j),lambdafix) # lambda
          Omega[j,supp_j]<-w_j
          Omega[supp_j,j]<-w_j
        } else if (length(supp_j)==1){
          Omega[j,j]<-1/S[j,j]
        }
      }
      Omega[Theta[i-1],]<-0
      Omega[,Theta[i-1]]<-0
    }
    Theta[i]<-setdiff(seq(p),Theta)[which.min(diag(Omega[-Theta[seq(i)],-Theta[seq(i)]]))]
  }
  # last two variables
  if (c(1,-1)%*%diag(S[setdiff(seq(p),Theta),setdiff(seq(p),Theta)])>0){
    Theta[c(p-1,p)]<-setdiff(seq(p),Theta)
  } else {
    Theta[c(p,p-1)]<-setdiff(seq(p),Theta)
  }
  return(rev(Theta))
}

# linear programming for rank-1 clime problems
clime_lp<-function(S,i,lambda){
  p<-dim(S)[1]
  ei<-rep(0,p); ei[i]<-1
  lsol<-Rglpk::Rglpk_solve_LP(obj = c(0,rep(0,p),rep(0,p),rep(1,p)),
                       mat = rbind(cbind(rep(-1,p),S,-S,matrix(0,p,p)),
                                   cbind(rep(-1,p),-S,S,matrix(0,p,p)),
                                   c(1,rep(0,p*3)),
                                   cbind(rep(0,p),diag(p),-diag(p),-diag(p)),
                                   cbind(rep(0,p),-diag(p),diag(p),-diag(p))),
                       dir = rep("<=",1+4*p),
                       rhs = c(ei,-ei,lambda,rep(0,p*2)),
                       control = list("tm_limit" = 500))
  return(lsol$solution[2:(p+1)]-lsol$solution[(p+2):(1+p*2)])
}

# clime solver
clime<-function(S,lambda){
  p<-dim(S)[1]
  Omega<-sapply(1:p, function(i)clime_lp(S,i,lambda))
  return(pmin(Omega,t(Omega)))
}



