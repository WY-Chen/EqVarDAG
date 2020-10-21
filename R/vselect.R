# Copyright (c) 2018 - 2019  Wenyu Chen [wenyuc@uw.edu]
# All rights reserved.  See the file COPYING for license terms.

###############
### vselect
###############
#' Variable selection toolbox for lasso,
#' @param X,Y: n x p and 1 x p matrix
#' @param alpha: desired selection significance level
#' @param p_total: total number of variables, used in bonferroni correction
#' @param selmtd: variable selection method: avaiable methods:
#'   "dlasso": debiased lasso (default with FCD=True and precmtd="sqrtlasso");
#'   "lasso": lasso with fixed lambda from [Penalized likelihood methods for estimation of sparse high-dimensional directed acyclic graphs. Shojaie and Michailidis. 2010];
#'   "adalasso": adaptive lasso with fixed lambda from [Shojaie and Michailidis. 2010];
#'   "cvlasso": cross-validated lasso from glmnet;
#'    "scallasso": scaled lasso.
#' @param FCD: for debiased lasso, use the FCD procedure [False Discovery Rate Control via Debiased Lasso. Javanmard and Montanari. 2018]
#' or use individual tests to select support.
#' @param precmtd: for debiased lasso, how to compute debiasing matrix
#'               "cv": node-wise lasso w/ joint 10 fold cv
#'               "sqrtlasso": square-root lasso (no tune, default)
#' @return one-hot model selection vector (1 x p) and fitted residuals


vselect<-function(X,Y,alpha=0.05,p_total=0,
                  selmtd="dlasso",
                  precmtd="sqrtlasso",
                  FCD=TRUE){
  p = dim(X)[2]
  n = dim(X)[1]
  # Debiased lasso
  if (selmtd=="dlasso"){
    # obtain debiased lasso estimate, recycle Mhat if possible
    dblt = debiased_lasso(X,Y,precmtd = precmtd,nlam=100)
    Ts = dblt$T_i
    if (FCD==FALSE){
      # select by individual tests
      return(list(selected = (2*(1-pnorm(abs(Ts)))<=alpha),
                  res = Y-X%*%dblt$coef))
    }
    # FCD procedure:
    t_p = sqrt(2*log(p)-2*log(log(p)))
    ts = sort(abs(Ts)[abs(Ts)<=t_p])
    ps = 2*p*(1-pnorm(ts))/
      pmax(1,sapply(ts, function(l)sum(abs(Ts)>l)))
    # compute FCD cutoff
    if (!length(which(ps<=(alpha)))){
      t0 <- sqrt(2*log(p))
    } else {
      t0 <- max(ts[which(ps<=(alpha))])
    }
    return(list(selected = (abs(Ts)>=t0),
                res = Y-X%*%dblt$coef))
  } else if (selmtd=="lasso"){
    # lasso with fixed penalty
    p_total=ifelse(p_total>0,p_total,p)
    lam = 1/sqrt(n)*qnorm(1-alpha/2/p_total/(p-1))
    lassom = coefficients(glmnet::glmnet(X,Y,lambda=lam))[-1]
    return(list(selected = (lassom!=0),res = Y-X%*%lassom))
  } else if (selmtd=="cvlasso"){
    # lasso with cv penalty
    lassom = coefficients(glmnet::cv.glmnet(X,Y))[-1]
    return(list(selected = (lassom!=0),res = Y-X%*%lassom))
  } else if (selmtd=="scallasso"){
    # scaledlasso
    lassom = coef(scalreg::scalreg(X,Y))
    return(list(selected = (lassom!=0),res = Y-X%*%lassom))
  } else if (selmtd=="adalasso"){
    # adaptive lasso with fixed penalty initialization
    p_total=ifelse(p_total>0,p_total,p)
    lam = 1/sqrt(n)*qnorm(1-alpha/2/p_total/(p-1))
    lassom = coefficients(glmnet::glmnet(X,Y,lambda=lam))[-1]
    if (all(lassom==0)){
      return(list(selected=rep(0,p), res=Y))
    } else {
      lassom = coefficients(
        glmnet::glmnet(X,Y,lambda=lam,
               penalty.factor =
                 pmax(1,1/abs(lassom)!=0)))[-1]
      return(list(selected = (lassom!=0),res = Y-X%*%lassom))
    }
  }
}

debiased_lasso<-function(X,Y,precmtd="cv",nlam=100){
  # Compute debiased lasso fit
  # using the algorithm from Van der Geer 2014
  # X: n x p Y: 1 x p
  # precmtd: chooses how to compute the debiasing matrix:
  #   "cv": joint cross-validated node-wise lasso
  #   "sqrtlasso": square-root lasso (no tune)
  # return coefficients, test statistics for variable selection
  p = dim(X)[2]
  n = dim(X)[1]
  # Compute debiasing matrix
  if (p>2){
    if (precmtd=="cv"){
      # lambda_max, implemented same as in glmnet
      lmax = max(sapply(1:p, function(i)
        max( abs(t(X[,i] - mean(X[,i])*
                     (1-mean(X[,i]))) %*% X[,-i] ) )/n))
      # lambda sequence
      lams =
        exp(seq(from = log(0.001*lmax),
                to = log(lmax),length.out = nlam))
      # CV for best lambda
      lamj = lams[which.min(rowSums(
        sapply(1:p, function(i)
          glmnet::cv.glmnet(X[,-i],X[,i],lambda=lams)$cvm)))]
      # matrix of coefs
      ghat = sapply(1:p, function(i)
        coefficients(glmnet::glmnet(X[,-i],X[,i],
                            lambda=lamj))[-1])
    } else if (precmtd == "sqrtlasso"){
      ghat = sapply(1:p, function(i)
        RPtests::sqrt_lasso(X[,-i],X[,i]))
    }
    # construct debiasing matrix as in VdG2014
    Chat = diag(rep(1,p))
    for (i in 1:p){Chat[i,-i]=-ghat[,i]}
    tauhat = sapply(1:p, function(i)
      t((X[,i]-X[,-i]%*%
           ghat[,i]))%*%X[,i]/n)
    Mhat = diag(1/tauhat)%*%Chat
  } else {
    Mhat = solve(cov(X))
  }

  # get lasso fit
  lmod= scalreg::scalreg(X,Y)
  # get estimated error noise level by refitting
  s = seq(p)[lmod$coefficients!=0]
  if (length(s)){
    theta_hat = sqrt(sum(resid(lm(Y~X[,s]))^2)/n)
  } else {
    theta_hat = sd(Y)
  }
  # debias
  dlfit = coefficients(lmod)+Mhat%*%t(X)%*%(Y-predict(lmod,X))/n
  # variable selection test stat
  T_i = sqrt(n)*dlfit/theta_hat/
    sqrt(diag(Mhat%*%cov(X)%*%t(Mhat)))
  return(list(coef=dlfit,T_i=T_i))
}

