# Copyright (c) 2018 - 2019  Wenyu Chen [wenyuc@uw.edu]
# All rights reserved.  See the file COPYING for license terms.
###############
### Fit DAG using topological ordering
###############
#' Infer  DAG using topological ordering
#' @param X: data in n x p matrix
#' @param TO: topological ordering
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
#' @param threshold: only used in rls and chol. the hard threshold level.
#' @param FCD: only used in debiased lasso,  the FCD procedure [False Discovery Rate Control via Debiased Lasso. Javanmard and Montanari. 2018]
#' or use individual tests to select support.
#' @param precmtd: only used in debiased lasso, how to compute debiasing matrix
#'               "cv": node-wise lasso w/ joint 10 fold cv
#'               "sqrtlasso": square-root lasso (no tune, default)
#' @return Adjacency matrix with ADJ[i,j]!=0 iff i->j
DAG_from_Ordering<-function(X,TO,mtd="ztest",alpha=0.05,
                            threshold=1e-1,FCD=NULL,precmtd=NULL){
  n=dim(X)[1]
  p=dim(X)[2]
  if (p!=length(TO)){stop("length mismatch")}
  if (mtd=="ztest"){
    # sidak
    C=cor(X)
    adj=matrix(0,p,p)
    for (i in 2:p){
      u=TO[i]
      for (j in 1:(i-1)){
        v = TO[j]
        s = setdiff(TO[seq(i-1)],v)
        pval = 1-(2*pnorm(abs(pcalg::zStat(u,v,s,C=C,n=n)))-1)^(p*(p-1)/2)
        adj[v,u]=ifelse(pval<alpha,1,0)
      }
    }
    return(adj!=0)
  }
  if (mtd=="chol"){
    Sigma=cov(X)
    B = solve(chol(Sigma[TO,TO])[order(TO),order(TO)])
    gm = diag(p)-B%*%diag(1/diag(B))
    return(gm*(abs(gm)>threshold)!=0)
  }
  if (mtd=="rls"){
    gm = upper.tri(matrix(0,p,p))[order(TO),order(TO)]
    colnames(gm)=rownames(gm)=colnames(X)
    return(abs(t(ggm::fitDag(gm,cov(X),dim(X)[1])$Ahat))-diag(p)>threshold)
  } else {
    # dblasso
    if (is.null(FCD)){FCD="T"}
    if (is.null(precmtd)){precmtd="sqrtlasso"}
    gm = matrix(0,p,p)
    gm[TO[1],TO[2]]=anova(lm(X[,TO[2]]~X[,TO[1]]))$`Pr(>F)`[1]<alpha
    if(p==2){return(gm)}
    for (i in 3:p){
      gm[TO[1:(i-1)],TO[i]]=
        vselect(X[,TO[1:i-1]],X[,TO[i]],alpha=alpha,p_total = p,
                selmtd = mtd,FCD = FCD,precmtd = precmtd)$selected
    }
    return(gm!=0)
  }
}
