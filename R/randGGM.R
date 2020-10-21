############
### Useful fns
############
plot_matrix<-function(M){
  # plot matrix as heatmap
  levelplot(t(M[ncol(M):1,1:ncol(M)]),
            xaxt='n',yaxt='n',
            col.regions=c(0,1),colorkey=FALSE,
            scales=list(alternating=0))
}

uchol <-function(S){
  # S = U%*%t(U), U upper triangular, instead of t(U)%*%U
  p=dim(S)[1]
  t(chol(S[p:1,p:1]))[p:1,p:1]
}

############
### Generate data
############

rand.ER.dag = function(p, s, l,u) {
  wFUN <- function(m,lB,uB) {runif(m,lB,uB) }
  g = pcalg::randDAG(p, s, method = "er",wFUN=list(wFUN,l,u))
  gs=as(as(g,"matrix")[RBGL::tsort(g),RBGL::tsort(g)]!=0,"graphNEL")
  igs=igraph::igraph.from.graphNEL(gs)
  if (length(igraph::E(igs))){
    igraph::E(igs)$weight<-igraph::E(igraph::igraph.from.graphNEL(g))$weight
  }
  igraph::igraph.to.graphNEL(igs)
}

rand.PL.dag = function(p, s, l,u) {
  wFUN <- function(m,lB,uB) {runif(m,lB,uB) }
  g = pcalg::randDAG(p, s, method = "power",wFUN=list(wFUN,l,u))
  gs=as(as(g,"matrix")[RBGL::tsort(g),RBGL::tsort(g)]!=0,"graphNEL")
  igs=igraph::igraph.from.graphNEL(gs)
  if (length(igraph::E(igs))){
    igraph::E(igs)$weight<-igraph::E(igraph::igraph.from.graphNEL(g))$weight
  }
  igraph::igraph.to.graphNEL(igs)
}

randGGM = function(n=500,p=10,s=2,
                   eqvar=T,gm=NULL,
                   Omega=NULL,bmin=0.1,pl=F){
  if (!is.null(gm)){
    g = as(gm,"graphNEL") # use given autoregression matrix
  } else {
    if (!pl){
      g=rand.ER.dag(p, s, bmin, 1)
    } else {
      g=rand.PL.dag(p, s, bmin, 1)
    }
    gm = as(g,"matrix")
    colnames(gm)=rownames(gm)=as.character(seq(p))
    g = as(gm,"graphNEL")
  }
  if (!is.null(Omega)){
    # do nothing
    errs = mvtnorm::rmvnorm(n,sigma = Omega)
  } else if (eqvar){
    # eq var
    Omega = diag(p)
    errs = mvtnorm::rmvnorm(n,sigma = Omega)
  } else {
    Omega = diag(runif(p,1,3))
    errs = mvtnorm::rmvnorm(n,sigma = Omega)
  }
  X = rmvDAG(n,g,errMat = errs)
  Theta=(diag(p)-gm)%*%solve(Omega)%*%t(diag(p)-gm)
  Sigma=solve(t(diag(p)-gm))%*%Omega%*%solve(diag(p)-gm)
  return(list(X=X,gm=gm,g=g,Theta=Theta,Sigma=Sigma,Omega=Omega))
}
