## This algorithm solves the dual problem of step (2) given a step (1) solution and a graph.
loc_metr_alg <- function(Gamma,graph,max_iter=1000,print_level=0,algorithm="NLOPT_LD_MMA"){
if(length(triangles(graph))==0){
  if(print_level!=0){cat("There are no triangles in the given graph, the algorithm returns the input. ")}
  return(list(mu_opt=c(),Gamma_opt=Gamma))}
else{
lm=loc_metr(Gamma,graph)  
if(lm$lm==TRUE){
  if(print_level!=0){cat("Gamma already satisfies the local metric property with respect to the given graph. ")}
  return(list(mu_opt=rep(0,length(lm$number_valid_ineq)),Gamma_opt=Gamma))
}
else{
m=ncol(Gamma)
Theta_hat=Gamma2Theta(Gamma)
Qmatrix=-Theta_hat
diag(Qmatrix)=rep(0,m)
edges=ends(graph,E(graph))
tri=as.vector(igraph::triangles(graph))
d=length(tri)
index = Atimesmu_index(graph)
Amulti= function(mu){
  a = rep(0,nrow(edges))
  for (i in (3*(0:(length(tri)/3-1))+1) ) {
    a[index[[i]]]=a[index[[i]]]+cbind(c(1,-1,-1),c(-1,1,-1),c(-1,-1,1))%*%mu[c(i,i+1,i+2)]
  }
  return(a)
}
tAmulti=function(Gam){
  G= Gam[edges]
  b=rep(0,length(tri))
  for (i in (3*(0:(length(tri)/3-1))+1) ) {
    b[i:(i+2)]=cbind(c(1,-1,-1),c(-1,1,-1),c(-1,-1,1))%*%G[index[[i]]]
  }
  return(b)
}

project_psd <- function(X, tol = 1e-13) {
  X <- 0.5 * (X + t(X))                
  ev <- eigen(X, symmetric = TRUE)
  vals <- pmax(ev$values, tol)        
  ev$vectors %*% diag(vals) %*% t(ev$vectors)
}
obj = function(mu){
      Q_upper=Qmatrix
      Q_upper[edges]=Qmatrix[edges]+Amulti(mu)
      Q=t(Q_upper)[lower.tri(Q_upper,diag=FALSE)]
      Theta=project_psd(Qvec2Theta(Q))
     # return(ifelse(is.positive.definite(Theta[-m,-m]),-log(determinant(x=Theta[-m,-m] ,logarithm = FALSE)$modulus),Inf))
      return(-log(determinant(x=Theta[-m,-m] ,logarithm = FALSE)$modulus))
      
}
grad = function(mu){
    Q_upper=Qmatrix
    Q_upper[edges]=Qmatrix[edges]+Amulti(mu)
    Q=t(Q_upper)[lower.tri(Q_upper,diag=FALSE)]
    Theta=project_psd(Qvec2Theta(Q))
    Gamma=Theta2Gamma(Theta)
    grad = -tAmulti(Gamma)
    return(grad)
}
mu_opt<- nloptr(x0 = rep(0,d) , eval_f = obj, eval_grad_f = grad,
                lb=rep(0,d),ub=rep(Inf,d),
                      opts = list("algorithm"=algorithm, xtol_rel=1e-12,
                                         maxeval=max_iter,print_level=print_level))
Qmatrix[edges]=Qmatrix[edges]+Amulti(mu_opt$solution)
Q_opt=t(Qmatrix)[lower.tri(Qmatrix,diag=FALSE)]
return(list(mu_opt=mu_opt$solution,Gamma_opt=Theta2Gamma(project_psd(Qvec2Theta(Q_opt)))))
}}}
