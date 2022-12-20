# Functions for learning locally metrical Huesler--Reiss graphical models

#Function computes the matrix A as described in Section 6 of the paper. 
Amatrix = function(graph){
  m = vcount(graph)
  triangles = as.vector(igraph::triangles(graph))
  A = matrix(0, nrow = length(triangles),ncol = choose(m,2))
  edges = combn(m,2,simplify = TRUE)
  for (i in (3*(0:(length(triangles)/3-1))+1) ) {
    A[c(i,i+1,i+2),which(edges[1,] %in% triangles[i:(i+2)] & edges[2,]%in% triangles[i:(i+2)])] = cbind(c(1,-1,-1),c(-1,1,-1),c(-1,-1,1))
  }
  return(t(A))
}

#Amatrix fast direct multiplication and indexing for algorithm
Atimesmu = function(mu,graph){
  edges = t(ends(graph,E(graph)))
  a = rep(0,ncol(edges))
  tri=as.vector(igraph::triangles(graph))
  for (i in (3*(0:(length(tri)/3-1))+1) ) {
    index = which(edges[1,] %in% tri[i:(i+2)] & edges[2,]%in% tri[i:(i+2)])
    a[index]=a[index]+cbind(c(1,-1,-1),c(-1,1,-1),c(-1,-1,1))%*%mu[c(i,i+1,i+2)]
  }
  return(a)
}

Atimesmu_index = function(graph){
  edges = t(ends(graph,E(graph)))
  index = list()
  tri=as.vector(igraph::triangles(graph))
  for (i in (3*(0:(length(tri)/3-1))+1) ) {
    index[[i]] = which(edges[1,] %in% tri[i:(i+2)] & edges[2,]%in% tri[i:(i+2)])
  }
  return(index)
}


#Function computes the matrix A without zero rows. 
smallAmatrix = function(graph){
  m = vcount(graph)
  triangles = as.vector(igraph::triangles(graph))
  A = matrix(0, nrow = length(triangles),ncol = ecount(graph))
  edges = t(ends(graph,E(graph)))
  for (i in (3*(0:(length(triangles)/3-1))+1) ) {
    A[c(i,i+1,i+2),which(edges[1,] %in% triangles[i:(i+2)] & edges[2,]%in% triangles[i:(i+2)])] = cbind(c(1,-1,-1),c(-1,1,-1),c(-1,-1,1))
  }
  return(t(A))
}

#Fast multiplication for transposed Amatrix
tAtimesGamma = function(Gamma,graph){
  tri=as.vector(igraph::triangles(graph))
  edges = t(ends(graph,E(graph)))
  Gam= Gamma[t(edges)]
  b=rep(0,length(tri))
  for (i in (3*(0:(length(tri)/3-1))+1) ) {
    b[i:(i+2)]=cbind(c(1,-1,-1),c(-1,1,-1),c(-1,-1,1))%*%Gam[which(edges[1,] %in% tri[i:(i+2)] & edges[2,]%in% tri[i:(i+2)])]
  }
  return(b)
}

#Simple check for local metric property
loc_metr = function(Gamma,graph){
  tri = as.vector(igraph::triangles(graph))
  m=0
  lm = FALSE
  failed=c()
  for (i in (3*(0:(length(tri)/3-1))+1) ) {
    if(Gamma[tri[i],tri[i+1]]<=
       Gamma[tri[i],tri[i+2]]+Gamma[tri[i+1],tri[i+2]]){m=m+1}else{failed=c(failed,i)}
    if(Gamma[tri[i],tri[i+2]]<=
       Gamma[tri[i],tri[i+1]]+Gamma[tri[i+1],tri[i+2]]){m=m+1}else{failed=c(failed,i)}
    if(Gamma[tri[i+1],tri[i+2]]<=
       Gamma[tri[i],tri[i+1]]+Gamma[tri[i],tri[i+2]]){m=m+1}else{failed=c(failed,i)}
  }
  if(length(failed)==0){lm = TRUE}
  return(list(lm=lm,failed_ineq=failed,number_valid_ineq=m))  
}

#Function returns product of matrix A and vectorized Gamma. Negative values mean that the triangle inequality holds.
local_metric_property = function(Gamma,graph){ 
  ineq = t(Amatrix(graph)) %*%Gamma2vec(Gamma)
  return(ineq)
}

#Transformation from Theta to Q
Theta2Qvec <-function(Theta){
  m = ncol(Theta)
  Q <- vector(length = choose(m,2))
  Q<- -Theta[lower.tri(Theta, diag=FALSE)]  
  return(Q)
}

#Function to vectorize Gamma
Gamma2vec <- function(Gamma){
  m = ncol(Gamma)
  G <- vector(length = choose(m,2))
  G<- Gamma[lower.tri(Gamma, diag=FALSE)]  
  return(G)
}

#Function to write Gamma vector in matrix
Gvec2Gamma <-function(Gvec){
  n=length(Gvec)
  m=sqrt(2*n+1/4)+1/2
  G = matrix(0,nrow = m,ncol = m)
  G[lower.tri(G, diag=FALSE)] <- Gvec
  return(G+t(G))
}

#Transformation from vectorized Q to Theta
Qvec2Theta <- function(Qvec){
  n=length(Qvec)
  m=sqrt(2*n+1/4)+1/2
  Theta = matrix(0,nrow = m,ncol = m)
  Theta[lower.tri(Theta, diag=FALSE)]=-Qvec
  Theta = Theta +t(Theta)
  diag(Theta)=-colSums(Theta)
  return(Theta)
}

#Transformation from vectorized Q to Theta^(m)
Qvec2Thetawithm <- function(Qvec,m){
  Theta = matrix(0,nrow = m,ncol = m)
  Theta[lower.tri(Theta, diag=FALSE)]=-Qvec
  Theta = Theta +t(Theta)
  diag(Theta)=-colSums(Theta)
  return(Theta)
}

