#This script contains functions that we used to measure polarization in data methods of archetypal analysis were used,
#therefore library archetypes is necessary to get inputs in the required format

library(conclust)
#We need function that performs constrained K-means

#This function evaluates coefficient of determination for a model of archetypal analysis. This is needed in order to
#determine suitable number of archetypes for multiple datasets.
#Arguments:
  #X: data (observations as rows) - matrix
  #archet: archetypes object (output of the function archetypes from the library archetypes)
archet_R_squared = function(X, archet){
  mu = apply(X, MARGIN = 2, FUN = mean) #sample mean
  X_hat = archet$alphas %*% archet$archetypes #data aproximations
  TSS = sum(apply(t(apply(X, MARGIN = 1,
                          FUN = function(x) x - mu)), MARGIN = 1, FUN = norm, type = "2")^2) #total sum of squares
  RSS = sum(apply(X - X_hat, MARGIN = 1, FUN = norm, type = "2")^2) #residual sum of squares
  R_squared = 1 - RSS/TSS #coefficient of determination
  return(R_squared)
}

################################################################################

#This function evaluates coefficient of determination for a K-means model.
#Arguments
  #X: data (observations as rows) - matrix
  #clust: clustering (vector where the ith element describes clustering membership of the ith observation) - array
cluster_R_squared = function(X, clust){
    k = length(unique(clust)) #number of clusters
    mu = apply(X, MARGIN = 2, FUN = mean) #sample mean
    mu_clust = list() #this list will contain means of the clusters
    WSS_clust = rep(0, k) #this vector will contain contributions of clusters to WSS
    
    clust_i = list() #ith element of this list will be a vector listing indices of all the observations belonging to the
                     #ith cluster
    #in this cycle, we evaluate contribution of each cluster to WSS
    for (i in 1:k){
      clust_i[[i]] = which(clust == i)
      if (length(clust_i[[i]]) == 1){
        mu_clust[[i]] = X[clust_i[[i]], ]
        WSS_clust[i] = 0
      } else{
        mu_clust[[i]] = apply(X[clust_i[[i]], ], MARGIN = 2, FUN = mean)
        WSS_clust[i] = sum(apply(t(apply(X[clust_i[[i]], ], MARGIN = 1,
                                       FUN = function(x) x - mu_clust[[i]])), MARGIN = 1, FUN = norm, type = "2")^2)
      }
    }
    TSS = sum(apply(t(apply(X, MARGIN = 1,
                            FUN = function(x) x - mu)), MARGIN = 1, FUN = norm, type = "2")^2) #total sum of squares
    WSS = sum(WSS_clust)
    R_squared = 1 - WSS/TSS #coefficient of determination
    return(R_squared)
}

################################################################################

#This function evaluates coefficient of determination for a K-means model. Unlike the previous one, it uses weighted versions
#of all the sums of squares. These weigh sumands in the sums but don't alter the means (of clusters or the entire dataset).
#Arguments
  #X: data (observations as rows) - matrix
  #clust: clustering (vector where the ith element describes clustering membership of the ith observation) - array
  #w: weights of observations - diagonal matrix or an array
cluster_R_squared_w = function(X, clust, w){
  #we will use w in the form of a diagonal matrix
  if (class(w)[1] != "matrix"){
    w = diag(w)
  }
  k = length(unique(clust)) #number of clusters
  mu = apply(X, MARGIN = 2, FUN = mean) #sample mean
  mu_clust = list() #this list will contain means of the clusters
  WSS_clust = rep(0, k) #this vector will contain contributions of clusters to WSS
  
  clust_i = list() #ith element of this list will be a vector listing indices of all the observations belonging to the
                   #ith cluster
  #in this cycle, we evaluate contribution of each cluster to WSS
  for (i in 1:k){
    clust_i[[i]] = which(clust == i)
    if (length(clust_i[[i]]) == 1){
      mu_clust[[i]] = X[clust_i[[i]], ]
      WSS_clust[i] = 0
    } else{
      mu_clust[[i]] = apply(X[clust_i[[i]], ], MARGIN = 2, FUN = mean)
      w_col = apply(w[clust_i[[i]], ], MARGIN = 2, FUN = function(x) length(unique(x))) == 2
      WSS_clust[i] = sum(w[clust_i[[i]], w_col] %*% apply(t(apply(X[clust_i[[i]], ], MARGIN = 1,
                                            FUN = function(x) x - mu_clust[[i]])), MARGIN = 1, FUN = norm, type = "2")^2)
    }
  }
  TSS = sum(w %*% apply(t(apply(X, MARGIN = 1,
                                   FUN = function(x) x - mu)), MARGIN = 1, FUN = norm, type = "2")^2) #total sum of squares
  WSS = sum(WSS_clust)
  R_squared = 1 - WSS/TSS #coefficient of determination
  return(R_squared)
}

################################################################################

#This function evaluates coefficient of determination for a K-means model. Unlike the previous one, it uses weighted versions
#of all the sums of squares. These weigh sumands in the sums but also alter the means (of clusters or the entire dataset),
#replacing them with their weighted counterparts.
#Arguments
  #X: data (observations as rows) - matrix
  #clust: clustering (vector where the ith element describes clustering membership of the ith observation) - array
  #w: weights of observations - diagonal matrix or an array
cluster_R_squared_w_2 = function(X, clust, w){
  #we will use w in the form of a diagonal matrix
  if (class(w)[1] != "matrix"){
    w = diag(w)
  }
  k = length(unique(clust)) #number of clusters
  mu = apply(X, MARGIN = 2, FUN = function(x) weighted.mean(x, diag(w))) #weighted sample mean
  mu_clust = list() #this list will contain means of the clusters
  WSS_clust = rep(0, k) #this vector will contain contributions of clusters to WSS
  
  clust_i = list() #ith element of this list will be a vector listing indices of all the observations belonging to the
                   #ith cluster
  #in this cycle, we evaluate contribution of each cluster to WSS
  for (i in 1:k){
    clust_i[[i]] = which(clust == i)
    if (length(clust_i[[i]]) == 1){
      mu_clust[[i]] = X[clust_i[[i]], ]
      WSS_clust[i] = 0
    } else{
      w_col = apply(w[clust_i[[i]], ], MARGIN = 2, FUN = function(x) length(unique(x))) == 2
      mu_clust[[i]] = apply(X[clust_i[[i]], ], MARGIN = 2,
                            FUN = function(x) weighted.mean(x, diag(w[w_col, w_col])))
      WSS_clust[i] = sum(w[clust_i[[i]], w_col] %*% apply(t(apply(X[clust_i[[i]], ], MARGIN = 1,
                                               FUN = function(x) x - mu_clust[[i]])), MARGIN = 1, FUN = norm, type = "2")^2)
    }
  }
  TSS = sum(w %*% apply(t(apply(X, MARGIN = 1,
                                   FUN = function(x) x - mu)), MARGIN = 1, FUN = norm, type = "2")^2) #total sum of squares
  WSS = sum(WSS_clust)
  R_squared = 1 - WSS/TSS #coefficient of determination
  return(R_squared)
}

################################################################################

#This function finds clustering using weighted constrained K-means. It uses the source code of the function ckmeans from the
#library conclust. We only made small alterations which we will highlight by comments.
#Arguments: Identical to the arguments of the function ckmeans with the addition of:
  #w: weights of observations - diagonal matrix or an array
ckmeans_w = function (data, k, mustLink, cantLink, w, maxIter = 100) {
  #we will use w in the form of a vector
  if (class(w)[1] == "matrix"){
    w = diag(w)
  }
  w = w/sum(w) #normalization
  mu = apply(data, MARGIN = 2, FUN = sum)/nrow(data)
  data = t(apply(data, MARGIN = 1, FUN = function(x) x - mu)) #data centering
  
  dist <- function(x, y) {
    tmp <- x - y
    sum(tmp * tmp)
  }
  violate <- function(i, j) {
    for (u in mlw[[i]]) {
      if (label[u] != 0 && label[u] != j) 
        return(1)
    }
    for (u in clw[[i]]) {
      if (label[u] == j) 
        return(1)
    }
    0
  }
  findMustLink <- function(i) {
    tmp = c()
    for (j in 1:n) {
      if (M[i, j] == 1) 
        tmp = c(tmp, j)
    }
    tmp
  }
  findCantLink <- function(i) {
    tmp = c()
    for (j in 1:n) {
      if (C[i, j] == 1) 
        tmp = c(tmp, j)
    }
    tmp
  }
  data = as.matrix(data)
  n <- nrow(data)
  d <- ncol(data)
  nm <- nrow(mustLink)
  nc <- nrow(cantLink)
  M = matrix(0, nrow = n, ncol = n)
  for (i in 1:nm) {
    if (i > nm) 
      break
    u = mustLink[i, 1]
    v = mustLink[i, 2]
    M[u, v] = 1
    M[v, u] = 1
  }
  for (u in 1:n) {
    for (i in 1:n) {
      for (j in 1:n) {
        if (M[i, u] == 1 && M[u, j] == 1) 
          M[i, j] = 1
      }
    }
  }
  tp = rep(0, n)
  ntp = 0
  for (i in 1:n) {
    if (tp[i] == 0) {
      ntp = ntp + 1
      tp[i] = ntp
      j = i + 1
      while (j <= n) {
        if (tp[j] == 0 && M[i, j] == 1) 
          tp[j] = ntp
        j = j + 1
      }
    }
  }
  findMember <- function(v) {
    tmp = c()
    for (u in 1:n) {
      if (tp[u] == v) 
        tmp = c(tmp, u)
    }
    tmp
  }
  tmpPi = lapply(1:ntp, findMember)
  C = matrix(0, nrow = n, ncol = n)
  for (i in 1:nc) {
    if (i > nc) 
      break
    u = cantLink[i, 1]
    v = cantLink[i, 2]
    x = tp[u]
    y = tp[v]
    if (x != y) {
      for (p in tmpPi[[x]]) {
        for (q in tmpPi[[y]]) {
          C[p, q] = 1
          C[q, p] = 1
        }
      }
    }
  }
  mlw <- lapply(1:n, findMustLink)
  clw <- lapply(1:n, findCantLink)
  tmp <- sample(1:n, k)
  C <- matrix(nrow = k, ncol = d)
  for (i in 1:k) {
    C[i, ] = data[tmp[i], ]
  }
  for (iter in 1:maxIter) {
    label <- rep(0, n)
    for (i in 1:n) {
      dd <- rep(1e+15, k)
      best <- -1
      for (j in 1:k) {
        if (violate(i, j) == 0) {
          dd[j] <- dist(data[i, ], C[j, ])
          if (best == -1 || dd[j] < dd[best]) {
            best = j
          }
        }
      }
      if (best == -1) 
        return(0)
      label[i] <- best
    }
    if (iter == maxIter) 
      return(label)
    C2 <- matrix(0, nrow = k, ncol = d)
    dem <- rep(0, k)
    for (i in 1:n) {
      j = label[i]
      C2[j, ] = C2[j, ] + w[i]*data[i, ] #here, we compute weighted means instead of the standard ones
      dem[j] = dem[j] + w[i]
    }
    for (i in 1:k) {
      if (dem[i] > 0) 
        C[i, ] = 1 * C2[i, ]/dem[i]
    }
  }
}

################################################################################

#Multistart version of the previous function. It uses the function ckmeans_w so that one needs to be in the workspace.
#This function uses one additional argument:
  #nstart: number of runs - integer
ckmeans_w_multistart = function(data, k, mustLink, cantLink, vahy, maxIter = 100, nstart = 10){
  #we will use w in both forms
  if (class(w)[1] == "matrix"){
    w = diag(w)
  }
  w_mat = diag(w)
  mu = apply(data, MARGIN = 2, FUN = sum)/nrow(data)
  data = t(apply(data, MARGIN = 1, FUN = function(x) x - mu)) #normalization
  
  WSS = rep(0, nstart) #we will evaluate weighted version of WSS for every run
  clust = matrix(0, ncol = nrow(data), nrow = nstart) #matrix with clusterings for runs as rows
  #nstart runs
  for (j in 1:nstart){
    clust[j, ] = ckmeans_w(data, k, mustLink, cantLink, w, maxIter)
    mu_clust = list()
    WSS_clust = rep(0, k)
    clust_i = list()
    for (i in 1:k){
      clust_i[[i]] = which(clust[j, ] == i)
      if (length(clust_i[[i]]) == 1){
        mu_clust[[i]] = data[clust_i[[i]], ]
        WSS_clust[i] = 0
      } else{
        w_col = apply(w_mat[clust_i[[i]], ], MARGIN = 2, FUN = function(x) length(unique(x))) == 2
        mu_clust[[i]] = apply(data[clust_i[[i]], ], MARGIN = 2,
                              FUN = function(x) weighted.mean(x, diag(w_mat[w_col, w_col])))
        WSS_clust[i] = sum(w_mat[clust_i[[i]], w_col] %*% apply(t(apply(data[clust_i[[i]], ], MARGIN = 1,
                                                FUN = function(x) x - mu_clust[[i]])), MARGIN = 1, FUN = norm, type = "2")^2)
      }
    }
    WSS[j] = sum(WSS_clust)
  }
  best = which(WSS == min(WSS))[1] #best clustering is selected based on the weighted version of WSS
  return(clust[best, ])
}