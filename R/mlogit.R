# ====================================================== #
# Mixed-effects logistric regression with a random slope #
# ====================================================== #

# Assume a nonpararametric distribution for the random slope

# k     number of groups
# gi    number of observations in each group
# ni    number of binomial trials for each observation
# pt    support points
# pr    proportions
# beta  
# x     given covariates

rmlogit = function(k, gi=2, ni=2, pt=0, pr=1, beta=1, X) {
  pr = rep(pr, len=length(pt))
  theta = sample(pt, k, replace=TRUE, prob=pr)
  gi = rep(gi, len=k)
  gic = c(0,cumsum(gi))
  n = gic[k+1]
  ni = rep(ni, len=n)
  r = length(beta)
  data = matrix(nrow=n, ncol=r+3)
  dimnames(data) = list(NULL, c("group","yi","ni",paste("x",1:r,sep="")))
  for(i in 1:k) {
    j = (gic[i]+1):gic[i+1]
    data[j,1] = i
    if(missing(X)) xi = matrix(rnorm(r*gi[i])-.5, nrow=gi[i])
    else xi = X[j,,drop=FALSE]
    data[j,4:(r+3)] = xi
    eta = drop(theta[i] + xi %*% beta)
    p = 1 / (1 + exp(-eta))
    data[j,2] = rbinom(rep(1,gi[i]), ni[j], prob=p)
    data[j,3] = ni[j]
  }
  class(data) = "mlogit"
  attr(data, "ui") = 1:k
  attr(data, "gi") = gi
  data
}

mlogit = function(x) {
  x1s = sort(x[,1])
  attr(x, "ui") = unique(x1s)
  attr(x, "gi") = tabulate(x1s)
  k = ncol(x) - 3
  dimnames(x) = list(NULL,c("group","yi","ni",paste("x",1:k,sep="")))
  class(x) = "mlogit"
  x
}

length.mlogit = function(x) length(attr(x, "ui"))

weights.mlogit = function(x, beta) 1

suppspace.mlogit = function(x, beta) c(-Inf,Inf)

gridpoints.mlogit = function(x, beta, grid=100) {
  xbeta = drop(x[,-(1:3),drop=FALSE] %*% beta)
  ui = attr(x, "ui")
  xbeta.max = xbeta.min = ybar = double(length(ui))
  for(i in 1:length(ui)) {
    xbeta.max[i] = max(xbeta[x[,1]==ui[i]])
    xbeta.min[i] = min(xbeta[x[,1]==ui[i]])
    ybar[i] = sum(x[x[,1]==ui[i],2]) / sum(x[x[,1]==ui[i],3])
  }
  yn = pmin(pmax(ybar, 0.0001), 0.9999)
  seq(min(log(yn/(1-yn))-xbeta.max), max(log(yn/(1-yn))-xbeta.min), len=grid)
}

initial.mlogit = function(x, beta=NULL, mix=NULL, kmax=NULL) {
  if(is.null(beta)) {
    resp = cbind(x[,2], x[,3]-x[,2])
    coef = glm(resp ~ x[,-(1:3)], family=binomial)$coef
    names(coef) = NULL
    beta = coef[-1]
  }
  beta = rep(beta, len=ncol(x)-3)
  if(is.null(kmax)) kmax = 10
  if(is.null(mix) || is.null(mix$pt)) 
    mix = disc(gridpoints(x, beta, grid=kmax), mix$pr)
  list(beta=beta, mix=mix)
}

valid.mlogit = function(x, beta, theta) TRUE

logd.mlogit = function(x, beta, pt, which, eps=1e-10) {
  dl = vector("list", 3)
  names(dl) = c("ld","db","dt")
  xij = x[,-(1:3),drop=FALSE]
  xbeta = xij %*% beta
  dim(xbeta) = NULL
  eta = rep(xbeta, length(pt)) + rep(pt, each=length(xbeta))
  # p = pmin(pmax(1 / (1 + exp(-eta)), eps), 1-eps)   # so no NA is produced
  p = 1 / (1 + exp(-eta))
  p[p < eps] = eps
  p[p > 1-eps] = 1 - eps               # To avoid producing NA's
  dim(p) = c(nrow(x), length(pt))
  ui = attr(x, "ui")
  if(which[1] == 1) {
    a = x[,2] * log(p) + (x[,3]-x[,2]) * log(1-p)
    dl$ld = matrix(nrow=length(ui), ncol=ncol(a))
    for(i in 1:length(ui)) dl$ld[i,] = colSums(a[x[,1]==ui[i],,drop=FALSE])
  }
  if(which[2] == 1) {
    d = sweep(array(x[,2] - x[,3]*p, dim=c(nrow(p), length(pt), length(beta))),
      c(1,3), xij, "*")
    dl$db = array(dim=c(length(ui), length(pt), length(beta)))
    for(i in 1:length(ui))
      dl$db[i,,] = apply(d[x[,1]==ui[i],,,drop=FALSE], 2:3, sum)
  }
  if(which[3] == 1) {
    a = x[,2] - x[,3] * p
    dl$dt = matrix(nrow=length(ui), ncol=ncol(a))
    for(i in 1:length(ui))
      dl$dt[i,] = colSums(a[x[,1]==ui[i],,drop=FALSE])
  }
  dl
}

