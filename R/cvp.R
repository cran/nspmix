# =========================== #
# The common variance problem #
# =========================== #

# ========== #
# Class cvps #
# ========== #

# The common variance problem

rcvp = function(k, ni=2, mu=0, pr=1, sd=1) {
  if(length(mu) == 1) mu = rep(mu, k)
  else mu = sample(mu, k, replace=TRUE, prob=pr)
  ni = rep(ni, len=k)
  nicum = c(0,cumsum(ni))
  n = nicum[k+1]
  x = matrix(nrow=n, ncol=2)
  dimnames(x) = list(NULL, c("group", "x"))
  for(i in 1:k) {
    j = (nicum[i]+1):(nicum[i+1])
    x[j,1] = i
    x[j,2] = rnorm(ni[i], mu[i], sd)
  }
  x
}

# The common variance problem with data stored in sufficient statistics
# (ni, mi, ri) for each group i, denoting the number of observations, the
# mean and the residual sum of squares, respectively.

cvps = function(x) {
  ni = tapply(x[,2], x[,1], length)
  group = as.numeric(names(ni))
  ni = as.vector(ni)
  mi = as.vector(tapply(x[,2], x[,1], mean))
  ri = as.vector(tapply(x[,2], x[,1], var) * (ni-1))
  ri[is.na(ri)] = 0
  names(ni) = names(mi) = names(ri) = NULL
  xs = list(group=group, ni=ni, mi=mi, ri=ri)
  class(xs) = "cvps"
  xs
}

rcvps = function(k, ni=2, mu=0, pr=1, sd=1)
  cvps( rcvp(k=k, ni=ni, mu=mu, pr=pr, sd=sd) )

print.cvps = function(x, ...) {
  print(cbind(group=x$group, ni=x$ni, mi=x$mi, ri=x$ri), ...)
  cat("attr(,\"class\")\n")
  print(class(x))
}

length.cvps = function(x) length(x$ni)

weights.cvps= function(x, beta) 1

suppspace.cvps = function(x, beta) c(-Inf,Inf)

gridpoints.cvps = function(x, beta, grid=100) {
  r = range(x$mi)
  seq(r[1], r[2], len=grid)
}

initial.cvps = function(x, beta=NULL, mix=NULL, kmax=NULL) {
  if(is.null(beta)) beta = sqrt(sum(x$ri) / (sum(x$ni) - length(x)))
  if(length(beta) != 1 || beta <= 0) stop("initial beta incorrect")
  if(is.null(kmax)) kmax = 10
  if(is.null(mix) || is.null(mix$pt)) 
    mix = disc(quantile(x$mi, p=seq(0,1,len=kmax), type=4))
  list(beta=beta, mix=mix)
}

valid.cvps = function(x, beta, theta) beta > 0

# Compute the log density and its derivatives, for each combination of x[i]
# and pt[j].

logd.cvps = function(x, beta, pt, which=c(1,0,0)) {
  dl = vector("list", 3)
  names(dl) = c("ld","db","dt")
  if(which[1] == 1)
    dl$ld = - x$ni * 0.5 * log(2*pi*beta^2) -
      (x$ri + x$ni * outer(x$mi, pt, "-")^2) / beta^2 * 0.5
  if(which[2] == 1)
    dl$db = array(- x$ni / beta +
      (x$ri + x$ni * outer(x$mi, pt, "-")^2) / beta^3,
      dim=c(length(x$mi), length(pt), 1))
  if(which[3] == 1)
    dl$dt = x$ni / beta^2 * outer(x$mi, pt, "-")
  dl
}
