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



##'Class `cvps'
##'
##'
##'These functions can be used to study a common variance problem (CVP), where
##'univariate observations fall in known groups. Observations in each group are
##'assumed to have the same mean, but different groups may have different
##'means. All observations are assumed to have a common variance, despite their
##'different means, hence giving the name of the problem. It is a
##'random-effects problem.
##'
##'Class \code{cvps} is used to store the CVP data in a summarized form.
##'
##'Function \code{cvps} creates an object of class \code{cvps}, given a matrix
##'that stores the values (column 2) and their grouping information (column 1).
##'
##'Function \code{rcvp} generates a random sample in the raw form for a common
##'variance problem, where the means follow a discrete distribution.
##'
##'Function \code{rcvps} generates a random sample in the summarized form for a
##'common variance problem, where the means follow a discrete distribution.
##'
##'Function \code{print.cvps} prints the CVP data given in the summarized form.
##'
##'
##'The raw form of the CVP data is a two-column matrix, where each row
##'represents an observation. The two columns along each row give,
##'respectively, the group membership (\code{group}) and the value (\code{x})
##'of an observation.
##'
##'The summarized form of the CVP data is a four-column matrix, where each row
##'represents the summarized data for all observations in a group. The four
##'columns along each row give, respectively, the group number (\code{group}),
##'the number of observations in the group (\code{ni}), the sample mean of the
##'observations in the group (\code{mi}), and the residual sum of squares of
##'the observations in the group (\code{ri}).
##'
##'@aliases cvp rcvp cvps rcvps print.cvps
##'@param x CVP data in the raw form as an argument in \code{cvps}, or an
##'object of class \code{cvps} in \code{print.cvps}.
##'@param k the number of groups.
##'@param ni a numeric vector that gives the sample size in each group.
##'@param mu a numeric vector for all the theoretical means.
##'@param pr a numeric vector for all the probabilities associated with the
##'theoretical means.
##'@param sd a scalar for the standard deviation that is common to all
##'observations.
##'@param ... arguments passed on to function \code{print}.
##'@author Yong Wang <yongwang@@auckland.ac.nz>
##'@seealso \code{\link{nnls}}, \code{\link{cnmms}}.
##'@references
##'
##'Neyman, J. and Scott, E. L. (1948). Consistent estimates based on partially
##'consistent observations. \emph{Econometrica}, \bold{16}, 1-32.
##'
##'Kiefer, J. and Wolfowitz, J. (1956). Consistency of the maximum likelihood
##'estimator in the presence of infinitely many incidental parameters.
##'\emph{Ann. Math. Stat.}, \bold{27}, 886-906.
##'
##'Wang, Y. (2010). Maximum likelihood computation for fitting semiparametric
##'mixture models. \emph{Statistics and Computing}, \bold{20}, 75-86.
##'@keywords class function
##'@examples
##'
##'
##'x = rcvps(k=50, ni=5:10, mu=c(0,4), pr=c(0.7,0.3), sd=3)
##'cnmms(x)              # CNM-MS algorithm
##'cnmpl(x)              # CNM-PL algorithm
##'cnmap(x)              # CNM-AP algorithm
##'
##'
##'@usage
##'cvps(x)
##'rcvp(k, ni=2, mu=0, pr=1, sd=1)
##'rcvps(k, ni=2, mu=0, pr=1, sd=1)
##'\method{print}{cvps}(x, ...)
##' 
##'@export cvps
##'@export rcvp
##'@export cvps
##'@export rcvps
##'@export print.cvps

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

weight.cvps= function(x, beta) 1

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
