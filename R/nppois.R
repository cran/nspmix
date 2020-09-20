# ============================== #
# Nonparametric Poisson Mixtures #
# ============================== #

# Class nppois #

# mix     Object of "disc"

rnppois = function(n, mix=disc(1)) {
  k = length(mix$pt)
  ma = max(mix$pt)
  x = 0:round(20+ma+sqrt(ma)*15)
  m = length(x)
  d = dpois(rep(x,k), rep(mix$pt,each=m)) * rep(mix$pr, each=m)
  dim(d) = c(m, k)
  d = rowSums(d)
  d = d / sum(d)
  w = drop(rmultinom(1,n,d))
  j = w != 0
  structure(list(v=x[j], w=w[j]), class="nppois")
}



##'Class `nppois'
##'
##'
##'Class \code{nppois} is used to store data that will be processed as those of
##'a nonparametric Poisson mixture.
##'
##'Function \code{nppois} creates an object of class \code{nppois}, given
##'values and weights/frequencies.
##'
##'Function \code{rnppois} generates a random sample from a Poisson mixture and
##'saves the data as an object of class \code{nppois}.
##'
##'Function \code{plot.nppois} plots the Poisson mixture.
##'
##'
##'When \code{components=TRUE}, the support points are shown on the horizontal
##'line of density 0. The component density curves, weighted appropriately, are
##'also shown.
##'
##'@aliases nppois rnppois plot.nppois
##'@param v a numeric vector that stores the values of a sample.
##'@param w a numeric vector that stores the corresponding weights/frequencies
##'of the observations.
##'@param n the sample size.
##'@param x an object of class \code{nppois}.
##'@param mix an object of class \code{disc}.
##'@param beta the structural parameter, which is not really needed for the
##'Poisson mixture.
##'@param col the color of the density curve to be plotted.
##'@param add if \code{FALSE}, creates a new plot; if \code{TRUE}, adds the
##'plot to the existing one.
##'@param components if \code{TRUE}, also show the support points and mixing
##'proportions.
##'@param main,lwd,lty,xlab,ylab arguments for graphical parameters (see
##'\code{par}).
##'@param ... arguments passed on to function \code{plot}.
##'@author Yong Wang <yongwang@@auckland.ac.nz>
##'@seealso \code{\link{nnls}}, \code{\link{cnm}}, \code{\link{cnmms}},
##'\code{\link{plot.nspmix}}.
##'@references
##'
##'Wang, Y. (2007). On fast computation of the non-parametric maximum
##'likelihood estimate of a mixing distribution. \emph{Journal of the Royal
##'Statistical Society, Ser. B}, \bold{69}, 185-198.
##'@keywords class function
##'@examples
##'
##'mix = disc(pt=c(1,4), pr=c(0.3,0.7))
##'x = rnppois(200, mix)
##'plot(x, mix)
##'
##'@usage
##'nppois(v, w=1)
##'rnppois(n, mix=disc(1))
##'\method{plot}{nppois}(x, mix, beta, col="red", add=FALSE,
##'     components=TRUE, main="nppois", lwd=1, lty=1, xlab="Data",
##'     ylab="Density", ...)
##' 
##'@export nppois
##'@export rnppois
##'@export plot.nppois

nppois = function(v, w=1) {
  if(class(v) == "nppois") {
    v$w = v$w * w
    v
  }
  else {
    if((is.data.frame(v) || is.matrix(v)) && ncol(v) == 2)
      r = list(v=v[,1], w=v[,2]*w)
    else r = list(v=v, w=w)
    structure(r, class="nppois")
  }
}

length.nppois = function(x) length(x$v)

weight.nppois = function(x, beta) x$w

# lower and upper bounds on theta

suppspace.nppois = function(x, beta) c(0,Inf)

gridpoints.nppois = function(x, beta, grid=100) {
  bs = suppspace(x, beta)
  r = sqrt(range(x$v))
  pt = seq(max(r[1],sqrt(bs[1])), min(r[2],sqrt(bs[2])), len=grid)^2
  pmin(pmax(pt,bs[1]), bs[2])
}

# beta        Not used
# mix         Discrete mixing distribution

initial.nppois = function(x, beta=NULL, mix=NULL, kmax=NULL) {
  if(is.null(mix) || is.null(mix$pt)) {
    mi = min(ceiling(sqrt(min(x$v))), 0)
    ma = max(ceiling(sqrt(max(x$v))), 1)
    if(is.null(kmax)) breaks = (mi:ma)^2
    else {
      if(kmax == 1)
        return(list(beta=beta,
                    mix=disc(sum(x$v*x$w) / sum(rep(x$w,len=length(x$v))))))
      breaks = seq(mi, ma, len=kmax+1)^2
    }
    r = nspmix::whist(x$v, x$w, breaks=breaks, freq=FALSE, plot=FALSE)
    i = r$density != 0 
    mix = disc(r$mids[i], r$density[i])
  }
  list(beta=beta, mix=mix)
}

valid.nppois = function(x, beta, mix) TRUE

# No beta

logd.nppois = function(x, beta, pt, which=c(1,0,0)) {
  dl = vector("list", 3)
  names(dl) = c("ld","db","dt")
  n = length(x$v)
  k = length(pt)
  rpt = rep(pmax(pt, 1e-100), rep(n,k))
  if(which[1] == 1) {
    dl$ld = x$v * log(rpt) - rpt - lfactorial(x$v)
    dim(dl$ld) = c(n, k)
  }
  if(which[3] == 1) {
    dl$dt = x$v / rpt - 1
    dim(dl$dt) = c(n, k)
  }
  dl
}

plot.nppois = function(x, mix, beta, col="red", add=FALSE,
                      components=TRUE, main="nppois", lwd=1,
                      lty=1, xlab="Data", ylab="Density", ...) {
  ptr = range(mix$pt)
  if(is.null(x))
    xlim = floor(c(max(ptr[1] - 4 * sqrt(ptr[1]), 0),
                   ptr[2] + 4 * sqrt(ptr[2])))
  else xlim = range(x$v)
  y = xlim[1]:xlim[2]
  if(!is.null(x) && !add) {
##    nspmix::whist(x$v, x$w, breaks=50, freq=FALSE)
    f2 = rep(0, xlim[2]+1)
    names(f2) = 0:xlim[2]
    f2[paste(x$v)] = x$w
    f2 = f2 / sum(f2)
    barplot(f2, 1, space=0, xlab=xlab, ylab=ylab, main=main, col="lightgrey",
            xlim=xlim, ...)
  }
  d = 0
  for(j in 1:length(mix$pt)) {
    dj = dpois(y, mix$pt[j]) * mix$pr[j]
    if(components) lines(y + 0.5, dj, lty=lty+1, col=col)
    d = d + dj
  }
  if(components) {
    points(mix$pt + 0.5, rep(0, length(mix$pt)), col=col)
  }
  lines(y + 0.5, d, col=col, lty=lty, lwd=lwd)
}


