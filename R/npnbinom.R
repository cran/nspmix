# ========================================= #
# Nonparametric Negagtive Binomial Mixtures #
# ========================================= #

# Class npnbinom #

# mix     Object of "disc"

rnpnbinom = function(n, size, mix=disc(0.5)) {
  k = length(mix$pt)
  mi = min(mix$pt)
  x = 0:(qnbinom(1e-3/n, size, mi, lower.tail=FALSE) + 10)
  m = length(x)
  d = dnbinom(rep(x,k), size, rep(mix$pt,each=m)) * rep(mix$pr, each=m)
  dim(d) = c(m, k)
  d = rowSums(d)
  d = d / sum(d)
  w = drop(rmultinom(1,n,d))
  j = w != 0
  structure(list(v=x[j], w=w[j], size=size), class="npnbinom")
}



##' @title Class \code{npnbinom}
##'
##' @usage
##' npnbinom(v, w=1, size, grouping=TRUE)
##' rnpnbinom(n, size, mix=disc(0.5))
##' dnpnbinom(x, mix=disc(0.5), size=NULL, log=FALSE)
##' pnpnbinom(x, mix=disc(0.5), size=NULL, lower.tail=TRUE, log.p=FALSE)
##'
##' @description Class \code{npnbinom} is used to store data that will be
##'   processed as those of a nonparametric negative binomial mixture.
##'
##' Function \code{npnbinom} creates an object of class
##' \code{npnbinom}, given values and weights/frequencies.
##'
##' Function \code{rnpnbinom} generates a random sample from a
##' negative binomial mixture and saves the data as an object of class
##' \code{npnbinom}.
##'
##' @aliases npnbinom rnpnbinom dnpnbinom pnpnbinom
##' @param v a numeric vector that stores the values of a sample.
##' @param w a numeric vector that stores the corresponding
##'   weights/frequencies of the observations.
##' @param size number of successful trials (ignored if \code{x} is
##'   an object of class \code{npnbinom}).
##' @param n the sample size.
##' @param x an object of class \code{npnbinom}, or a numeric vector
##'   (then value of \code{size} must be provided).
##' @param mix an object of class \code{disc}.
##' @param grouping logical, to use frequencies (w) for identical
##'   values
##' @param log, =FALSE, if log-values are to be returned.
##' @param lower.tail, =FALSE, if lower.tail values are to be
##'   returned.
##' @param log.p, =FALSE, if log probability values are to be
##'   returned.
##'
##' @author Yong Wang <yongwang@@auckland.ac.nz>
##' 
##' @seealso \code{\link[lsei]{nnls}}, \code{\link{cnm}},
##'   \code{\link{cnmms}}, \code{\link{plot.nspmix}}.
##' 
##' @references
##'
##'   Wang, Y. (2007). On fast computation of the non-parametric
##'   maximum likelihood estimate of a mixing
##'   distribution. \emph{Journal of the Royal Statistical Society,
##'   Ser. B}, \bold{69}, 185-198.
##' 
##' @examples
##' 
##' mix = disc(pt=c(0.2,0.5), pr=c(0.3,0.7))
##' (x = rnpnbinom(200, size=10, mix))
##' dnpnbinom(x, mix, size=10)
##' pnpnbinom(x, mix, size=10)
##'
##' @importFrom stats dnbinom pnbinom qnbinom rbinom
##' 
##' @export npnbinom
##' @export rnpnbinom
##' @export dnpnbinom
##' @export pnpnbinom

npnbinom = function(v, w=1, size, grouping=TRUE) {
  if("npnbinom" %in% class(v)) v$w = v$w * w
  else {
    if((is.data.frame(v) || is.matrix(v)) && ncol(v) == 2)
      r = list(v=v[,1], w=v[,2]*w)
    else r = list(v=v, w=w, size=size)
    v = structure(r, class="npnbinom")
  }
  if(! grouping) return(v)
  w = rep(v$w, len=length(v$v))
  r = aggregate(w, by=list(group=v$v), sum)
  structure(list(v=r$group, w=r$x, size=v$size), class="npnbinom")
}

##' @export 

length.npnbinom = function(x) length(x$v)

##' @export 

weight.npnbinom = function(x, beta) x$w

# lower and upper bounds on theta

##' @export 

suppspace.npnbinom = function(x, beta) c(0, 1)

##' @export 

gridpoints.npnbinom = function(x, beta, grid=100) {
  log.ph = log(x$size / (x$size + x$v))
  exp(seq(min(log.ph), max(log.ph), len=grid))      # may include 1
}

# beta        Not used
# mix         Discrete mixing distribution

##' @export 

initial.npnbinom = function(x, beta=NULL, mix=NULL, kmax=NULL) {
  if(is.null(mix) || is.null(mix$pt)) 
    mix = disc(gridpoints.npnbinom(x, beta, min(20,kmax)))
  list(beta=beta, mix=mix)
}

##' @export 

valid.npnbinom = function(x, beta, theta) TRUE

# No beta

##' @export 

logd.npnbinom = function(x, beta, pt, which=c(1,0,0)) {
  dl = vector("list", 3)
  names(dl) = c("ld","db","dt")
  n = length(x$v)
  k = length(pt)
  rpt = rep(pmin(pmax(pt, 1e-100), 1-1e-10), rep(n,k))
  if(which[1] == 1) {
    dl$ld = lchoose(x$size + x$v - 1, x$size - 1) +
      x$size * log(rpt) + x$v * log(1 - rpt)
    dim(dl$ld) = c(n, k)
  }
  if(which[3] == 1) {
    dl$dt = x$size / rpt - x$v / (1 - rpt)
    dim(dl$dt) = c(n, k)
  }
  dl
}

##' @title Plotting a nonparametric negative binomial mixture
##' 
##' @description Function \code{plot.npnbinom} plots a negative
##'   binomial mixture, along with data.
##' 
##' @param x an object of class \code{npnbinom}.
##' @param mix an object of class \code{disc}.
##' @param beta the structural parameter (not used for a negative
##'   binomial mixture).
##' @param col the color of the density curve to be plotted.
##' @param add if \code{FALSE}, creates a new plot; if \code{TRUE},
##'   adds the plot to the existing one.
##' @param components if \code{TRUE}, also show the support points and
##'   mixing proportions (with vertical lines in proportion).
##' @param main,lwd,lty,xlab,ylab arguments for graphical parameters
##'   (see \code{par}).
##' @param ... arguments passed on to function \code{plot}.
##'  
##' @author Yong Wang <yongwang@@auckland.ac.nz>
##'  
##' @seealso \code{\link[lsei]{nnls}}, \code{\link{cnm}},
##'   \code{\link{cnmms}}, \code{\link{plot.nspmix}}.
##'  
##' @references
##' 
##' Wang, Y. (2007). On fast computation of the non-parametric maximum
##' likelihood estimate of a mixing distribution. \emph{Journal of the
##' Royal Statistical Society, Ser. B}, \bold{69}, 185-198.
##'  
##' @keywords class function
##'  
##' @examples
##' 
##' mix = disc(pt=c(0.2,0.5), pr=c(0.3,0.7))  # a discrete distribution
##' x = rnpnbinom(200, 10, mix)
##' plot(x, mix)
##'  
##' @export 

plot.npnbinom = function(x, mix, beta, col="red", add=FALSE,
                         components=TRUE, main="npnbinom", lwd=1,
                         lty=1, xlab="Data", ylab="Density", ...) {
  pt.min = min(mix$pt)
  if(is.null(x)) xlim = c(0, round(log(1e-3/pt.min) / log(pt.min)))
  else xlim = range(x$v)
  y = xlim[1]:xlim[2]
  if(!is.null(x) && !add) {
    breaks = 0:max(xlim[2], max(x$v)+1) - 0.5
    f2 = whist(x$v, x$w, breaks=breaks, freq=FALSE, plot=FALSE)$density
    names(f2) = breaks[-1] - 0.5
    barplot(f2, 1, space=0, xlab=xlab, ylab=ylab, main=main,
            col="lightgrey", xlim=xlim, ...)
    # f2 = rep(0, xlim[2]+1)
    # names(f2) = 0:xlim[2]
    # f2[paste(x$v)] = x$w
    # f2 = f2 / sum(f2)
    # barplot(f2, 1, space=0, xlab=xlab, ylab=ylab, main=main, col="lightgrey",
    #         xlim=xlim, ...)
  }
  # d = 0
  # for(j in 1:length(mix$pt)) {
  #   dj = dnbinom(y, x$size, mix$pt[j]) * mix$pr[j]
  #   if(components) lines(y + 0.5, dj, lty=lty+1, col=col)
  #   d = d + dj
  # }
  dj = outer(y, mix$pt, dnbinom, size=x$size) * rep(mix$pr, each=length(y))
  if(components) matplot(y + 0.5, dj, type="l", lty=lty+1, col=col, add=TRUE)
  d = rowSums(dj)
  
  ## if(components) points(mix$pt + 0.5, rep(0, length(mix$pt)), col=col)
  lines(y + 0.5, d, col=col, lty=lty, lwd=lwd)
}


## Additional functions

dnpnbinom = function(x, mix=disc(0.5), size=NULL, log=FALSE) {
  if("npnbinom" %in% class(x)) {size = x$size; x = x$v}
  if(is.null(size)) stop("Value of 'size' not provided.")
  logd = outer(x, mix$pt, dnbinom, size=size, log=TRUE) +
    rep(log(mix$pr), rep(length(x), length(mix$pr)))
  if(log) {
    ma = matMaxs(logd)
    ma + log(rowSums(exp(logd - ma)))
  }
  else rowSums(exp(logd))
}

pnpnbinom = function(x, mix=disc(0.5), size=NULL, lower.tail=TRUE,
                     log.p=FALSE) {
  if("npnbinom" %in% class(x)) {size = x$size; x = x$v}
  if(is.null(size)) stop("Value of 'size' not provided.")
  logp = outer(x, mix$pt, pnbinom, size=size,
               lower.tail=lower.tail, log=TRUE) +
    rep(log(mix$pr), rep(length(x), length(mix$pr)))
  if(log.p) {
    ma = matMaxs(logp)
    ma + log(rowSums(exp(logp - ma)))
  }
  else rowSums(exp(logp))
}

