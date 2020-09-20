# ============================= #
# Nonparametric normal mixtures #
# ============================= #

# ============ #
# Class npnorm #
# ============ #

# mix     Object of "disc"

rnpnorm = function(n, mix=disc(0), sd=1) {
  if (n == 0) return(numeric(0))
  k = length(mix$pt)
  suppressWarnings(i <- sample.int(k, n, prob = mix$pr, replace = TRUE))
  x = list(v=mix$pt[i] + rnorm(n, sd=sd), w=1)
  structure(x, class= "npnorm")
}



##'Class `npnorm'
##'
##'
##' Class \code{npnorm} can be used to store data that will be
##' processed as those of a nonparametric normal mixture. There are
##' several functions associated with the class.
##'
##' Function \code{npnorm} creates an object of class \code{npnorm},
##' given values and weights/frequencies.
##'
##' Function \code{rnpnorm} generates a random sample from a normal
##' mixture and saves the data as an object of class \code{npnorm}.
##'
##' Function \code{plot.npnorm} plots the normal mixture.
##'
##'
##' When \code{components="proportions"}, the component means are
##' shown on the horizontal line of density 0. The vertical lines
##' going upwardly at the support points are proportional to the
##' mixing proportions at these points.
##'
##'@aliases npnorm rnpnorm plot.npnorm
##'@param v a numeric vector that stores the values of a sample.
##'@param w a numeric vector that stores the corresponding weights/frequencies
##'of the observations.
##'@param n the sample size.
##'@param mix an object of class \code{disc}, for a discrete distribution.
##'@param beta the structural parameter.
##'@param sd a scalar for the component standard deviation that is common to
##'all components.
##'@param x an object of class \code{npnorm}.
##'@param breaks the rough number bins used for plotting the histogram.
##'@param col the color of the density curve to be plotted.
##'@param len the number of points roughly used to plot the density curve over
##'the interval of length 8 times the component standard deviation around each
##'component mean.
##'@param add if \code{FALSE}, creates a new plot; if \code{TRUE}, adds the
##'plot to the existing one.
##'@param border.col color for the border of histogram boxes.
##'@param border.lwd line width for the border of histogram boxes.
##'@param fill color to fill in the histogram boxes.
##'@param components if \code{proportions} (default), also show the support
##'points and mixing proportions; if \code{curves}, also show the component
##'density curves; if \code{null}, components are not shown.
##'@param lty.components,lwd.components line type and width for the component
##'curves.
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
##'mix = disc(pt=c(0,4), pr=c(0.3,0.7))  # a discrete distribution
##'x = rnpnorm(200, mix, sd=1)
##'plot(x, mix, beta=1)
##'
##'@usage
##'npnorm(v, w=1)
##'rnpnorm(n, mix=disc(0), sd=1)
##'\method{plot}{npnorm}(x, mix, beta, breaks=NULL, col="red", len=100,
##'     add=FALSE, border.col=NULL, border.lwd=1,
##'     fill="lightgrey", main, lwd=2, lty=1, xlab="Data",
##'     ylab="Density", components=c("proportions","curves","null"),
##'     lty.components=2, lwd.components=2, ...) 
##' 
##'@export npnorm
##'@export rnpnorm
##'@export plot.npnorm

npnorm = function(v, w=1) {
  if("npnorm" %in% class(v)) { v$w = v$w * w; v }
  else structure(list(v=v, w=w), class="npnorm")
}

length.npnorm = function(x) length(x$v)

weight.npnorm = function(x, beta) x$w

gridpoints.npnorm = function(x, beta, grid=100) {
  rx = range(x$v)
  breaks = max(ceiling(diff(rx) / (5*beta)), 5)   # number of breaks
  r = whist(x$v, x$w, breaks=breaks, freq=FALSE, plot=FALSE)
  i = r$density != 0
  i = i | c(i[-1],FALSE) | c(FALSE,i[-length(i)])  # include neighbours
  m = sum(i)
  k = max(ceiling(grid / m), 10)           # at least 10 in each interval
  d = r$breaks[2] - r$breaks[1]
  s = r$breaks[-length(r$breaks)][i]
  c(rx[1], rep(s, rep(k,m)) + d * (1:k-0.5)/k, rx[2])
}

# beta        Standard deviation, with default value = 1
# mix         Discrete mixing distribution

initial.npnorm = function(x, beta=NULL, mix=NULL, kmax=NULL) {
  if(is.null(beta)) beta = 1
  if(is.null(mix) || is.null(mix$pt)) {
    if(is.null(kmax)) breaks = max(ceiling(diff(range(x$v)) / (5*beta)), 5)
    else {
      if(kmax == 1)
        return(list(beta=beta,
                    mix=disc(sum(x$v*x$w) / sum(rep(x$w,len=length(x$v))))))
      breaks = seq(min(x$v), max(x$v), len=min(20, kmax))
    }
    r = whist(x$v, x$w, breaks=breaks, freq=FALSE, plot=FALSE)
    i = r$density != 0 
    mix = disc(r$mids[i], r$density[i])
  }
  list(beta=beta, mix=mix)
}

valid.npnorm = function(x, beta, theta) beta > 0

suppspace.npnorm = function(x, beta) c(-Inf,Inf)

logd.npnorm = function(x, beta, pt, which=c(1,0,0)) {
  dl = vector("list", 3)
  names(dl) = c("ld","db","dt")
  x.pt = x$v - rep(pt, rep(length(x$v), length(pt)))
  dim(x.pt) = c(length(x$v), length(pt))
  beta2 = beta * beta
  if(any(which[1:2] == 1)) xpt2beta2 = x.pt^2 / beta2
  if(which[1] == 1)
    dl$ld = - .5 * (log(2*base::pi*beta2) + xpt2beta2)
  if(which[2] == 1) {
    dl$db = (xpt2beta2 - 1) / beta
    dim(dl$db) = c(length(x$v), length(pt), 1)
  }
  if(which[3] == 1)
    dl$dt = x.pt / beta2
  dl
}

plot.npnorm = function(x, mix, beta, breaks=NULL, col="red", len=100,
                       add=FALSE, border.col=NULL, border.lwd=1,
                       fill="lightgrey", main, lwd=2, lty=1, xlab="Data",
                       ylab="Density", 
                       components=c("proportions","curves","null"),
                       lty.components=2, lwd.components=2, ...) {
  components = match.arg(components)
  m = length(mix$pt)
  z = sort(unique(round(mix$pt +
                        rep(beta*c(-10*exp(10*10:0), seq(-4,4,len=len),
                                   10*exp(10*10:0)),
                            each=m),
                        ceiling(-log10(beta/20)))))
  nz = length(z)
  dj = outer(z, mix$pt, dnorm, sd=beta) * rep(mix$pr, rep(nz,length(mix$pr)))
  d = rowSums(dj)
  if(missing(main))
    main = substitute("npnorm (" * sigma ~ "=" ~ xxx * ")",
                      list(xxx=signif(beta,3)))
  if(add || missing(x) || length(x) == 0) {
    if(add) lines(z, d, col=col, lwd=lwd, lty=lty)
    else {
      plot(0, 0, type="n", xlim=range(mix$pt) + beta * c(-3,3), ylim=range(d),
           xlab=xlab, ylab=ylab, frame.plot=FALSE,
           main=main, ...)
      lines(range(z), rep(0,2), col="darkgrey")
      lines(z, d, col=col, lwd=lwd, lty=lty) 
    }
  }
  else {
    if(is.null(breaks)) breaks = 10 + round(sqrt(length(x$v)))
    whist(x$v, x$w, breaks=breaks, freq=FALSE, 
          xlab=xlab, ylab=ylab, main=main, col=fill, border=border.col, 
          lwd=border.lwd, ylim=range(d), ...)
    lines(z, d, col=col, lwd=lwd, lty=lty)
  }
  if(components != "null") {
    points(mix$pt, rep(0,length(mix$pt)), col=col)
    switch(components,
           proportions = {
             segments(mix$pt, rep(0,m), y1=mix$pr*max(d), col=col, lwd=3)
           },
           curves = {
             if(ncol(dj) > 1) {
               j2 = 1:floor(ncol(dj)/2) * 2
               for(j in 1:ncol(dj))     # slow if too many components
                 lines(z, dj[,j], col=col, lty=lty.components,
                       lwd=lwd.components)
             }
           }
           )
  }
}

## Additional functions

## x     vector

dnpnorm = function(x, mix=disc(0), sd=1, log=FALSE) {
  if("npnorm" %in% class(x)) x = x$v
  logd = outer(x, mix$pt, dnorm, sd=sd, log=TRUE) +
    rep(log(mix$pr), rep(length(x), length(mix$pr)))
  if(log) {
    ma = matMaxs(logd)
    ma + log(rowSums(exp(logd - ma)))
  }
  else rowSums(exp(logd))
}

## x     vector

pnpnorm = function(x, mix=disc(0), sd=1, lower.tail=TRUE, log.p=FALSE) {
  if("npnorm" %in% class(x)) x = x$v
  logp = outer(x, mix$pt, pnorm, sd=sd, lower.tail=lower.tail, log=TRUE) +
    rep(log(mix$pr), rep(length(x), length(mix$pr)))
  if(log.p) {
    ma = matMaxs(logp)
    ma + log(rowSums(exp(logp - ma)))
  }
  else rowSums(exp(logp))
}

