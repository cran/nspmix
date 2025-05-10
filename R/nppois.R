# ============================== #
# Nonparametric Poisson Mixtures #
# ============================== #

# Class nppois #

# n       number of observations
# mix     mixing distribution, an object of "disc", or "function" for
#         generating random support points (useful for a continuous
#         mixing distribution) 
# ...     arguments to be passed to function 'mix()'

rnppois = function(n, mix=disc(1), ...) {
  switch(class(mix)[1],
         disc = {
           k = length(mix$pt)
           ma = max(mix$pt)
           ## x = 0:round(20+ma+sqrt(ma)*15)
           x = 0:(qpois(1e-3/n, ma, lower.tail=FALSE) + 10)
           m = length(x)
           d = dpois(rep(x,k), rep(mix$pt,each=m)) * rep(mix$pr, each=m)
           dim(d) = c(m, k)
           d = rowSums(d)
           d = d / sum(d)
           w = drop(rmultinom(1,n,d))
           j = w != 0
           structure(list(v=x[j], w=w[j]), class="nppois")
         },
         "function" = {
           pt = mix(n, ...)
           x = rpois(n, pt)
           w = table(x)
           v = as.integer(names(w))
           names(w) = NULL
           j = w != 0
           structure(list(v=v[j], w=w[j]), class="nppois")
         },
         stop("`class(mix)` incompatible")
         )
}

##' @title Class \code{nppois}
##' 
##' @usage
##' nppois(v, w=1, grouping=TRUE)
##' rnppois(n, mix=disc(1), ...)
##' dnppois(x, mix=disc(1), log=FALSE)
##' pnppois(x, mix=disc(1), lower.tail=TRUE, log.p=FALSE)
##' 
##' @description Class \code{nppois} is used to store data that will
##'   be processed as those of a nonparametric Poisson mixture.
##' 
##' Function \code{nppois} creates an object of class \code{nppois},
##'  given values and weights/frequencies.
##' 
##' Function \code{rnppois} generates a random sample from a Poisson
##'  mixture and saves the data as an object of class \code{nppois}.
##' 
##' Function \code{dnppois} is the density function of a Poisson
##'  mixture.
##' 
##' Function \code{pnppois} is the distribution function of a Poisson
##'  mixture.
##' 
##' @aliases nppois rnppois dnppois pnppois
##'  
##' @param v a numeric vector that stores the values of a sample.
##' @param w a numeric vector that stores the corresponding
##'   weights/frequencies of the observations.
##' @param n the sample size.
##' @param x an object of class \code{nppois}.
##' @param mix an object of class \code{disc}.
##' @param grouping logical, to use frequencies (w) for identical
##'   values
##' @param log logical, to compute the log-values or not.
##' @param lower.tail, =FALSE, if lower.tail values are to be
##'   returned.
##' @param log.p, =FALSE, if log probability values are to be
##'   returned.
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
##' @examples
##' mix = disc(pt=c(0,4), pr=c(0.3,0.7))  # a discrete distribution
##' x = rnppois(200, mix)
##' dnppois(0:10, mix)
##' pnppois(0:10, mix)
##' dnppois(nppois(0:10), mix)
##' pnppois(nppois(0:10), mix)
##'  
##' @importFrom stats dpois ppois qpois rpois
##' 
##' @export nppois
##' @export rnppois
##' @export dnppois
##' @export pnppois

nppois = function(v, w=1, grouping=TRUE) {
  if("nppois" %in% class(v)) v$w = v$w * w
  else {
    if((is.data.frame(v) || is.matrix(v)) && ncol(v) == 2)
      r = list(v=v[,1], w=v[,2]*w)
    else r = list(v=v, w=w)
    v = structure(r, class="nppois")
  }
  if(! grouping) return(v)
  w = rep(v$w, len=length(v$v))
  r = aggregate(w, by=list(group=v$v), sum)
  structure(list(v=r$group,w=r$x), class="nppois")
}

##' @export 

length.nppois = function(x) length(x$v)

##' @export 

weight.nppois = function(x, beta) x$w

# lower and upper bounds on theta

##' @export 

suppspace.nppois = function(x, beta) c(0,Inf)

##' @export 

gridpoints.nppois = function(x, beta, grid=100) {
  bs = suppspace(x, beta)
  r = sqrt(range(x$v))
  pt = seq(max(r[1],sqrt(bs[1])), min(r[2],sqrt(bs[2])), len=grid)^2
  pmin(pmax(pt,bs[1]), bs[2])
}

# beta        Not used
# mix         Discrete mixing distribution

##' @export 

initial.nppois = function(x, beta=NULL, mix=NULL, kmax=NULL) {
  if(is.null(mix) || is.null(mix$pt)) {
    mi = max(floor(sqrt(min(x$v)))-1, 0)
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

##' @export 

valid.nppois = function(x, beta, theta) TRUE

# No beta

##' @export 

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

##' @title Plotting a nonparametric Poisson mixture
##' 
##' @description Function \code{plot.nppois} plots a Poisson mixture,
##'    along with data.
##' 
##' @param x an object of class \code{nppois}.
##' @param mix an object of class \code{disc}.
##' @param beta the structural parameter (not used for a Poisson
##'   mixture).
##' @param col the color of the density curve to be plotted.
##' @param add if \code{FALSE}, creates a new plot; if \code{TRUE},
##'   adds the plot to the existing one.
##' @param components if \code{proportions} (default), also show the
##'   support points and mixing proportions (in proportional vertical
##'   lines); if \code{curves}, also show the component density
##'   curves; if \code{null}, components are not shown.
##' @param main,lwd,lty,xlab,ylab,xlim arguments for graphical
##'   parameters (see \code{par}).
##' @param ... arguments passed on to function \code{barplot}.
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
##' mix = disc(pt=c(0,4), pr=c(0.3,0.7))  # a discrete distribution
##' x = rnppois(200, mix)
##' plot(x, mix)
##'  
##' @export 

plot.nppois = function(x, mix, beta, col=2, add=FALSE,
                       components=c("proportions","curves","null"),
                       main="nppois", lwd=1,
                       lty=1, xlab="Data", ylab="Density", xlim=NULL, 
                       ...) {
  components = match.arg(components)
  if(missing(x) && missing(mix))
    stop('Both arguments "x" and "mix" are missing')

  if(is.null(xlim)) {
    if(missing(x))
      xlim = floor(c(max(ptr[1] - 4 * sqrt(ptr[1]), 0),
                     ptr[2] + 4 * sqrt(ptr[2])))
    else xlim = c(0,max(x$v)+1)
  }
  if(! missing(x) && ! add) {
    maxb = max(xlim[2], max(x$v)+1)
    if(maxb < 1000) { by =1; breaks = 0:maxb - 0.5 }
    else { by = ceiling(maxb / 1000); breaks = seq(0, maxb+by, by=by) - 0.5 }
    f2 = whist(x$v, x$w, breaks=breaks, freq=FALSE, plot=FALSE)$density
    names(f2) = round((breaks[-1] + breaks[-length(breaks)]) * 0.5)
    barplot(f2, by, space=0, xlab=xlab, ylab=ylab, main=main,
            col="lightgrey", xlim=xlim, ...)
  }
  if(! missing(mix)) {
    ptr = range(mix$pt)
    if(xlim[2] - xlim[1] < 1000) y = xlim[1]:xlim[2]
    else y = unique(round(seq(xlim[1], xlim[2], len=1000)))
    dj = outer(y, mix$pt, dpois) * rep(mix$pr, each=length(y))
    d = rowSums(dj)
    if(! missing(x) || add) lines(y + 0.5, d, col=col, lty=lty, lwd=lwd)
    else plot(y + 0.5, d, type="l", col=col, lty=lty, lwd=lwd,
              xlim=xlim, ...)
    switch(components,
           proportions = {
             segments(mix$pt + 0.5, rep(0,length(mix$pt)),
                      y1=mix$pr*max(d)*0.8,
                      col=col, lwd=3)
             points(mix$pt + 0.5, rep(0, length(mix$pt)), col=col)
           },
           curves = {
             matplot(y+0.5, dj, type="l", lty=lty+1, col=col, add=TRUE)
             points(mix$pt + 0.5, rep(0, length(mix$pt)), col=col)
           }
           )
    ## if(components) {
    ## matplot(y+0.5, dj, type="l", lty=lty+1, col=col, add=TRUE)
    ## points(mix$pt + 0.5, rep(0, length(mix$pt)), col=col)
  }
}


## Additional functions

## x     vector

dnppois = function(x, mix=disc(1), log=FALSE) {
  if("nppois" %in% class(x)) x = x$v
  logd = outer(x, mix$pt, dpois, log=TRUE) +
    rep(log(mix$pr), rep(length(x), length(mix$pr)))
  if(log) {
    ma = matMaxs(logd)
    ma + log(rowSums(exp(logd - ma)))
  }
  else rowSums(exp(logd))
}

## x     vector

pnppois = function(x, mix=disc(1), lower.tail=TRUE, log.p=FALSE) {
  if("nppois" %in% class(x)) x = x$v
  logp = outer(x, mix$pt, ppois, lower.tail=lower.tail, log=TRUE) +
    rep(log(mix$pr), rep(length(x), length(mix$pr)))
  if(log.p) {
    ma = matMaxs(logp)
    ma + log(rowSums(exp(logp - ma)))
  }
  else rowSums(exp(logp))
}

##' @title Sorting of an Object of Class \code{nppois}
##' 
##' @description Function \code{sort.nppois} sorts an object of class
##'   \code{nppois} in the order of the obsersed values.
##' 
##' @param x an object of class \code{nppois}.
##' @param decreasing logical, in the decreasing (default) or
##'    increasing order.
##' @param ... arguments passed to function \code{order}.
##'  
##' @examples
##' 
##' mix = disc(pt=c(0,4), pr=c(0.3,0.7))  # a discrete distribution
##' x = rnppois(20, mix)
##' sort(x)
##' 
##'  
##' @export 

sort.nppois = function(x, decreasing=FALSE, ...) {
  i = order(x$v, decreasing=decreasing, ...)
  structure(list(v=x$v[i], w=if(length(x$w)==1) x$w else x$w[i]),
            class="nppois")
}

