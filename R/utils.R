## weighted histogram



##' @title Weighted Histograms
##' 
##' Plots or computes the histogram with observations with
##' multiplicities/weights.
##' 
##' Just like \code{hist}, \code{whist} can either plot the histogram
##' or compute the values that define the histogram, by setting
##' \code{plot} to \code{TRUE} or \code{FALSE}.
##' 
##' The histogram can either be the one for frequencies or density, by
##' setting \code{freq} to \code{TRUE} or \code{FALSE}.
##' 
##' @param x a vector of values for which the histogram is desired.
##' @param w a vector of multiplicities/weights for the values in
##'   \code{x}.
##' @param
##'   breaks,plot,freq,xlim,ylim,xlab,ylab,main,add,col,border,lwd
##'   These arguments have similar functionalities to their namesakes
##'   in function \code{hist}.
##' @param ... arguments passed on to function \code{plot}.
##' 
##' @return
##' 
##' \item{breaks}{the break points.}
##' 
##' \item{counts}{weighted counts over the intervals determined by
##' \code{breaks}}
##' 
##' \item{density}{density values over the intervals determined by
##' \code{breaks}}
##' 
##' \item{mids}{midpoints of the intervals determined by
##' \code{breaks}}
##' 
##' @author Yong Wang <yongwang@@auckland.ac.nz>
##' 
##' @seealso \code{\link{hist}}.
##' 
##' @export

whist = function(x, w=1, breaks="Sturges", plot=TRUE, freq=NULL, 
                 xlim=NULL, ylim=NULL, xlab="Data", ylab=NULL, main=NULL,
                 add=FALSE, col="lightgray", border=NULL, lwd=1, ...) {
  w = rep(w, len=length(x))
  o = order(x)
  x = x[o]
  w = w[o]
  r = hist(x, breaks=breaks, plot=FALSE)
  breaks = r$breaks
  nb = length(breaks)
  i = cut(x, breaks, include.lowest=TRUE)
  ii = as.integer(i)
  f = double(nb - 1)
  k = which(c(TRUE, diff(ii) > 0, TRUE))
  cw = c(0, cumsum(w))
  f[ii[k-1]] = diff(cw[k])                               # frequency
  
  # f2 = tapply(w, i, sum)                               # frequency
  # f2[is.na(f2)] = 0
  # dimnames(f2)[[1]] = NULL
  # print(cbind(f, f2))
  
  d = f / sum(f) / diff(breaks)                       # density
  if(! is.null(freq) && ! freq) {
    y = d
    if(is.null(ylab)) ylab = "Density"
  }
  else {
    y = f
    if(is.null(ylab)) ylab = "Frequency"
  }
  ny = length(y)
  if (is.null(xlim)) xlim = range(breaks)
  if(is.null(ylim)) ylim = range(0, y, finite=TRUE)
  else {
    ymax = max(y)
    ylim = c(0, pmin(pmax(ymax, max(ylim)), 2*ymax))
  }
  if(plot) {
    if(!add) plot(r$mids, y, xlim=xlim, ylim=ylim, type="n", frame.plot=FALSE,
                  xlab=xlab, ylab=ylab, main=main, ...)
    rect(breaks[-(ny+1)], 0, breaks[-1], y, col=col, border=border, lwd=lwd)
    lines(range(breaks), c(0,0), col=border)
  }
  else list(breaks=breaks, counts=f, density=d,
            mids=breaks[-(ny+1)] + diff(breaks) * .5)
}
