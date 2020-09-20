## weighted histogram



##'Weighted Histograms
##'
##'
##'Plots or computes the histogram with observations with
##'multiplicities/weights.
##'
##'
##'Just like \code{hist}, \code{whist} can either plot the histogram or compute
##'the values that define the histogram, by setting \code{plot} to \code{TRUE}
##'or \code{FALSE}.
##'
##'The histogram can either be the one for frequencies or density, by setting
##'\code{freq} to \code{TRUE} or \code{FALSE}.
##'
##'@param x a vector of values for which the histogram is desired.
##'@param w a vector of multiplicities/weights for the values in \code{x}.
##'@param breaks,plot,freq,xlim,ylim,xlab,ylab,main,add,col,border,lwd These
##'arguments have similar functionalities to their namesakes in function
##'\code{hist}.
##'@param ... arguments passed on to function \code{plot}.
##'@author Yong Wang <yongwang@@auckland.ac.nz>
##'@seealso \code{\link{hist}}.
##'@keywords function
##'@export whist
whist = function(x, w=1, breaks="Sturges", plot=TRUE, freq=NULL, 
                 xlim=NULL, ylim=NULL, xlab="Data", ylab=NULL, main=NULL,
                 add=FALSE, col=NULL, border=NULL, lwd=1, ...) {
  r = hist(x, breaks=breaks, plot=FALSE)
  breaks = r$breaks
  i = cut(x, breaks, include.lowest=TRUE)
  f = tapply(rep(w,len=length(i)), i, sum)            # frequency
  f[is.na(f)] = 0
  dimnames(f)[[1]] = NULL
  d = f / sum(f) / (breaks[2] - breaks[1])            # density
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
