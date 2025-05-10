# ===================== #
# Discrete distribution #
# ===================== #

## pt    Points
## pr    Probabilities at the points

## a "disc" object always has its support points sorted in ascending order



##' Class \code{disc}
##' 
##' Class \code{disc} is used to represent an arbitrary univariate
##' discrete distribution with a finite number of support points.
##' 
##' Function \code{disc} creates an object of class \code{disc}, given
##' the support points and probability values at these points.
##' 
##' Function \code{print.disc} prints the discrete distribution.
##' 
##' @param pt a numeric vector for support points.
##' @param pr a numeric vector for probability values at the support
##'   points.
##' @param sort =TRUE, by default. If TRUE, support points are sorted
##'   (in increasing order).
##' @param collapse =TRUE, by default. If TRUE, identical support
##'   points are collapsed, with their masses aggregated.
##' 
##' @author Yong Wang <yongwang@@auckland.ac.nz>
##' 
##' @seealso \code{\link{cnm}}, \code{\link{cnmms}}.
##' 
##' @examples
##' 
##' (d = disc(pt=c(0,4), pr=c(0.3,0.7)))
##' 
##' @usage
##' 
##' disc(pt, pr=1, sort=TRUE, collapse=FALSE)
##'  
##' @export

disc = function(pt, pr=1, sort=TRUE, collapse=FALSE) {
  if(collapse) sort = TRUE
  if (length(pt) == 0) d = list(pt=numeric(0), pr=numeric(0))
  else {
    k = max(length(pt), length(pr), na.rm=TRUE)
    pt = rep(pt, len=k)
    if(length(pr) == 0) pr = 1
    pr = rep(pr, len=k)
    if(sort) {
      o = order(pt)
      d = list(pt=pt[o], pr=pr[o]/sum(pr))
    }
    else d = list(pt=pt, pr=pr/sum(pr))
  }
  if(collapse) {
    j = d$pr != 0
    pt = d$pt[j]
    pr = d$pr[j]
    jd0 = which(! duplicated(pt))    # which distinct
    jd1 <- c(jd0, length(pr)+2)
    cumpr <- c(0,cumsum(pr),1)[jd1]
    pr = diff(cumpr)
    pt = pt[jd0]
    structure(list(pt=pt, pr=pr), class="disc")
  }
  else structure(d, class="disc")
}

##' @title Prints a discrete distribution function
##' 
##' @description Class \code{disc} is used to represent an arbitrary
##'   univariate discrete distribution with a finite number of support
##'   points.
##' 
##' Function \code{disc} creates an object of class \code{disc}, given
##'  the support points and probability values at these points.
##' 
##' Function \code{print.disc} prints the discrete distribution.
##' 
##' @param x an object of class \code{disc}.
##' @param ... arguments passed on to function \code{print}.
##' 
##' @author Yong Wang <yongwang@@auckland.ac.nz>
##' 
##' @seealso \code{\link{cnm}}, \code{\link{cnmms}}.
##' 
##' @examples
##' 
##' (d = disc(pt=c(0,4), pr=c(0.3,0.7)))
##'  
##' @export

print.disc = function (x, ...) {
  b = cbind(x$pt, x$pr)
  dimnames(b) = list(NULL, c("pt", "pr"))
  print(b, ...)
}


##' @title Plot a discrete distribution function
##' 
##' @description Class \code{disc} is used to represent an arbitrary
##'   univariate discrete distribution with a finite number of support
##'   points.
##' 
##' Function \code{disc} creates an object of class \code{disc}, given
##'  the support points and probability values at these points.
##' 
##' Function \code{plot.disc} plots the discrete distribution.
##' 
##' @param x an object of class \code{disc}.
##' @param type plot its pdf or cdf.
##' @param add add the plot or not.
##' @param col colour to be used.
##' @param lwd,ylim,xlab,ylab graphical parameters.
##' @param ... arguments passed on to function \code{plot}.
##' 
##' @author Yong Wang <yongwang@@auckland.ac.nz>
##' 
##' @seealso \code{\link{disc}}, \code{\link{cnm}}, \code{\link{cnmms}}.
##' 
##' @examples
##' 
##' plot(disc(pt=c(0,4), pr=c(0.3,0.7)))
##' plot(disc(rnorm(5), 1:5))
##' for(i in 1:5)
##'    plot(disc(rnorm(5), 1:5), type="cdf", add=(i>1), xlim=c(-3,3))
##' 
##' @export

plot.disc = function (x, type=c("pdf","cdf"), add=FALSE, col=4, lwd=1,
                      ylim,
                      xlab="", ylab="Probability", ...) {
  type = match.arg(type)
  if(missing(ylim)) ylim = switch(type, pdf=c(0, max(x$pr)), cdf=c(0,1))
  type2 = switch(type, pdf="h", cdf="s")
  x1 = switch(type, pdf=x$pt, cdf=c(x$pt[1], x$pt))
  y1 = switch(type, pdf=x$pr, cdf=c(0,cumsum(x$pr)))
  if(add) lines(x1, y1, type=type2, col=col, lwd=lwd, ...)
  else {
    plot(x1, y1, type=type2, col=col, lwd=lwd, ylim=ylim,
         xlab=xlab, ylab=ylab, ...)
    if(type == "cdf") abline(h=c(0,1), col="grey")
  }
}

##' @export

sort.disc = function(x, decreasing=FALSE, ...) {
  if( is.null(x$pt) ) return(x)
  o = order(x$pt, decreasing=decreasing)
  x$pt = x$pt[o]
  x$pr = x$pr[o]
  x
}

## unique: remove support points with zero mass and collapse similar support
## points

## x    object of class "disc"

## Detail: Find similar pairs and collapse them

##' @export 

unique.disc = function(x, incomparables=FALSE, prec=0, ...) {
  if( length(x$pt) == 1 ) return(x)
  prec = rep(prec, len=2)
  j0  = x$pr <= prec[2]
  pt = x$pt[!j0]
  pr = x$pr[!j0]
  repeat {
    i = diff(pt) <= prec[1]                # will use all odd indexes of TRUEs
    if(sum(i) == 0) break
    noti = c(!i, TRUE)                              # not i
    even = seq(2, by=2, len=length(i))
    i[even] = i[even] & noti[even-1] & noti[even+1] # use evens if only isolated
    ind = which(i)
    pt[ind] = (pt[ind] * pr[ind] + pt[ind+1] * pr[ind+1]) /
      (pr2 <- pr[ind] + pr[ind+1])                  # collapse neighbouring pairs
    pr[ind] = pr2
    pt = pt[-(ind+1)]                               # remove the second of a pair
    pr = pr[-(ind+1)]
  }
  disc(pt=pt, pr=pr)
}



