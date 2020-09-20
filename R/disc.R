# ===================== #
# Discrete distribution #
# ===================== #

## pt    Points
## pr    Probabilities at the points

## a "disc" object always has its support points sorted in ascending order



##'Class `disc'
##'
##'
##'Class \code{disc} is used to represent an arbitrary univariate discrete
##'distribution with a finite number of support points.
##'
##'Function \code{disc} creates an object of class \code{disc}, given the
##'support points and probability values at these points.
##'
##'Function \code{print.disc} prints the discrete distribution.
##'
##'
##'@aliases disc print.disc
##'@param pt a numeric vector for support points.
##'@param pr a numeric vector for probability values at the support points.
##'@param x an object of class \code{disc}.
##'@param ... arguments passed on to function \code{print}.
##'@author Yong Wang <yongwang@@auckland.ac.nz>
##'@seealso \code{\link{cnm}}, \code{\link{cnmms}}.
##'@keywords class function
##'@examples
##'
##'
##'(d = disc(pt=c(0,4), pr=c(0.3,0.7)))
##'
##'@usage
##'disc(pt, pr=1)
##'\method{print}{disc}(x, ...)
##' 
##'@export disc
##'@export print.disc

disc = function(pt, pr=1) {
  if (length(pt)== 0) d = list(pt=numeric(0), pr=numeric(0))
  else {
    k = max(length(pt), length(pr), na.rm=TRUE)
    pt = rep(pt, len=k)
    if(length(pr) == 0) pr = 1
    pr = rep(pr, len=k)
    o = order(pt)
    d = list(pt=pt[o], pr=pr[o]/sum(pr))
  }
  structure(d, class="disc")
}

# is.null.disc = function (d) is.null(d$pt)

# is.disc = function (d) any(class(d) == "disc")

print.disc = function (x, ...) {
  b = cbind(x$pt, x$pr)
  dimnames(b) = list(NULL, c("pt", "pr"))
  print(b, ...)
}

plot.disc = function (x, type=c("pdf","cdf"), add=FALSE, col="blue", lwd=1, ylim,
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

# Sort

## sort.disc = function(x) {
##   if( is.null(x$pt) ) return(x)
##   o = order(x$pt)
##   x$pt = x$pt[o]
##   x$pr = x$pr[o]
##   x
## }

# is.unsorted.disc = function(d) is.unsorted(d$pt)

# var.disc = function(d) {
#   if(length(d$pt) <= 1) 0
#   else {
#     mu = sum(d$pr * d$pt)
#     sum( d$pr * (d$pt - mu)^2 )
#   }
# }

## unique: remove support points with zero mass and collapse similar support
## points

## x    object of class "disc"

## Detail: Find similar pairs and collapse them

unique.disc = function(x, prec=0) {
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


## OLD.unique.disc = function(x, prec=0) {
##   if( length(x$pt) == 1 ) return(x)
## #   x = sort.disc(x)
##   prec = rep(prec, len=2)
##   ## if ( all(prec < 0) ) return(x)
##   pt2 = pt = x$pt
##   pr2 = pr = x$pr
##   j  = pr <= prec[2]
##   pt = pt[!j]
##   pr = pr[!j]
##   index = 0
##   repeat {
##     if( length(pt) == 0 ) break
##     j = abs(pt[1] - pt) <=  prec[1]
##     index = index + 1
##     pt2[index] = weighted.mean(pt[j], pr[j])
##     pr2[index] = sum( pr[j] )
##     pt = pt[!j]
##     pr = pr[!j]
##   }
##   disc(pt=pt2[1:index], pr=pr2[1:index])
## }


