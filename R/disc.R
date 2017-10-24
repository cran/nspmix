# ===================== #
# Discrete distribution #
# ===================== #

## pt    Points
## pr    Probabilities at the points

## a "disc" object must have its support points sorted in ascending order

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

plot.disc = function (x, ...) {
  ylim = c(0, max(x$pr))
  plot(x$pt, x$pr, type="h", col="blue", lwd=2, ylim=ylim,
       xlab="", ylab="Probability", ...)
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


