# ===================== #
# Discrete distribution #
# ===================== #

# pt    Points
# pr    Probabilities at the points

disc = function(pt, pr=1) {
  if (length(pt)== 0) d = list(pt=numeric(0), pr=numeric(0))
  else {
    k = max(length(pt), length(pr), na.rm=TRUE)
    pt = rep(pt, len=k)
    if(length(pr) == 0) pr = 1
    pr = rep(pr, len=k)
    d = list(pt=pt, pr=pr/sum(pr))
  }
  class(d) = "disc"
  d
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

sort.disc = function(x) {
  if( is.null(x$pt) ) return(x)
  index = order(x$pt)
  x$pt = x$pt[index]
  x$pr = x$pr[index]
  x
}

# is.unsorted.disc = function(d) is.unsorted(d$pt)

# var.disc = function(d) {
#   if(length(d$pt) <= 1) 0
#   else {
#     mu = sum(d$pr * d$pt)
#     sum( d$pr * (d$pt - mu)^2 )
#   }
# }

# Unique
# apt     Attraction points

unique.disc = function(x, prec=0) {
  if( length(x$pt) == 1 ) return(x)
  x = sort.disc(x)
  prec = rep(prec, len=2)
  if ( all(prec < 0) ) return(x)
  pt2 = pt = x$pt
  pr2 = pr = x$pr
  j  = pr <= prec[2]
  pt = pt[!j]
  pr = pr[!j]
  index = 0
  repeat {
    if( length(pt) == 0 ) break
    j = abs(pt[1] - pt) <=  prec[1]
    index = index + 1
    pt2[index] = weighted.mean(pt[j], pr[j])
    pr2[index] = sum( pr[j] )
    pt = pt[!j]
    pr = pr[!j]
  }
  disc(pt=pt2[1:index], pr=pr2[1:index])
}


