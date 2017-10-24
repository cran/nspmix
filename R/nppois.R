# ============================== #
# Nonparametric Poisson Mixtures #
# ============================== #

# Class nppois #

# mix     Object of "disc"

# rnppois = function(n=100, lambda=c(1,4), pr=c(.3,.7), mix) {
#   if(missing(mix)) mix = disc(lambda, pr)
#   if (n == 0) return(numeric(0))
#   k = length(mix$pt)
#   suppressWarnings(i <- sample.int(k, n, prob = mix$pr, replace = TRUE))
#   x = rpois(n, mix$pt[i])
#   w = table(x)
#   r = list(v=as.integer(names(w)), w=as.integer(w))
#   class(r) = "nppois"
#   r
# }

rnppois = function(n, lambda=1, pr=1) {
  mix = disc(lambda, pr)
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

weights.nppois = function(x, beta) x$w

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
    xlim = floor(c(pmax(ptr[1] - 4 * sqrt(ptr[1]), 0),
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


