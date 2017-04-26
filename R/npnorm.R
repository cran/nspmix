# ============================= #
# Nonparametric normal mixtures #
# ============================= #

# ============ #
# Class npnorm #
# ============ #

# mix     Object of "disc"

rnpnorm = function(n, mu=0, pr=1, sd=1) {
  mix = disc(mu, pr)
  if (n == 0) return(numeric(0))
  k = length(mix$pt)
  suppressWarnings(i <- sample.int(k, n, prob = mix$pr, replace = TRUE))
  x = list(v=mix$pt[i] + rnorm(n, sd=sd), w=1)
  class(x) = "npnorm"
  x
}

npnorm = function(v, w=1) {
  if(class(v) == "npnorm") {
    v$w = v$w * w
    v
  }
  else {
    x = list(v=v, w=w)
    class(x) = "npnorm"
    x
  }
}

length.npnorm = function(x) length(x$v)

weights.npnorm = function(x, beta) x$w

gridpoints.npnorm = function(x, beta, grid=100) {
  breaks = pmax(ceiling(diff(range(x$v)) / (5*beta)), 5)   # number of breaks
  r = whist(x$v, x$w, breaks=breaks, freq=FALSE, plot=FALSE)
  i = r$density != 0
  i = i | c(i[-1],FALSE) | c(FALSE,i[-length(i)])  # include neighbours
  m = sum(i)
  k = pmax(ceiling(grid / m), 10)
  d = r$breaks[2] - r$breaks[1]
  s = r$breaks[-length(r$breaks)][i]
  unique(rep(s, rep(k,m)) + d * 0:(k-1)/(k-1))
}

# beta        Standard deviation, with default value = 1
# mix         Discrete mixing distribution

initial.npnorm = function(x, beta=NULL, mix=NULL, kmax=NULL) {
  if(is.null(beta)) beta = 1
  if(is.null(mix) || is.null(mix$pt)) {
    if(is.null(kmax)) breaks = pmax(ceiling(diff(range(x$v)) / (5*beta)), 5)
    else {
      if(kmax == 1)
        return(list(beta=beta,
                    mix=disc(sum(x$v*x$w) / sum(rep(x$w,len=length(x$v))))))
      breaks = seq(min(x$v), max(x$v), len=3)
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
  x.pt = x$v - rep(pt, each=length(x$v))
  dim(x.pt) = c(length(x$v), length(pt))
  beta2 = beta * beta
  if(which[1] == 1)
    dl$ld = - .5 * (log(2*base::pi*beta2) + x.pt^2 / beta2)
  if(which[2] == 1) {
    dl$db = (-1 + x.pt^2 / beta2) / beta
    dim(dl$db) = c(length(x$v), length(pt), 1)
  }
  if(which[3] == 1)
    dl$dt = x.pt / beta2
  dl
}

plot.npnorm = function(x, mix, beta, breaks=NULL, col="red", add=FALSE,
                       border.col=NULL, border.lwd=1, fill="lightgrey",
                       components=TRUE, main, lwd=2,
                       lty=1, xlab="Data", ylab="Density", ...) {
  m = length(mix$pt)
  z = sort(unique(round(mix$pt + rep(beta*seq(-3,3,len=100), each=m),
                        ceiling(-log10(beta/20)))))
  nz = length(z)
  dj = outer(z, mix$pt, dnorm, sd=beta) * rep(mix$pr, each=nz)
  d = rowSums(dj)
  if(missing(main))
    main = substitute("npnorm (" * sigma ~ "=" ~ xxx * ")",
                      list(xxx=signif(beta,3)))
  if(add || missing(x) || length(x) == 0) {
    if(add) lines(z, d, col=col, lwd=lwd, lty=lty)
    else {
      plot(0, 0, type="n", xlim=range(z), ylim=range(d),
           xlab=xlab, ylab=ylab, frame.plot=FALSE,
           main=main, ...)
      lines(range(z), rep(0,2), col="darkgrey")
      lines(z, d, col=col, lwd=lwd, lty=lty) 
    }
  }
  else {
    if(is.null(breaks)) breaks = round(sqrt(length(x$v)))
    whist(x$v, x$w, breaks=breaks, freq=FALSE, 
          xlab=xlab, ylab=ylab, main=main, col=fill, border=border.col, 
          lwd=border.lwd, ylim=range(d), ...)
    lines(z, d, col=col, lwd=lwd, lty=lty)
  }
  if(components) {
    points(mix$pt, rep(0,length(mix$pt)), col=col)
    segments(mix$pt, rep(0,m), y1=mix$pr*max(d), col=col, lwd=3)
  }
}

## Additional functions

pnpnorm = function(x, mix, beta) {
  rowSums( outer(x, mix$pt, pnorm, sd=beta) *
          rep(mix$pr, rep(length(x), length(mix$pr))) )
}

