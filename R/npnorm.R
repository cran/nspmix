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
  structure(x, class= "npnorm")
}

npnorm = function(v, w=1) {
  if(class(v) == "npnorm") { v$w = v$w * w; v }
  else structure(list(v=v, w=w), class="npnorm")
}

length.npnorm = function(x) length(x$v)

weights.npnorm = function(x, beta) x$w

gridpoints.npnorm = function(x, beta, grid=100) {
  rx = range(x$v)
  breaks = pmax(ceiling(diff(rx) / (5*beta)), 5)   # number of breaks
  r = whist(x$v, x$w, breaks=breaks, freq=FALSE, plot=FALSE)
  i = r$density != 0
  i = i | c(i[-1],FALSE) | c(FALSE,i[-length(i)])  # include neighbours
  m = sum(i)
  k = pmax(ceiling(grid / m), 10)           # at least 10 in each interval
  d = r$breaks[2] - r$breaks[1]
  s = r$breaks[-length(r$breaks)][i]
  c(rx[1], rep(s, rep(k,m)) + d * (1:k-0.5)/k, rx[2])
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
      breaks = seq(min(x$v), max(x$v), len=pmin(20, kmax))
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

pnpnorm = function(x, mix, beta) {
  rowSums( outer(x, mix$pt, pnorm, sd=beta) *
          rep(mix$pr, rep(length(x), length(mix$pr))) )
}

