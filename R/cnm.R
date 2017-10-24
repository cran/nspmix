# ==================================================== #
# Maximum likelihood computation for nonparametric and #
# semiparametric mixture model                         # 
# ==================================================== #

valid = function(x, beta, mix) UseMethod("valid")
suppspace = function(x, beta) UseMethod("suppspace")
gridpoints = function(x, beta, grid) UseMethod("gridpoints")
initial = function(x, beta, mix, kmax) UseMethod("initial")
logd = function(x, beta, pt, which) UseMethod("logd")
llex = function(x, beta, mix) UseMethod("llex")
llexdb = function(x, beta, mix) UseMethod("llexdb")
# Other generic functions needed: length, weights
# print?

valid.default = function(x, beta, mix) TRUE
suppspace.default = function(x, beta) c(-Inf, Inf)
llex.default = function(x, beta, mix) 0
llexdb.default = function(x, beta, mix) 0

# Compute the NPMLE or mixing proportions only

# model   = prop,    computes the mixing proportions
#         = npmle,   computes the NPMLE

cnm = function(x, init=NULL, model=c("npmle","proportions"), maxit=100, tol=1e-6,
               grid=100, plot=c("null", "gradient", "probability"),
               verbose=0) {
  if(is.null(class(x))) stop("Data belongs to no class")
  plot = match.arg(plot)
  model = match.arg(model)
  k = length(x)
  if(model == "proportions") {
    if(is.null(init$mix))
      stop("init$mix must be provided for model=='proportions'")
    pt0 = sort(unique(init$mix$pt))
    m0 = length(pt0)
    ind = unique(round(seq(1, length(pt0), len=100)))
    init$mix = disc(pt0[ind])
  }
  init = initial.snpmle(x, init)
  beta = init$beta
  nb = length(beta)
  mix = init$mix
  attr(mix, "ll") = init$ll
  convergence = 1
  gp = gridpoints(x, beta, grid)          # does not depend on beta
  w = weights(x, beta)
  wr = sqrt(w)
  if(maxit == 0) 
    return( structure(list(family=class(x)[1], mix=mix, beta=beta,
                           num.iterations=0, ll=ll[1],
                           max.gradient=max(maxgrad(x, beta, mix, tol=-5)$grad),
                           convergence=1),
                      class = "nspmix")  )
  for(i in 1:maxit) {
    ll = attr(mix, "ll")
    dmix = attr(ll, "dmix")
    ma = attr(ll, "logd") - log(dmix)
    switch(plot, 
           "gradient" = plotgrad(x, mix, beta, pch=1), 
           "probability" = plot(x, mix, beta) )
    mix1 = mix
    if(model == "npmle") {
      g = maxgrad(x, beta, mix, grid=gp, tol=-5)
      if(plot=="gradient") points(g$pt, g$grad, col="red", pch=20)
      gradient = max(g$grad)
      mix = disc(c(mix$pt,g$pt), c(mix$pr,rep(0,length(g$pt))))
    }
    else {
      gpt0 = grad(pt0, x, beta, mix)$d0
      ind = indx(pt0, mix$pt)
      ind[ind == 0] = 1
      s = split(gpt0, ind)
      imax = sapply(s, which.max)
      gradient = max(gpt0[imax])
      imax = imax[imax > 1]
      csum = c(0,cumsum(sapply(s, length)))
      ipt = csum[as.numeric(names(imax))] + imax
      mix = disc(c(mix$pt, pt0[ipt]), c(mix$pr, double(length(ipt))))
    }
    lpt = logd(x, beta, mix$pt, which=c(1,0,0))$ld
    D = pmin(exp(lpt - ma), 1e100)
    if(verbose > 0) 
      print.snpmle(verbose, attr(mix, "ll")[1], mix1, beta, gradient, i-1)
    r = hcnm(D, mix$pr, w=w, maxit=3)
    mix$pr = r$pf
    # attr(mix, "ll") = r$ll + sum(w * ma)
    mix = collapse.snpmle(mix, beta, x, tol=0)
    if( attr(mix, "ll")[1] <= attr(mix1, "ll")[1] + tol ) {convergence = 0; break}
  }
  if(verbose > 0)
    print.snpmle(verbose, attr(mix, "ll")[1], mix, beta, gradient, i)
  if(model == "npmle")
    mix = collapse.snpmle(mix, beta, x,
                          tol=max(tol*.1, abs(attr(mix, "ll")[1])*1e-16))
  ll = attr(mix,"ll")[1]
  attr(mix,"ll") = NULL
  structure(list(family=class(x)[1], num.iterations=i, max.gradient=gradient,
                 convergence=convergence, ll=ll, mix=mix, beta=beta),
            class="nspmix")
}

## Hierarchical CNM

## S = D / P

hcnm = function(D, p0, w=1, maxit=3, tol=1e-6,
                blockpar=NULL, recurs.maxit=2, depth=1) {
  nx = nrow(D)
  w = rep(w, length=nx)
  wr = sqrt(w)
  n = sum(w)
  converge = FALSE
  m = ncol(D)
  j = 1:m
  nblocks = 1
  ## maxdepth = depth
  if(length(p <- p0) != m) stop("Argument 'p0' is the wrong length.")
  p = p / sum(p)
  P = drop(D %*% p)
  ll = sum(w * log(P))     # log-likelihood relative to D
  evenstep = FALSE
  
  for(iter in 1:maxit) {
    p.old = p
    ll.old = ll
    S = D / pmax(P, 1e-100)
    g = colSums(w * S)
    dmax = max(g) - n
    sj = m
    if(is.null(blockpar) || is.na(blockpar))
      iter.blockpar = ifelse(sj < 30, 0,
                             1 - log(max(20,10*log(sj/100)))/log(sj))
    else iter.blockpar = blockpar
    if(iter.blockpar==0 | sj < 30) {
      nblocks = 1
      BW = matrix(1, nrow=sj, ncol=1)
    }
    else {
      nblocks = max(1, if(iter.blockpar>1) round(sj/iter.blockpar)
                       else floor(min(sj/2, sj^iter.blockpar)))
      i = seq(0, nblocks, length=sj+1)[-1]
      if(evenstep) {
        nblocks = nblocks + 1
        BW = outer(round(i)+1, 1:nblocks, "==")
      }
      else BW = outer(ceiling(i), 1:nblocks, "==")
      storage.mode(BW) = "numeric"
    }

    for(block in 1:nblocks) {
      jj = logical(m)
      jj[j] = BW[,block] > 0
      sjj = sum(jj)
      if (sjj > 1 && (delta <- sum(p.old[jj])) > 0) {
        Sj = S[,jj,drop=FALSE]
        res = pnnls(wr * Sj, wr * drop(Sj %*% p.old[jj]) + wr, sum=delta)
        ## if (res$mode > 1) warning("Problem in pnnls(a,b)")
        p[jj] = p[jj] +  BW[jj[j],block] *
          (res$x * (delta / sum(res$x)) - p.old[jj])
      }
    }
    
    p.gap = p - p.old 
    ll.rise.gap = sum(g * p.gap) 
    alpha = 1
    p.alpha = p
    ll.rise.alpha = ll.rise.gap
    repeat {
      P = drop(D %*% p.alpha)
      ll = sum(w * log(P))
      if(ll >= ll.old && ll + ll.rise.alpha <= ll.old) {
        p = p.alpha
        converge = TRUE
        break
      }
      if(ll > ll.old && ll >= ll.old + ll.rise.alpha * .33) {
        p = p.alpha 
        break
      }
      if((alpha <- alpha * 0.5) < 1e-10) {
        p = p.old
        P = drop(D %*% p)
        ll = ll.old
        converge = TRUE
        break
      }
      p.alpha = p.old + alpha * p.gap
      ll.rise.alpha = alpha * ll.rise.gap
    }
    if(converge) break

    if (nblocks > 1) {
      Q = sweep(BW,1,p[j],"*")
      q = colSums(Q) 
      Q = sweep(D[,j] %*% Q, 2, q, "/")  
      if( any(q == 0) ) {
        ## warning("A block has zero probability!")
      }
      else {
        res = hcnm(D=Q, p0=q, w=w, blockpar=iter.blockpar,
                   maxit=recurs.maxit, recurs.maxit=recurs.maxit,
                   depth=depth+1)
        ## maxdepth = max(maxdepth, res$maxdepth)
        if (res$ll > ll) {
          p[j] = p[j] * (BW %*% (res$pf / q))
          P = drop(D %*% p)
          ll = sum(w * log(P))
        }
      }
    }
    if(iter > 2) if( ll <= ll.old + tol ) {converge=TRUE; break}
    evenstep = !evenstep
  }
  list(pf=p, convergence=converge, method="hcnm", ll=ll,
       maxgrad=max(crossprod(w/P, D))-n, numiter=iter)
}

# Updates mix using CNM and then (pi, theta, beta) using BFGS.

# Modifying the support set

cnmms = function(x, init=NULL, maxit=1000,
      model=c("spmle","npmle"), tol=1e-6, grid=100, kmax=Inf,
      plot=c("null", "gradient", "probability"), verbose=0) {
  if(is.null(class(x))) stop("Data belongs to no class")
  plot = match.arg(plot)
  model = match.arg(model)
  init = initial.snpmle(x, init, kmax=kmax)
  beta = init$beta
  if(is.null(beta) || is.na(beta)) model = "npmle"
  nb = length(beta)
  mix = init$mix
  attr(mix, "ll") = init$ll
  gradient = "Not computed"
  switch(plot,
         "gradient" = plotgrad(x, mix, beta, pch=1),
         "probability" = plot(x, mix, beta) )
  ## ll = -Inf       # logLik.snpmle(x, beta, mix, attr=FALSE)
  if(maxit == 0) 
    return( structure(list(family=class(x)[1], mix=mix, beta=beta,
                           num.iterations=0, ll=attr(mix, "ll")[1],
                           max.gradient=max(maxgrad(x, beta, mix, tol=-5)$grad),
                           convergence=1), 
                      class="nspmix") )
  convergence = 1
  for(i in 1:maxit) {
    mix1 = mix
    if(length(mix$pt) < kmax) {
      gp = gridpoints(x, beta, grid)
      g = maxgrad(x, beta, mix, grid=gp, tol=-Inf)
      gradient = max(g$grad)
      kpt = min(kmax - length(mix$pt), length(g$pt))
      jpt = order(g$grad, decreasing=TRUE)
      mix = disc(c(mix$pt,g$pt[jpt][1:kpt]), c(mix$pr,rep(0,kpt)))
      if(plot=="gradient") points(g$pt, g$grad, col="red", pch=20)
    }
    if(verbose > 0)
      print.snpmle(verbose, attr(mix1, "ll")[1], mix1, beta, gradient, i-1)
    lpt = logd(x, beta, mix$pt, which=c(1,0,0))$ld
    ## ll1 = attr(mix1, "ll")
    ## logd.mix1 = log(attr(ll1,"dmix")) + attr(ll1,"ma")
    ## S = pmin(exp(lpt - logd.mix1), 1e100)
    S = pmin(exp(lpt - attr(attr(mix1,"ll"),"logd")), 1e100)
    w = weights(x, beta)
    wr = sqrt(w)
    grad.support = colSums(w * (S - 1))
    r = pnnls(wr * S, wr * 2, sum=1)
    sol = r$x / sum(r$x)
    r = lsch(mix, beta, disc(mix$pt,sol), beta, x, which=c(1,0,0))
    mix = r$mix
    attr(mix, "ll") = logLik.snpmle(x, beta, mix)
    mix = collapse.snpmle(mix, beta, x)
    if(max(grad.support) < 1e10) {    # 1e20
      r = switch(model,
        spmle = bfgs(mix, beta, x, which=c(1,1,1)),
        npmle = bfgs(mix, beta, x, which=c(1,1,0))   )
      beta = r$beta
      mix = collapse.snpmle(r$mix, beta, x)
    }
    switch(plot,
           "gradient" = plotgrad(x, mix, beta, pch=1),
           "probability" = plot(x, mix, beta) )
    if(attr(mix,"ll")[1] <= attr(mix1,"ll")[1] + tol) {convergence = 0; break}
  }
  if(verbose > 0)
    print.snpmle(verbose, attr(mix,"ll")[1], mix, beta, gradient, i)
  mix = collapse.snpmle(mix, beta, x, tol=tol) # , tol=1e-12)
  m = length(mix$pt)
  if(m < length(r$mix$pt)) {
    d = dll.snpmle(x, mix, beta, which=c(0,1,1,1))
    grad = c(d$dp, d$dt, d$db)
    names(grad) = c(paste0("pr",1:m), paste0("pt",1:m), 
           if(is.null(d$db)) NULL else paste0("beta",1:length(beta)) )
  }
  else grad = r$grad
  grad[1:m] = grad[1:m] - sum(rep(weights(x, beta), len=length(x)))
  ll = attr(mix,"ll")[1]
  attr(mix,"ll") = NULL
  structure(list(family=class(x)[1], num.iterations=i, grad=grad,
                 convergence=convergence, ll=ll, mix=mix, beta=beta),
            class="nspmix")
}

# Use the BFGS method to maximize the profile likelihood function

cnmpl = function(x, init=NULL, tol=1e-6,
         tol.npmle=tol*1e-4, grid=100, maxit=1000,
         plot=c("null", "gradient", "probability"), verbose=0) {
  if(is.null(class(x))) stop("Data belongs to no class")
  plot = match.arg(plot)
  init = initial.snpmle(x, init)
  beta = init$beta
  mix = init$mix
  nb = length(beta)
  if(is.null(beta) || is.na(beta))
    stop("No initial value for beta. Use cnm() for NPMLE problems.")
  switch(plot,
         "gradient" = plotgrad(x, mix, beta, pch=1),
         "probability" = plot(x, mix, beta) )
  r = pll(beta, mix, x, tol=tol.npmle, grid=grid)
  D = diag(-1, nrow=nb)
  convergence = 1
  num.npmle = 1
  if(maxit == 0) {
    r = list(family=class(x)[1], mix=mix, beta=beta, num.iterations=0,
        ll=logLik.snpmle(x, beta, mix),
        max.gradient=max(maxgrad(x, beta, mix, tol=-5)$grad), convergence=1)
    class(r) = "nspmix"
    return(r)
  }
  for(i in 1:maxit) {
    old.r = r
    beta2 = drop(r$beta - D %*% r$grad)
    r = lsch.pll(old.r, beta2, x, tol=tol, tol.npmle=tol.npmle, 
      grid=grid, brkt=TRUE)
    switch(plot,
           "gradient" = plotgrad(x, r$mix, r$beta, pch=1),
           "probability" = plot(x, r$mix, r$beta) )
    num.npmle = num.npmle + r$num.npmle
    if( r$ll <= old.r$ll + tol ) {convergence = 0; break}
    g = r$grad - old.r$grad
    d = r$beta - old.r$beta
    dg = sum(d * g)
    if( dg < 0 ) {
      nd = length(d)
      dod = d * rep(d, rep(nd, nd))
      dim(dod) = c(nd,nd)
      dog = d * rep(g, rep(nd, nd))
      dim(dog) = c(nd,nd)
      dogD = dog %*% D
      D = D + (1 + drop(t(g) %*% D %*% g) / dg) * dod / dg -
        (dogD + t(dogD)) / dg
    }
    if(verbose > 0) print.snpmle(verbose, r$ll[1], r$mix, r$beta, r$max.grad, i)
  }
  structure(list(family=class(x)[1], num.iterations=i, grad=r$grad,
                  max.gradient=r$max.gradient, # num.cnm.calls=num.npmle,
                  convergence=convergence, ll=r$ll, mix=r$mix, beta=r$beta),
            class="nspmix")
}

# Updates mix using one iteration of CNM and then updates beta usin BFGS.
# It is almost the standard cyclic method, except only one iteration of CNM
# is used.

cnmap = function(x, init=NULL, maxit=1000, tol=1e-6, grid=100, 
      plot=c("null", "gradient"), verbose=0) {
  if(is.null(class(x))) stop("Data belongs to no class")
  plot = match.arg(plot)
  k = length(x)
  init = initial.snpmle(x, init)
  beta = init$beta
  if(is.null(beta) || is.na(beta))
    stop("No initial value for beta. Use cnm() for NPMLE problems.")
  nb = length(beta)
  mix = init$mix
  attr(mix, "ll") = init$ll
  convergence = 1
  for(i in 1:maxit) {
    mix1 = mix
    switch(plot, "gradient" = plotgrad(x, mix, beta) )
    gp = gridpoints(x, beta, grid)
    g = maxgrad(x, beta, mix, grid=gp, tol=-Inf)
    gradient = max(g$grad)
    if(verbose > 0)
      print.snpmle(verbose, attr(mix, "ll")[1], mix, beta, gradient, i-1)
    mix = disc(c(mix$pt,g$pt), c(mix$pr,rep(0,length(g$pt))))
    lpt = logd(x, beta, mix$pt, which=c(1,0,0))$ld
    ## ll1 = attr(mix1, "ll")
    # logd.mix1 = log(attr(ll1,"dmix")) + attr(ll1,"ma")
    ## S = pmin(exp(lpt - logd.mix1), 1e100)
    S = pmin(exp(lpt - attr(attr(mix1,"ll"),"logd")), 1e100)
    wr = sqrt(weights(x, beta))
    r = pnnls(wr * S, wr * 2, sum=1)
    sol = r$x / sum(r$x)
    r = lsch(mix, beta, disc(mix$pt,sol), beta, x, which=c(1,0,0))
    mix = collapse.snpmle(r$mix, beta, x)
    r = bfgs(mix, beta, x, which=c(0,0,1))
    beta = r$beta
    mix = collapse.snpmle(r$mix, beta, x)
    if( attr(mix, "ll")[1] <= attr(mix1, "ll")[1] + tol ) {convergence = 0; break} 
  }
  mix = collapse.snpmle(mix, beta, x, tol=tol) # , tol=1e-12)
  m = length(mix$pt)
  d = dll.snpmle(x, mix, beta, which=c(0,1,1,1))
  grad = c(d$dp, d$dt, d$db)
  names(grad) = c(paste0("pr",1:m), paste0("pt",1:m), 
                  if(is.null(d$db)) NULL else paste0("beta",1:length(beta)) )
  grad[1:m] = grad[1:m] - sum(rep(weights(x, beta), len=length(x)))
  ll = attr(mix,"ll")[1]
  attr(mix,"ll") = NULL
  structure(list(family=class(x)[1], num.iterations=i, grad=grad,
                 max.gradient=gradient, convergence=convergence,
                 ll=ll, mix=r$mix, beta=r$beta),
            class="nspmix")
}

# The BFGS quasi-Newton method for updating pi, theta and beta

# The default negative identity matrix for D may be badly scaled for a
# given problem. It can result in failures to convergence. 
# In such a case, it'd better use cnmpl().

# beta      a real- or vector-valued parameter that has no constraints 

# which     which of pi, theta and beta are to be updated

bfgs = function(mix, beta, x, tol=1e-15, maxit=1000, which=c(1,1,1), D=NULL) {
  k1 = if(which[1]) length(mix$pr) else 0
  k2 = if(which[2]) length(mix$pt) else 0
  k3 = if(which[3]) length(beta) else 0
  if( sum(which) == 0 ) stop("No parameter specified for updating in bfgs()")
  dl = dll.snpmle(x, mix, beta, which=c(1,which), ind=TRUE)
  w = weights(x, beta)
  ll = llex(x, beta, mix) + sum(w * dl$ll)
  grad = c(if(which[1]) colSums(w * dl$dp) else NULL,
    if(which[2]) colSums(w * dl$dt) else NULL,
    if(which[3]) llexdb(x,beta,mix) + colSums(w * dl$db) else NULL)
  r = list(mix=mix, beta=beta, ll=ll, grad=grad, convergence=1)
  par = c(if(which[1]) r$mix$pr else NULL,
    if(which[2]) r$mix$pt else NULL,
    if(which[3]) r$beta else NULL)
  # if(is.null(D)) D = diag(c(rep(-1,k1+k2),rep(-1,k3)))
  if(is.null(D)) D = diag(-1, nrow=k1+k2+k3)
  else if(nrow(D) != k1+k2+k3) stop("D provided has incompatible dimensions")
  # D = - solve(crossprod(sqrt(w) * cbind(dl$dp, dl$dt, dl$db)))   # Gauss-Newton
  # D = - diag(1/colSums(w * cbind(dl$dp, dl$dt, dl$db)^2), nrow=k1+k2+k3)
  # Note: it may also depend on llexdb
  bs = suppspace(x, beta)
  for(i in 1:maxit) {
    old.r = r
    old.par = par
    d1 = - drop(D %*% r$grad)
    if(which[1]) {      # correction for unity restriction (Wu, 1978)
      D1 = D[,1:k1,drop=FALSE]
      D11 = D1[1:k1,,drop=FALSE]
      lambda = - sum(r$grad * rowSums(D1)) / sum(D11)
      d1 = d1 - lambda * rowSums(D1)
    }
    par2 = par + d1

    if(which[2]) {
      theta1 = par[k1+1:k2]
      theta2 = par2[k1+1:k2]
      dtheta = d1[k1+1:k2]
      if(any(theta2 < bs[1], theta2 > bs[2])) { # New support points outside
        out = pmax(c(bs[1]-theta2, theta2-bs[2]), 0)
        j = out > 0 & (theta1 == bs[1] | theta1 == bs[2])
        if(any(j)) {      # When old ones on boundary
          jj = (which(j)-1) %% k2 + 1
          D[k1+jj,] = D[,k1+jj] = 0    # Remove new ones outside 
          d1 = - drop(D %*% r$grad)    # Re-calculate
          if(which[1]) {
            D1 = D[,1:k1,drop=FALSE]
            D11 = D1[1:k1,,drop=FALSE]
            lambda = - sum(r$grad * rowSums(D1)) / sum(D11)
            d1 = d1 - lambda * rowSums(D1)
          }
          par2 = par + d1
        }
        else {         # Backtrack to boundary when all old ones inside
          ratio = out / abs(dtheta)
          ratio[dtheta==0] = 0
          jmax = which.max(ratio)
          alpha = 1 - ratio[jmax]
          d1 = alpha * d1
          par2 = par + d1
          if(jmax <= k2) par2[k1 + jmax] = bs[1]
          else par2[k1 + jmax - k2] = bs[2]
        }
      }
    }
    if(which[1]) {
      if(any(par2[1:k1] < 0)) {
        step = d1[1:k1]
        ratio = pmax(-par2[1:k1], 0) / abs(step)
        jmax = which.max(ratio)
        alpha = 1 - ratio[jmax]
        d1 = alpha * d1
        par2 = par + d1
        par2[jmax] = 0
      }
    }

    mix2 = disc(if(which[2]) par2[k1+1:k2] else r$mix$pt,
      if(which[1]) par2[1:k1] else r$mix$pr)
    beta2 = if(which[3]) par2[k1+k2+1:k3] else r$beta
    r = lsch(r$mix, r$beta, mix2, beta2, x, which=which, brkt=TRUE)
    if( any(r$mix$pr == 0) ) {
      j = r$mix$pr == 0
      r$mix = disc(r$mix$pt[!j], r$mix$pr[!j])
      j2 = which(j)
      if(which[2]) j2 = c(j2,k1+j2)
      D = D[-j2,-j2]
      return( bfgs(r$mix, r$beta, x, tol, maxit, which, D=D) )
    }
    if( r$conv != 0 ) {
      # cat("lsch() failed in bfgs(): convergence =", r$conv, "\n");
      break}
    # if(r$ll >= old.r$ll - tol * abs(r$ll) && r$ll <= old.r$ll + tol * abs(r$ll))
    if(r$ll >= old.r$ll - tol && r$ll <= old.r$ll + tol)
      {convergence = 0; break}
    g = r$grad - old.r$grad
    par = c(if(which[1]) r$mix$pr else NULL,
      if(which[2]) r$mix$pt else NULL,
      if(which[3]) r$beta else NULL)
    d = par - old.par
    dg = sum(d * g)
    if( dg < 0 ) {
      nd = length(d)
      dod = d * rep(d, rep(nd, nd))
      dim(dod) = c(nd,nd)
      dog = d * rep(g, rep(nd, nd))
      dim(dog) = c(nd,nd)
      dogD = dog %*% D
      D = D + (1 + drop(t(g) %*% D %*% g) / dg) * dod / dg -
        (dogD + t(dogD)) / dg
    }
  }
  r$num.iterations = i
  # print(-diag(D))
  r
}

# Line search for beta and mix

# If brkt is TRUE, (mix1, beta1) must be an interior point

# Two conditions must be satisfied at termination: (a) the Armijo rule; (b) 

# mix1       Current disc object
# beta1      Current structural parameters
# mix2       Next disc object, of the same size as mix1
# beta2      Next structural parameters
# x          Data
# maxit      Maximum number of iterations
# which      Indicators for pi, theta, beta, which their first derivatives
#            of the log-likelihood should be computed and examined
# brkt       Bracketing first?

# Output:

# mix          Estimated mixing distribution
# beta         Estimated beta
# ll           Log-likelihood value at the estimate
# convergence  = 0, converged successfully; = 1, failed

lsch = function(mix1, beta1, mix2, beta2, x, 
          maxit=100, which=c(1,1,1), brkt=FALSE ) {
  k = length(mix1$pt)
  je = (mix1$pt == mix2$pt) | (mix1$pt == mix2$pt)
  convergence = 1
  dl1 = dll.snpmle(x, mix1, beta1, which=c(1,which))
  lla = ll1 = dl1$ll
  names.grad = c(if(which[1]) paste0("pr",1:k) else NULL,
    if(which[2]) paste0("pt",1:k) else NULL,
    if(which[3]) paste0("beta",1:length(beta1)) else NULL)
  grad1 = c(if(which[1]) dl1$dp else NULL,
    if(which[2]) dl1$dt else NULL,
    if(which[3]) dl1$db else NULL)
  names(grad1) = names.grad
  d1 = c(if(which[1]) mix2$pr - mix1$pr else NULL,
      if(which[2]) mix2$pt - mix1$pt else NULL,
      if(which[3]) beta2 - beta1 else NULL)
  d1.norm = sqrt(sum(d1*d1))
  s = d1 / d1.norm
  g1d1 = sum(grad1 * d1)
  dla = g1s = g1d1 / d1.norm
  if(d1.norm == 0 || g1s <= 0) {
    # cat("Warning in lsch(): d1.norm =", d1.norm, " g1s =", g1s, "\n")
    return( list(mix=mix1, beta=beta1, grad=grad1, ll=ll1, convergence=3) )
  }
  a = 0
  b = 1
  if(which[1] && any(mix2$pr == 0)) brkt = FALSE
  for(i in 1:maxit) {
    for(j in 1:1000) {
      m = disc( (1-b) * mix1$pt + b * mix2$pt, (1-b) * mix1$pr + b * mix2$pr)
      m$pt[je] = mix1$pt[je]
      # beta = if(is.null(beta1)) NULL else (1-b) * beta1 + b * beta2
      beta = if(which[3]) (1-b) * beta1 + b * beta2 else beta1
      if( valid.snpmle(x, beta, m) ) break
      brkt = FALSE
      b = 0.5 * a + 0.5 * b
    }
    if(j == 1000) stop("Can not produce valid interior point in lsch()")
    dl = dll.snpmle(x, m, beta, which=c(1,which))
    ll = dl$ll
    grad = c(if(which[1]) dl$dp else NULL,
      if(which[2]) dl$dt else NULL,
      if(which[3]) dl$db else NULL)
    gs = sum(grad * s)
    # With the following condition, "ll(a)" stays above the Armijo line and
    # "ll(b)" either fails to meet the gradient condition or falls below
    # the Armijo line.  Therefore, bracketing must have succeeded
    if( brkt && gs > g1s * .5 && ll >= ll1 + g1d1 * b * .33 )
      {a = b; b = 2 * b; lla = ll; dla = gs}
    else break
  }
  if(i == maxit) brkt = FALSE
  alpha = b; llb = ll; dlb = gs
  for(i in 1:maxit) {
    g1d = g1d1 * alpha
    if(is.na(ll)) ll = -Inf
    if( ll >= ll1 - 1e-15 * abs(ll1) && g1d <= 1e-15 * abs(ll))   # Flat land
      { # cat("Warning in lsch(): flat land reached.\n");
        convergence=2; break} 
    if( brkt ) {          # Armijo rule + slope test
      if( ll >= ll1 + g1d * .33 && abs(gs) <= g1s * .5)
        {convergence=0; break}
      if( ll >= ll1 + g1d * .33 && gs > g1s * .5 )
        { a = alpha; lla = ll; dla = gs }
      else { b = alpha; llb = ll; dlb = gs }
    }
    else {                  # Armijo rule
      if( ll >= ll1 + g1d * .33 ) {convergence=0; break}
      else { b = alpha; llb = ll; dlb = gs }
    }
    alpha = (a + b) * .5
    m = disc( (1-alpha) * mix1$pt + alpha * mix2$pt,
      (1-alpha) * mix1$pr + alpha * mix2$pr)
    m$pt[je] = mix1$pt[je]
    beta = if(which[3]) (1-alpha) * beta1 + alpha * beta2 else beta1
#     beta = if(is.null(beta1)) NULL else (1-alpha) * beta1 + alpha * beta2
    dl = dll.snpmle(x, m, beta, which=c(1,which))
    ll = dl$ll
    grad = c(if(which[1]) dl$dp else NULL,
      if(which[2]) dl$dt else NULL,
      if(which[3]) dl$db else NULL)
    gs = sum(grad * s)
  }
  # if(i == 100) { print("i = 100 in lsch()"); stop() }   # should never happen
  names(grad) = names.grad
  beta = if(which[3]) (1-alpha) * beta1 + alpha * beta2 else beta1
  mix = disc((1-alpha)*mix1$pt+alpha*mix2$pt, (1-alpha)*mix1$pr+alpha*mix2$pr)
  mix$pt[je] = mix1$pt[je]
  list(mix=mix, beta=beta, grad=grad,
       ll=ll, convergence=convergence, num.iterations=i)
}

# The Wolfe-Powell conditions are satisfied after line search.

# An exact line search is perhaps unnecessary, since BFGS has superlinear
# convergence

# Input:
#
# r       A list containing "mix" and "beta", for the current estimates
#         of th mixing distribution and beta, respectively 
# beta2   Tentative new beta
# x       Data
# tol     Tolerance level 
# ...     Arguments passed to pll

# Output:
#
# outputs from function pll, plus
# convergence     0, converged successfully;
#                 1, failed (likely due to precision limit reached)
# num.npmle       0, number of the NPMLEs computed

lsch.pll = function(r, beta2, x, tol=1e-15, maxit=100, tol.npmle=tol*1e-4,
          brkt=FALSE, ...) {
  convergence = 0
  num.npmle = 0
  lla = ll1 = r$ll
  d1 = beta2 - r$beta
  d1.norm = sqrt(sum(d1*d1))
  s = d1 / d1.norm
  dla = g1s = sum(r$grad * s)
  g1d1 = sum(r$grad * d1)
  # if(d1.norm == 0 || g1s <= 0) return(r)
  a = 0
  alpha = 1
  for(i in 1:maxit) {
    repeat {
      beta = (1-alpha) * r$beta + alpha * beta2
      if( valid.snpmle(x, beta, r$mix) ) break
      alpha = a + 0.9 * (alpha - a)
      brkt = FALSE
    }
    r.alpha = pll(beta, r$mix, x, tol=tol.npmle, ...)
    num.npmle = num.npmle + 1
    d = beta - r$beta
    gs = sum(r.alpha$grad * s)
    if( brkt && gs > g1s * .5 && r.alpha$ll >= ll1 + g1d1 * alpha * .33)
      {a = alpha; alpha = 2 * alpha; lla = r.alpha$ll; dla = gs}
    else break
    # if( sum(r.alpha$grad * r.alpha$grad) < 1e-12 ) {brkt=FALSE; break}
  }
  if(i == maxit) brkt = FALSE
  b = alpha; llb = r.alpha$ll; dlb = gs
  for(i in 1:maxit) {
    if(b - a <= 1e-10) {convergence=1; break}
#    if(r.alpha$ll >= r$ll - tol * abs(r.alpha$ll) &&
#       r.alpha$ll <= r$ll  + tol * abs(r.alpha$ll)) break
    if(r.alpha$ll >= r$ll - tol && r.alpha$ll <= r$ll  + tol) break
    g1d = g1d1 * alpha
    if( brkt ) {
      if( r.alpha$ll >= r$ll + g1d * .33 & abs(gs) <= g1s * .5 ) break
      if( r.alpha$ll >= r$ll + g1d * .33 & gs > g1s  * .5 )
        {a = alpha; lla = r.alpha$ll; dla = gs}
      else {b = alpha; llb = r.alpha$ll; dlb = gs}
    }
    else {
      if( r.alpha$ll >= r$ll + g1d * .33 ) {convergence=0; break}
      else {b = alpha; llb = r.alpha$ll; dlb = gs}
    }
    alpha = (a + b) * .5
    beta = (1-alpha) * r$beta + alpha * beta2
    r.alpha = pll(beta, r$mix, x, tol=tol.npmle, ...)
    num.npmle = num.npmle + 1
    d = beta - r$beta
    gs = sum(r.alpha$grad * s)
  }
  r.alpha$convergence = convergence
  r.alpha$num.npmle = num.npmle
  r.alpha
}

logLik.snpmle = function(x, beta, mix) {
  ld = logd(x, beta, mix$pt, which=c(1,0,0))$ld
  ma = matMaxs(ld)
  dmix = drop(exp(ld - ma) %*% mix$pr) + 1e-100
  logd = log(dmix) + ma
  ll = llex(x, beta, mix) + sum( weights(x, beta) * logd )
  attr(ll, "dmix") = dmix
  attr(ll, "logd") = logd    # log(mixture density)
  ll
}

density.snpmle = function(x, beta, mix) {
  rowSums(exp(logd(x, beta, mix$pt, which=c(1,0,0))$ld +
      rep(log(mix$pr), rep(length(x), length(mix$pr)))))
}

# Profile log-likelihood

pll = function(beta, mix0, x, tol=1e-15, grid=100, ...) {
  r = cnm(x, init=list(beta=beta, mix=mix0), tol=tol, grid=grid, maxit=50, ...)
  mix = r$mix
  dl = logd(x, beta, mix$pt, which=c(1,1,0))
  lpt = dl$ld
  ma = matMaxs(lpt)
  dmix = drop(exp(lpt - ma) %*% mix$pr) + 1e-100
  D = pmin(exp(lpt - ma), 1e100)
  p = weights(x, beta) * D/dmix * rep(mix$pr, rep(nrow(D), length(mix$pr)))
  db = apply(sweep(dl$db, c(1,2), p, "*"), c(1,3), sum)
  g = llexdb(x,beta,mix) + colSums(db)
  list(ll=r$ll, beta=r$beta, mix=mix, grad=g, max.gradient=r$max.gradient,
       num.iterations=r$num.iter)
}

## Compute all local maxima of the gradient function using a hybrid
## secant-bisection method. It does not need the second derivative of the
## gradient function.

## O(n*m) + #iteration * O(n*mhat)
## n - sample size
## m - grid points
## mhat - number of local maxima

maxgrad = function(x, beta, mix, grid=100, tol=-Inf, maxit=100) {
  if( length(grid) == 1 ) grid = gridpoints(x, beta, grid)
  dg = grad(grid, x, beta, mix, order=1)$d1
  j = is.na(dg)
  grid = grid[!j]
  tol.pt = 1e-14 * diff(range(grid))
  dg = dg[!j]
  np = length(grid)
  jmax = which(dg[-np] > 0 & dg[-1] < 0)
  if( length(jmax) < 1 ) return
  left = grid[jmax]
  right = grid[jmax+1]
  pt = (left + right) * .5
  if(length(pt) != 0) {
    pt.old = left
    d1.old = dg[jmax]
    d2 = rep(-1, length(pt))
    for(i in 1:maxit) {
      d1 = grad(pt, x, beta, mix, order=1)$d1
      d2t = (d1 - d1.old) / (pt - pt.old)
      jd = !is.na(d2t) & d2t < 0
      d2[jd] = d2t[jd]
      left[d1>0] = pt[d1>0]
      right[d1<0] = pt[d1<0]
      pt.old = pt
      d1.old = d1
      pt = pt - d1 / d2
      j = is.na(pt) | pt < left | pt > right     # those out of brackets
      pt[j] = (left[j] + right[j]) * .5
      if( max(abs(pt - pt.old)) <= tol.pt) break
    }
  }
  else i = 0
  if(dg[np] >= 0) pt = c(pt, grid[np])
  if(dg[1] <= 0) pt = c(grid[1], pt)
  if(length(pt) == 0) stop("no new support point found") # should never happen
  g = grad(pt, x, beta, mix, order=0)$d0
  names(pt) = names(g) = NULL
  j = g >= tol
  list(pt=pt[j], grad=g[j], num.iterations=i)
}

# gradient function, with O(n*m)

grad = function(pt, x, beta, mix, order=0) {
  w = weights(x, beta)
  if(length(w) != length(x)) w = rep(w, length=length(x))
  if(is.null(attr(mix, "ll")) || is.null(attr(attr(mix, "ll"), "dmix")))
    attr(mix, "ll") = logLik.snpmle(x, beta, mix)
  ## dmix = attr(attr(mix, "ll"), "dmix")
  ## ma = attr(attr(mix, "ll"), "ma")
  logd.mix = attr(attr(mix, "ll"), "logd")
  g = vector("list", length(order))
  names(g) = paste0("d", order)
  which = c(1,0,0)
  if(any(order >= 1)) which[3:max(order+2)] = 1
  dl = logd(x, beta, pt, which=which)
  ws = pmin(exp(dl$ld + (log(w) - logd.mix)), 1e100)
  if(0 %in% order) g$d0 = colSums(ws) - sum(w)
  if(1 %in% order) g$d1 = colSums(ws * dl$dt)
  g
}

plotgrad = function(x, mix, beta, len=500, order=0, col="blue",
                    col2="red", add=FALSE,
                    main=paste0("Class: ",class(x)),
                    xlab=expression(theta),
                    ylab=paste0("Gradient (order = ",order,")"), cex=1,
                    pch=1, lwd=1, xlim, ylim, ...) {
  if(order > 1) stop("order > 1")
  if( missing(xlim) ) xlim = range(gridpoints(x, beta, grid=2))
  pt = seq(xlim[1], xlim[2], len=len)
  if(class(mix) == "disc") pt = sort(c(mix$pt, pt))
  g = grad(pt, x, beta, mix, order=order)[[1]]
  if(missing(ylim)) ylim = range(0,g,na.rm=TRUE)
  if(!add) plot(xlim, c(0,0), type="l", col="black", ylim=ylim, main=main,
       xlab=xlab, ylab=ylab, cex = cex, cex.axis = cex, cex.lab = cex, ...)
  lines(pt, g, col=col, lwd=lwd)
  if(class(mix) == "disc") {
    j = mix$pr > 0
    points(mix$pt[j], rep(0,length(mix$pt[j])), pch=pch, col=col2)
    lines(mix$pt[j], mix$pr[j]*min(g), type="h", col=col2)
  }
  invisible(list(pt=pt,g=g))
}

# The functions examines whether the initial values are proper. If not,
# proper ones are provided, by employing the user's function "initial."

initial.snpmle = function(x, init=NULL, kmax=NULL) {
  if(!is.null(kmax) && kmax == Inf) kmax = NULL
  init = initial(x, init$beta, init$mix, kmax=kmax)
  init$ll = logLik.snpmle(x, init$beta, init$mix)
    # ll = logLik.snpmle(x, init0$beta, init$mix)
#     if(ll > init0$ll) init = list(beta=init0$beta, mix=init$mix, ll=ll)
#    else init = init0
  if(! valid.snpmle(x, init$beta, init$mix)) stop("Invalid initial values!")
  init
}

# Compute the first derivatives of the log-likelihood function of pi,
# theta and beta

# INPUT:
# 
# which     A 4-vector having values either 0 and 1, indicating which
#           derivatives to be computed, in the order of pi, theta and beta
# ind       = TRUE, return the derivatives with respect to individual
#                   observations (with weights)
#           = FALSE, return the derivatives
#
# OUTPUT:   a list with ll, dp, dt, db, as specified by which

dll.snpmle = function(x, mix, beta, which=c(1,0,0,0), ind=FALSE) {
  w = weights(x, beta)
  r = list()
  dl = logd(x, beta, mix$pt, which=c(1,which[4:3]))
  lpt = dl$ld
  ma = matMaxs(lpt)
  if(which[1] == 1) {
    r$ll = log(drop(exp(lpt - ma) %*% mix$pr)) + ma
    if(!ind) r$ll = llex(x, beta, mix) + sum(w * r$ll)
  }
  if(sum(which[2:4]) == 0) return(r)
  dmix = drop(exp(lpt - ma) %*% mix$pr) + 1e-100
  S = pmin(exp(lpt - (ma + log(dmix))), 1e100)
  # dp = D / dmix
  if(which[2] == 1) {
    if(ind) r$dp = S
    else r$dp = colSums(w * S)
  }
  if(sum(which[3:4]) == 0) return(r)
  p = S * rep(mix$pr, rep(nrow(S), length(mix$pr)))
  if(which[3] == 1) {
    r$dt = p * dl$dt
    if(!ind) r$dt = colSums(w * r$dt)
  }
  if(which[4] == 0) return(r)
  dl1 = dl$db
  if(is.null(dl1)) r$db = NULL
  else {
    r$db = apply(sweep(dl1, c(1,2), p, "*"), c(1,3), sum)
    if(!ind) r$db = llexdb(x, beta, mix) + colSums(w * r$db)
  }
  r
}

valid.snpmle = function(x, beta, mix) {
  bs = suppspace(x, beta)
  bs = bs + 1e-14 * c(-1,1) * abs(bs)             # tolerance
  valid(x, beta, mix) && all(mix$pr >= 0, mix$pt >= bs[1], mix$pt <= bs[2])
}

collapse.snpmle = function(mix, beta, x, tol=1e-15) {
  j = mix$pr == 0
  if(any(j)) {
    mix = disc(mix$pt[!j], mix$pr[!j])
    attr(mix, "ll") = logLik.snpmle(x, beta, mix)
  }
  else { if( is.null(attr(mix, "ll")) || is.null(attr(attr(mix, "ll"), "dmix")) )
           attr(mix, "ll") = logLik.snpmle(x, beta, mix) }
##   if( any(mix$pr < tol.p) ) {
##     j = mix$pr >= tol.p
##     mixt = mix
##     mixt$pt = mixt$pt[j]
##     mixt$pr = mixt$pr[j] / sum(mix$pr[j])
##     attr(mixt, "ll") = logLik.snpmle(x, beta, mixt)
##     if(attr(mixt,"ll")[1] >= attr(mix,"ll")[1] - tol) mix = mixt
##   }
  
##   p = attr(attr(mix, "ll"), "p")
##   pj = colSums(weights(x, beta) * p)
##   if( any(pj == 0) ) {
##     j = pj != 0
##     mix$pt = mix$pt[j]
##     mix$pr = mix$pr[j] / sum(mix$pr[j])
##     attr(mix, "ll") = logLik.snpmle(x, beta, mix)
##   }
  
  if(tol > 0) {
    repeat {
      if(length(mix$pt) <= 1) break
      prec = 10 * min(diff(mix$pt))   ## one digit at a time
      ## if(prec > diff(range(gridpoints(x, beta, grid=2))) * .01) break
      mixt = unique(mix, prec=c(prec,0))
      attr(mixt, "ll") = logLik.snpmle(x, beta, mixt)
      # if(attr(mixt,"ll")[1] >= attr(mix,"ll")[1] - tol * abs(attr(mix,"ll")[1]))
      if(attr(mixt,"ll")[1] >= attr(mix,"ll")[1] - tol) mix = mixt  
      else break
    }
  }
  mix
}

print.snpmle = function(verbose, ll, mix, beta, maxdg, i) {
  if(verbose > 0)
    cat("## Iteration ", i, ":\n   Log-likelihood = ", as.character(ll), "\n", sep="")
  if(verbose > 1) cat("   Max gradient =", maxdg, "\n")
  if(verbose > 2) cat("   theta =", signif(mix$pt,6), "\n")
  if(verbose > 3) cat("   Proportions =", signif(mix$pr,6), "\n")
  if(verbose > 4)
    cat("   beta =", if(is.null(beta)) "NULL" else signif(beta,6), "\n")
}

plot.nspmix = function(x, data, type=c("probability","gradient"), ...) {
  type = match.arg(type)
  if(missing(data)) data = NULL
  else {
    if(is.null(class(data))) stop("Data belongs to no class")
    if(class(data)[1] != x$family)
      stop(paste0("Class of data does not match the mixture family: ",
                 class(data)[1], " != ", x$family))
  }
  plotf = getFunction(paste0("plot.",x$family))
  switch(type,
         "probability" = plotf(data, x$mix, x$beta, ...),
         "gradient" = plotgrad(data, x$mix, x$beta, pch=1, ...) )
  invisible(x)
}

### old functions

# gradient function

OLD.grad = function(pt, x, beta, dmix, ma, order=0) {
  w = weights(x, beta)
  if(length(w) != length(x)) w = rep(w, length=length(x))
  if(class(dmix) == "disc") {
    l = logd(x, beta, dmix$pt, which=c(1,0,0))$ld
    ma = matMaxs(l)
    dmix = drop(exp(l - ma) %*% dmix$pr) + 1e-100
  }
  g = vector("list", length(order))
  names(g) = paste0("d", order)
  which = c(1,0,0)
  if(any(order >= 1)) which[3:max(order+2)] = 1
  dl = logd(x, beta, pt, which=which)
  ws = w * pmin(exp(dl$ld - ma), 1e100) / dmix
  if(0 %in% order) g$d0 = colSums(ws) - sum(w)
  if(1 %in% order) g$d1 = colSums(ws * dl$dt)
  g
}

NEW.maxgrad = function(x, beta, mix, grid=100, tol=-Inf, maxit=100) {
  if( length(grid) == 1 ) grid = gridpoints(x, beta, grid)
  dg = grad(grid, x, beta, mix, order=1)$d1
  j = is.na(dg)
  grid = grid[!j]
  tol.pt = 1e-14 * diff(range(grid))
  dg = dg[!j]
  np = length(grid)
  jmax = which(dg[-np] > 0 & dg[-1] < 0)
  if( length(jmax) < 1 ) return
  left = grid[jmax]
  right = grid[jmax+1]
  pta = pt = (left + right) * .5
  if(length(pt) != 0) {
    pt.old = left
    d1.old = dg[jmax]
    d2 = rep(-1, length(pt))
    ju = 1:length(pt)            # indexes for unsolved points
    for(i in 1:maxit) {
      d1 = grad(pt, x, beta, mix, order=1)$d1
      d2t = (d1 - d1.old) / (pt - pt.old)
      jd = !is.na(d2t) & d2t < 0
      d2[jd] = d2t[jd]
      left[d1>0] = pt[d1>0]
      right[d1<0] = pt[d1<0]
      pt.old = pt
      d1.old = d1
      pt = pt - d1 / d2
      j = is.na(pt) | pt < left | pt > right     # those out of brackets
      pt[j] = (left[j] + right[j]) * .5
      jj = abs(pt - pt.old) > tol.pt             # unsolved cases
      sjj = sum(jj)
      if(any(!jj)) pta[ju[!jj]] = pt[!jj]        # cases solved
      if(sjj == 0) break
      if(sjj < length(pt)) {
        pt = pt[jj]
        pt.old = pt.old[jj]
        d1.old = d1.old[jj]
        d2 = d2[jj]
        left = left[jj]
        right = right[jj]
        ju = ju[jj]
      }
    }
  }
  else i = 0
  if(dg[np] >= 0) pta = c(grid[np], pta)
  if(dg[1] <= 0) pta = c(grid[1], pta)
  if(length(pt) == 0) stop("no new support point found") # should never happen
  g = grad(pta, x, beta, mix, order=0)$d0
  names(pta) = names(g) = NULL
  j = g >= tol
  list(pt=pta[j], grad=g[j], num.iterations=i)
}

