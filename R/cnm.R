# ==================================================== #
# Maximum likelihood computation for nonparametric and #
# semiparametric mixture model                         # 
# ==================================================== #

##' @title Valid parameter values
##'
##' @description A generic method used to return \code{TRUE} if the
##'   values of the paramters use for a nonparametric/semiparametric
##'   mixture are valid, or \code{FALSE} if otherwise.
##'
##' @param x an object of a class for data.
##' @param beta instrumental parameter in a semiparametric mixture.
##' @param theta values of the mixing variable.
##' @author Yong Wang <yongwang@@auckland.ac.nz>
##' 
##' @return
##' 
##'   A logical value.
##' 
##' @export

valid = function(x, beta, theta) UseMethod("valid")

##' @title Support space
##'
##' @description Range of the mixing variable (theta). 
##'
##' @param x an object of a class for data.
##' @param beta instrumental parameter in a semiparametric mixture.
##' @author Yong Wang <yongwang@@auckland.ac.nz>
##'
##' @return
##' 
##'  A vector of length 2.
##'
##' @export

suppspace = function(x, beta) UseMethod("suppspace")

##' @title Grid points
##'
##' @description A generic method used to return a vector of grid
##'   points used for searching local maxima of the gradient funcion.
##'
##' @param x an object of a class for data.
##' @param beta instrumental parameter in a semiparametric mixture.
##' @param grid number of grid points to be generated.
##' 
##' @author Yong Wang <yongwang@@auckland.ac.nz>
##' 
##' @return
##' 
##'   A numeric vector containing grid points.
##'
##' @export

gridpoints = function(x, beta, grid) UseMethod("gridpoints")

##' @title Initialization for a nonparametric/semiparametric mixture
##'
##' @description A generic method used to return an initialization for
##'   a nonparametric/semiparametric mixture.
##'
##' @param x an object of a class for data.
##' @param beta instrumental parameter in a semiparametric mixture.
##' @param mix an object of class \code{disc} for the mixing
##'   distribution.
##' @param kmax the maximum allowed number of support points used.
##' 
##' @author Yong Wang <yongwang@@auckland.ac.nz>
##' 
##' @return
##'
##' \item{beta}{an initialized value of beta}
##' 
##' \item{mix}{an initialised or updated object of class \code{disc}
##' for the mixing distribution.}
##'
##' @export

initial = function(x, beta, mix, kmax) UseMethod("initial")

##' @title Log-density and its derivative values
##'
##' @description A generic method to compute the log-density values
##'   and possibly their first derivatives with respec to theta and
##'   beta.
##'
##' @param x an object of a class for data.
##' @param beta instrumental parameter in a semiparametric mixture.
##' @param pt a vector of values for the mixing variable theta.
##' @param which an integer vector of length 3, indicating if,
##'   respectively, the log-density values, the derivatives wrt beta
##'   and the derivatives wrt theta are to be computed and returned if
##'   being 1 (\code{TRUE}).
##' 
##' @author Yong Wang <yongwang@@auckland.ac.nz>
##' 
##' @return 
##'
##' \item{ld}{a matrix, storing the log-density values for each (x[i],
##' beta, pt[j], or NULL if not asked for.}
##' 
##' \item{db}{a matrix, storing the log-density derivatives wrt beta
##' for each (x[i], beta, pt[j], or NULL if not asked for.}
##' 
##' \item{dt}{a matrix, storing the log-density derivatives wrt theta
##' for each (x[i], beta, pt[j], or NULL if not asked for.}
##'
##' @export

logd = function(x, beta, pt, which) UseMethod("logd")

##' @title Weights
##'
##' @description Weights or frequencis of observations.
##'
##' @param x an object of a class for data.
##' @param beta instrumental parameter in a semiparametric mixture.
##' 
##' @author Yong Wang <yongwang@@auckland.ac.nz>
##' 
##' @return 
##'
##' a numeric vector of the weights.
##' 
##' @export

weight = function(x, beta) UseMethod("weight")

##' @title Log-likleihood Extra Term. 
##'
##' @description Value of possibly an extra term in the log-likleihood
##'   function for the instrumental parameter beta
##'
##' @param x an object of a class for data.
##' @param beta instrumental parameter in a semiparametric mixture.
##' @param mix an object of class \code{disc} for the mixing
##'   distribution.
##' 
##' @author Yong Wang <yongwang@@auckland.ac.nz>
##' 
##' @return 
##'
##' a scalar value
##' 
##' @export

llex = function(x, beta, mix) UseMethod("llex")

##' @title Derivative of the log-likleihood Extra Term
##'
##' @description Derivative of the log-likleihood extra term wrt beta
##'
##' @param x an object of a class for data.
##' @param beta instrumental parameter in a semiparametric mixture.
##' @param mix an object of class \code{disc} for the mixing
##'   distribution.
##' 
##' @author Yong Wang <yongwang@@auckland.ac.nz>
##' 
##' @return 
##'
##' a scalar value
##' 
##' @export

llexdb = function(x, beta, mix) UseMethod("llexdb")

##' @export

valid.default = function(x, beta, theta) TRUE

##' @export

suppspace.default = function(x, beta) c(-Inf, Inf)

##' @export

llex.default = function(x, beta, mix) 0

##' @export

llexdb.default = function(x, beta, mix) 0








##' @title Maximum Likelihood Estimation of a Nonparametric Mixture Model
##'
##' @description Function \code{cnm} can be used to compute the maximum
##'   likelihood estimate of a nonparametric mixing distribution
##'   (NPMLE) that has a one-dimensional mixing parameter, or simply
##'   the mixing proportions with support points held fixed.
##'
##'
##' A finite mixture model has a density of the form
##'
##' \deqn{f(x; \pi, \theta, \beta) = \sum_{j=1}^k \pi_j f(x; \theta_j,
##' \beta).}{f(x; pi, theta, beta) = sum_{j=1}^k pi_j f(x; theta_j,
##' beta),}
##'
##' where \eqn{\pi_j \ge 0}{pi_j >= 0} and \eqn{\sum_{j=1}^k \pi_j
##' }{sum_{j=1}^k pi_j =1}\eqn{ =1}{sum_{j=1}^k pi_j =1}.
##'
##' A nonparametric mixture model has a density of the form
##'
##' \deqn{f(x; G) = \int f(x; \theta) d G(\theta),}{f(x; G) = Integral
##' f(x; theta) d G(theta),} where \eqn{G} is a mixing distribution
##' that is completely unspecified. The maximum likelihood estimate of
##' the nonparametric \eqn{G}, or the NPMLE of \eqn{G}, is known to be
##' a discrete distribution function.
##'
##' Function \code{cnm} implements the CNM algorithm that is proposed
##' in Wang (2007) and the hierarchical CNM algorithm of Wang and
##' Taylor (2013). The implementation is generic using S3
##' object-oriented programming, in the sense that it works for an
##' arbitrary family of mixture models defined by the user.  The user,
##' however, needs to supply the implementations of the following
##' functions for their self-defined family of mixture models, as they
##' are needed internally by function \code{cnm}:
##'
##' \code{initial(x, beta, mix, kmax)}
##'
##' \code{valid(x, beta, theta)}
##'
##' \code{logd(x, beta, pt, which)}
##'
##' \code{gridpoints(x, beta, grid)}
##'
##' \code{suppspace(x, beta)}
##'
##' \code{length(x)}
##'
##' \code{print(x, ...)}
##'
##' \code{weight(x, ...)}
##'
##' While not needed by the algorithm for finding the solution, one
##' may also implement
##'
##' \code{plot(x, mix, beta, ...)}
##'
##' so that the fitted model can be shown graphically in a
##' user-defined way.  Inside \code{cnm}, it is used when
##' \code{plot="probability"} so that the convergence of the algorithm
##' can be graphically monitored.
##'
##' For creating a new class, the user may consult the implementations
##' of these functions for the families of mixture models included in
##' the package, e.g., \code{npnorm} and \code{nppois}.
##'
##' @param x a data object of some class that is fully defined by the
##'   user. The user needs to supply certain functions as described
##'   below.
##' @param init list of user-provided initial values for the mixing
##'   distribution \code{mix} and the structural parameter
##'   \code{beta}.
##' @param model the type of model that is to estimated: the
##'   non-parametric MLE (if \code{npmle}), or mixing proportions only
##'   (if \code{proportions}).
##' @param maxit maximum number of iterations.
##' @param tol a tolerance value needed to terminate an
##'   algorithm. Specifically, the algorithm is terminated, if the
##'   increase of the log-likelihood value after an iteration is less
##'   than \code{tol}.
##' @param grid number of grid points that are used by the algorithm to
##'   locate all the local maxima of the gradient function. A larger
##'   number increases the chance of locating all local maxima, at the
##'   expense of an increased computational cost. The locations of the
##'   grid points are determined by the function \code{gridpoints}
##'   provided by each individual mixture family, and they do not have
##'   to be equally spaced. If needed, a \code{gridpoints} function
##'   may choose to return a different number of grid points than
##'   specified by \code{grid}.
##' @param plot whether a plot is produced at each iteration. Useful
##'   for monitoring the convergence of the algorithm. If
##'   \code{="null"}, no plot is produced. If \code{="gradient"}, it
##'   plots the gradient curves and if \code{="probability"}, the
##'   \code{plot} function defined by the user for the class is used.
##' @param verbose verbosity level for printing intermediate results in
##'   each iteration, including none (= 0), the log-likelihood value
##'   (= 1), the maximum gradient (= 2), the support points of the
##'   mixing distribution (= 3), the mixing proportions (= 4), and if
##'   available, the value of the structural parameter beta (= 5).
##' 
##' @return
##'
##' \item{family}{the name of the mixture family that is used to fit
##' to the data.}
##'
##' \item{num.iterations}{number of iterations required by the
##' algorithm}
##'
##' \item{max.gradient}{maximum value of the gradient function,
##' evaluated at the beginning of the final iteration}
##'
##' \item{convergence}{convergence code. \code{=0} means a success,
##' and \code{=1} reaching the maximum number of iterations}
##'
##' \item{ll}{log-likelihood value at convergence}
##'
##' \item{mix}{MLE of the mixing distribution, being an object of the
##' class \code{disc} for discrete distributions.}
##'
##' \item{beta}{value of the structural parameter, that is held fixed
##' throughout the computation.}
##' 
##' @author Yong Wang <yongwang@@auckland.ac.nz>
##' 
##' @seealso \code{\link[lsei]{nnls}}, \code{\link{npnorm}},
##'   \code{\link{nppois}}, \code{\link{cnmms}}.
##' 
##' @references
##'
##' Wang, Y. (2007). On fast computation of the non-parametric maximum
##' likelihood estimate of a mixing distribution. \emph{Journal of the
##' Royal Statistical Society, Ser. B}, \bold{69}, 185-198.
##'
##' Wang, Y. (2010). Maximum likelihood computation for fitting
##' semiparametric mixture models. \emph{Statistics and Computing},
##' \bold{20}, 75-86
##'
##' Wang, Y. and Taylor, S. M. (2013). Efficient computation of
##' nonparametric survival functions via a hierarchical mixture
##' formulation. \emph{Statistics and Computing}, \bold{23}, 713-725.
##' 
##' @examples
##'
##' ## Simulated data
##' x = rnppois(200, disc(c(1,4), c(0.7,0.3))) # Poisson mixture
##' (r = cnm(x))
##' plot(r, x)
##'
##' x = rnpnorm(200, disc(c(0,4), c(0.3,0.7)), sd=1) # Normal mixture
##' plot(cnm(x), x)                        # sd = 1
##' plot(cnm(x, init=list(beta=0.5)), x)   # sd = 0.5
##'
##' ## Real-world data
##' data(thai)
##' plot(cnm(x <- nppois(thai)), x)     # Poisson mixture
##'
##' data(brca)
##' plot(cnm(x <- npnorm(brca)), x)     # Normal mixture
##'
##'
##' @importFrom graphics barplot hist lines plot points rect segments abline
##' @importFrom graphics matplot
##' @importFrom methods getFunction
##' @importFrom stats aggregate binomial glm quantile rmultinom var
##' @importFrom stats weighted.mean
##' @import lsei
##' 
##' @export

cnm = function(x, init=NULL, model=c("npmle","proportions"), maxit=100,
               tol=1e-6,
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
  init = initial0(x, init)
  beta = init$beta
  nb = length(beta)
  mix = init$mix
  attr(mix, "ll") = init$ll
  convergence = 1
  gp = gridpoints(x, beta, grid)          # does not depend on beta
  w = weight(x, beta)
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
      if(plot=="gradient") points(g$pt, g$grad, col=2, pch=20)
      gmax = g$gmax
      mix = disc(c(mix$pt,g$pt), c(mix$pr,rep(0,length(g$pt))))
    }
    else {
      gpt0 = grad(pt0, x, beta, mix)$d0
      ind = indx(pt0, mix$pt)
      ind[ind == 0] = 1
      s = split(gpt0, ind)
      imax = sapply(s, which.max)
      gmax = max(gpt0[imax])
      imax = imax[imax > 1]
      csum = c(0,cumsum(sapply(s, length)))
      ipt = csum[as.numeric(names(imax))] + imax
      mix = disc(c(mix$pt, pt0[ipt]), c(mix$pr, double(length(ipt))))
    }
    lpt = logd(x, beta, mix$pt, which=c(1,0,0))$ld
    D = pmin(exp(lpt - ma), 1e100)
    if(verbose > 0) 
      verb(verbose, attr(mix1, "ll")[1], mix1, beta, gmax, i-1)
    r = hcnm(D, mix$pr, w=w, maxit=3, compact=FALSE)
    mix$pr = r$p
    # attr(mix, "ll") = r$ll + sum(w * ma)
    mix = collapse0(mix, beta, x, tol=0)
    if( attr(mix, "ll")[1] <= attr(mix1, "ll")[1] + tol )
      {convergence = 0; break}
  }
  if(model == "npmle")
    mix = collapse0(mix, beta, x,
                          tol=max(tol*.1, abs(attr(mix, "ll")[1])*1e-16))
  ll = attr(mix,"ll")[1]
  attr(mix,"ll") = NULL
  structure(list(family=class(x)[1], num.iterations=i, max.gradient=gmax,
                 convergence=convergence, ll=ll, mix=mix, beta=beta),
            class="nspmix")
}

##' @title Hierarchical Constrained Newton method
##' 
##' @description Function \code{hcnm} can be used to compute the MLE
##'   of a finite discrete mixing distribution, given the component
##'   density values of each observation. It implements the
##'   hierarchical CNM algorithm of Wang and Taylor (2013).
##' 
##' @param D A numeric matrix, each row of which stores the component
##'   density values of an observation.
##' @param p0 Initial mixture component proportions.
##' @param w Duplicity of each row in matrix \code{D} (i.e., that of a
##'   corresponding observation).
##' @param maxit Maximum number of iterations.
##' @param tol A tolerance value to terminate the
##'   algorithm. Specifically, the algorithm is terminated, if the
##'   increase of the log-likelihood value after an iteration is less
##'   than \code{tol}.
##' @param blockpar Block partitioning parameter. If > 1, the number
##'   of blocks is roughly \code{nrol(D)/blockpar}. If < 1, the number
##'   of blocks is roughly \code{nrol(D)^blockpar}.
##' @param recurs.maxit Maximum number of iterations in recursions.
##' @param compact Whether iteratively select and use a compact subset
##'   (which guarantees convergence), or not (if already done so
##'   before calling the function).
##' @param depth Depth of recursion/hierarchy.
##' @param verbose Verbosity level for printing intermediate results.
##'  
##' @return
##' 
##' \item{p}{Computed probability vector.}
##' 
##' \item{convergence}{convergence code. \code{=0} means a success,
##'    and \code{=1} reaching the maximum number of iterations}
##' 
##' \item{ll}{log-likelihood value at convergence}
##' 
##' \item{maxgrad}{Maximum gradient value.}
##' 
##' \item{numiter}{number of iterations required by the algorithm}
##' 
##' @author Yong Wang <yongwang@@auckland.ac.nz>
##'  
##' @seealso \code{\link{cnm}}, \code{\link{nppois}}, \code{\link{disc}}.
##'  
##' @references
##' 
##' Wang, Y. (2007). On fast computation of the non-parametric maximum
##'  likelihood estimate of a mixing distribution. \emph{Journal of
##'  the Royal Statistical Society, Ser. B}, \bold{69}, 185-198.
##' 
##' Wang, Y. and Taylor, S. M. (2013). Efficient computation of
##'  nonparametric survival functions via a hierarchical mixture
##'  formulation. \emph{Statistics and Computing}, \bold{23}, 713-725.
##'  
##' @examples
##' 
##' x = rnppois(1000, disc(0:50))    # Poisson mixture
##' D = outer(x$v, 0:1000/10, dpois)
##' (r = hcnm(D, w=x$w))
##' disc(0:1000/10, r$p, collapse=TRUE)
##' 
##' cnm(x, init=list(mix=disc(0:1000/10)), model="p")
##' 
##' @export

hcnm = function(D, p0=NULL, w=1, maxit=1000, tol=1e-6,
                blockpar=NULL, recurs.maxit=2, compact=TRUE,
                depth=1, verbose=0) {
  nx = nrow(D)
  w = rep(w, length=nx)
  wr = sqrt(w)
  n = sum(w)
  converge = FALSE
  m = ncol(D)
  m1 = 1:m
  j = 1:m
  sj = m
  
  # nblocks = 1                      
  i = rowSums(D) == 1
  if(is.null(p0)) {
    ## Derive an initial p vector.
    j0 = colSums(D[i,,drop=FALSE]) > 0
    while(any(c(FALSE,(i <- rowSums(D[,j0,drop=FALSE])==0)))) {
      j0[which.max(colSums(D[i,,drop=FALSE]))] = TRUE
    }
    p = colSums(w * D) * j0
  }
  else { if(length(p <- p0) != m) stop("Argument 'p0' is the wrong length.") }
  p = p / sum(p)                   # mixing proportions
  P = drop(D %*% p)                # mixture probability of an observation
  ll = sum(w * log(P))             # log-likelihood given D
  evenstep = FALSE
  converge = FALSE
  
  for(iter in 1:maxit) {
    p.old = p
    ll.old = ll
    S = D / pmax(P, 1e-100)        # S = D / P
    g = colSums(w * S)
    if(verbose > 0) {
      cat("##### Iteration", iter-1, "#####\n")
      cat("Log-likelihood: ", signif(ll, 6), "\n")
    }
    if(verbose > 1) cat("Maximum gradient: ", signif(max(g) - n, 6), "\n")
    if(verbose > 2) {cat("Probability vector:\n"); print(p)} 

    if(compact) {
      j = p > 0
      if(depth == 1) {
        s = unique(c(1,m1[j],m))
        if (length(s) > 1)
          for(l in 2:length(s))
            j[s[l-1] + which.max(g[s[l-1]:s[l]]) - 1] = TRUE
      }
      sj = sum(j)
    }

    ## BW: matrix of block weights: sj rows, nblocks columns
    if(is.null(blockpar) || is.na(blockpar))
      iter.blockpar = ifelse(sj < 30, 0,
                             1 - log(max(20,10*log(sj/100)))/log(sj))
    else iter.blockpar = blockpar
    if(iter.blockpar == 0 | sj < 30) {   # find the number of blocks
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

    for(block in 1:nblocks) {            # updating inside each block
      jj = logical(m)
      jj[j] = BW[,block] > 0
      sjj = sum(jj)
      if (sjj > 1 && (delta <- sum(p.old[jj])) > 0) { # only block with pos. mass
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
      # if(ll >= ll.old && ll + ll.rise.alpha <= ll.old) {
      if(ll >= ll.old + ll.rise.alpha * .33) {  # backtracking successful
        p = p.alpha 
        break
      }
      if((alpha <- alpha * 0.5) < 1e-10) {      # backtracking failed
        p = p.old
        P = drop(D %*% p)
        ll = ll.old
        # converge = TRUE
        break
      }
      p.alpha = p.old + alpha * p.gap
      ll.rise.alpha = alpha * ll.rise.gap
    }
    # if(ll <= ll.old + tol) {
      ## p = p.alpha
      # converge = TRUE
      # break
    # }

    if (nblocks > 1) {            # updating masses of all blocks
      Q = sweep(BW,1,p[j],"*")
      q = colSums(Q)
      # print(dim(Q))
      # print(q)
      Q = sweep(D[,j] %*% Q, 2, q, "/")  
      if( any(q == 0) ) {         # no block can have zero probability
        ## print("A block has zero probability!")
        warning("A block has zero probability!")
      }
      else {
        res = hcnm(D=Q, p0=q, w=w, tol=tol, blockpar=iter.blockpar,
                   maxit=recurs.maxit, recurs.maxit=recurs.maxit, depth=depth+1)
        if (res$ll > ll) {
          p[j] = p[j] * (BW %*% (res$p / q))
          P = drop(D %*% p)
          ll = sum(w * log(P))
        }
      }
    }
    if(iter > 2) if( ll <= ll.old + tol ) {converge = TRUE; break}
    evenstep = !evenstep
  }
  list(p=p, convergence=converge, ll=ll,
       maxgrad=max(crossprod(w/P, D))-n, numiter=iter)
}

##' @title Maximum Likelihood Estimation of a Semiparametric Mixture Model
##' 
##' @usage
##' cnmms(x, init=NULL, maxit=1000, model=c("spmle","npmle"), tol=1e-6,
##'       grid=100, kmax=Inf, plot=c("null", "gradient", "probability"),
##'       verbose=0)
##' cnmpl(x, init=NULL, tol=1e-6, tol.npmle=tol*1e-4, grid=100, maxit=1000,
##'       plot=c("null", "gradient", "probability"), verbose=0)
##' cnmap(x, init=NULL, maxit=1000, tol=1e-6, grid=100, plot=c("null",
##'       "gradient"), verbose=0)
##'  
##' @description Functions \code{cnmms}, \code{cnmpl} and \code{cnmap}
##'   can be used to compute the maximum likelihood estimate of a
##'   semiparametric mixture model that has a one-dimensional mixing
##'   parameter. The types of mixture models that can be computed
##'   include finite, nonparametric and semiparametric ones.
##' 
##' Function \code{cnmms} can also be used to compute the maximum
##' likelihood estimate of a finite or nonparametric mixture model.
##' 
##' A finite mixture model has a density of the form
##' 
##' \deqn{f(x; \pi, \theta, \beta) = \sum_{j=1}^k \pi_j f(x; \theta_j,
##' \beta).}{f(x; pi, theta, beta) = sum_{j=1}^k pi_j f(x; theta_j, beta),}
##' 
##' where \eqn{pi_j \ge 0}{pi_j >= 0} and \eqn{\sum_{j=1}^k pi_j }{sum_{j=1}^k
##' pi_j =1}\eqn{ =1}{sum_{j=1}^k pi_j =1}.
##' 
##' A nonparametric mixture model has a density of the form
##' 
##' \deqn{f(x; G) = \int f(x; \theta) d G(\theta),}{f(x; G) = Integral
##' f(x; theta) d G(theta),} where \eqn{G} is a mixing distribution
##' that is completely unspecified. The maximum likelihood estimate of
##' the nonparametric \eqn{G}, or the NPMLE of $\eqn{G}, is known to
##' be a discrete distribution function.
##' 
##' A semiparametric mixture model has a density of the form
##' 
##' \deqn{f(x; G, \beta) = \int f(x; \theta, \beta) d G(\theta),}{f(x; G, beta)
##' = Int f(x; theta, beta) d G(theta),}
##' 
##' where \eqn{G} is a mixing distribution that is completely
##' unspecified and \eqn{\beta}{beta} is the structural parameter.
##' 
##' Of the three functions, \code{cnmms} is recommended for most
##' problems; see Wang (2010).
##' 
##' Functions \code{cnmms}, \code{cnmpl} and \code{cnmap} implement
##' the algorithms CNM-MS, CNM-PL and CNM-AP that are described in
##' Wang (2010).  Their implementations are generic using S3
##' object-oriented programming, in the sense that they can work for
##' an arbitrary family of mixture models that is defined by the
##' user. The user, however, needs to supply the implementations of
##' the following functions for their self-defined family of mixture
##' models, as they are needed internally by the functions above:
##' 
##' \code{initial(x, beta, mix, kmax)}
##' 
##' \code{valid(x, beta)}
##' 
##' \code{logd(x, beta, pt, which)}
##' 
##' \code{gridpoints(x, beta, grid)}
##' 
##' \code{suppspace(x, beta)}
##' 
##' \code{length(x)}
##' 
##' \code{print(x, ...)}
##' 
##' \code{weight(x, ...)}
##' 
##' While not needed by the algorithms, one may also implement
##' 
##' \code{plot(x, mix, beta, ...)}
##' 
##' so that the fitted model can be shown graphically in a way that the user
##' desires.
##' 
##' For creating a new class, the user may consult the implementations
##' of these functions for the families of mixture models included in
##' the package, e.g., \code{cvp} and \code{mlogit}.
##' 
##' @aliases cnmms cnmpl cnmap
##'  
##' @param x a data object of some class that can be defined fully by
##'   the user
##' @param init list of user-provided initial values for the mixing
##'   distribution \code{mix} and the structural parameter \code{beta}
##' @param model the type of model that is to estimated:
##'   non-parametric MLE (\code{npmle}) or semi-parametric MLE
##'   (\code{spmle}).
##' @param maxit maximum number of iterations
##' @param tol a tolerance value that is used to terminate an
##'   algorithm.  Specifically, the algorithm is terminated, if the
##'   relative increase of the log-likelihood value after an iteration
##'   is less than \code{tol}. If an algorithm converges rapidly
##'   enough, then \code{-log10(tol)} is roughly the number of
##'   accurate digits in log-likelihood.
##' @param tol.npmle a tolerance value that is used to terminate the
##'   computing of the NPMLE internally.
##' @param grid number of grid points that are used by the algorithm
##'   to locate all the local maxima of the gradient function. A
##'   larger number increases the chance of locating all local maxima,
##'   at the expense of an increased computational cost. The locations
##'   of the grid points are determined by the function
##'   \code{gridpoints} provided by each individual mixture family,
##'   and they do not have to be equally spaced. If needed, an
##'   individual \code{gridpoints} function may return a different
##'   number of grid points than specified by \code{grid}.
##' @param kmax upper bound on the number of support points. This is
##'   particularly useful for fitting a finite mixture model.
##' @param plot whether a plot is produced at each iteration. Useful
##'   for monitoring the convergence of the algorithm. If \code{null},
##'   no plot is produced. If \code{gradient}, it plots the gradient
##'   curves and if \code{probability}, the \code{plot} function
##'   defined by the user of the class is used.
##' @param verbose verbosity level for printing intermediate results
##'   in each iteration, including none (= 0), the log-likelihood
##'   value (= 1), the maximum gradient (= 2), the support points of
##'   the mixing distribution (= 3), the mixing proportions (= 4), and
##'   if available, the value of the structural parameter beta (= 5).
##'  
##' @return
##' 
##' \item{family}{the class of the mixture family that is used to fit
##' to the data.}
##' 
##' \item{num.iterations}{Number of iterations required by the algorithm}
##' 
##' \item{grad}{For \code{cnmms}, it contains the values of the
##' gradient function at the support points and the first derivatives
##' of the log-likelihood with respect to theta and beta. For
##' \code{cnmpl}, it contains only the first derivatives of the
##' log-likelihood with respect to beta. For \code{cnmap}, it contains
##' only the gradient of \code{beta}.}
##' 
##' \item{max.gradient}{Maximum value of the gradient function,
##' evaluated at the beginning of the final iteration. It is only
##' given by function \code{cnmap}.}
##' 
##' \item{convergence}{convergence code. \code{=0} means a success,
##' and \code{=1} reaching the maximum number of iterations}
##' 
##' \item{ll}{log-likelihood value at convergence}
##' 
##' \item{mix}{MLE of the mixing distribution, being an object of the
##' class \code{disc} for discrete distributions}
##' 
##' \item{beta}{MLE of the structural parameter}
##'  
##' @author Yong Wang <yongwang@@auckland.ac.nz>
##'  
##' @seealso \code{\link[lsei]{nnls}}, \code{\link{cnm}},
##'    \code{\link{cvp}}, \code{\link{cvps}}, \code{\link{mlogit}}.
##'  
##' @references
##' 
##' Wang, Y. (2007). On fast computation of the non-parametric maximum
##' likelihood estimate of a mixing distribution. \emph{Journal of the Royal
##' Statistical Society, Ser. B}, \bold{69}, 185-198.
##' 
##' Wang, Y. (2010). Maximum likelihood computation for fitting semiparametric
##' mixture models. \emph{Statistics and Computing}, \bold{20}, 75-86
##'  
##' @examples
##' 
##' ## Compute the MLE of a finite mixture
##' x = rnpnorm(100, disc(c(0,4), c(0.7,0.3)), sd=1)
##' for(k in 1:6) plot(cnmms(x, kmax=k), x, add=(k>1), comp="null", col=k+1,
##'                    main="Finite Normal Mixtures")
##' legend("topright", 0.3, leg=paste0("k = ",1:6), lty=1, lwd=2, col=2:7)
##' 
##' ## Compute a semiparametric MLE
##' # Common variance problem 
##' x = rcvps(k=50, ni=5:10, mu=c(0,4), pr=c(0.7,0.3), sd=3)
##' cnmms(x)              # CNM-MS algorithm
##' cnmpl(x)              # CNM-PL algorithm
##' cnmap(x)              # CNM-AP algorithm
##' 
##' # Logistic regression with a random intercept
##' x = rmlogit(k=30, gi=3:5, ni=6:10, pt=c(0,4), pr=c(0.7,0.3),
##'             beta=c(0,3))
##' cnmms(x)
##' 
##' data(toxo)            # k = 136
##' cnmms(mlogit(toxo))
##' 
##' @export cnmms
##' @export cnmpl
##' @export cnmap

cnmms = function(x, init=NULL, maxit=1000,
      model=c("spmle","npmle"), tol=1e-6, grid=100, kmax=Inf,
      plot=c("null", "gradient", "probability"), verbose=0) {
  if(is.null(class(x))) stop("Data belongs to no class")
  plot = match.arg(plot)
  model = match.arg(model)
  init = initial0(x, init, kmax=kmax)
  beta = init$beta
  if(is.null(beta) || any(is.na(beta))) model = "npmle"
  nb = length(beta)
  mix = init$mix
  attr(mix, "ll") = init$ll
  gradient = "Not computed"
  switch(plot,
         "gradient" = plotgrad(x, mix, beta, pch=1),
         "probability" = plot(x, mix, beta) )
  ## ll = -Inf
  if(maxit == 0) 
    return( structure(list(family=class(x)[1], mix=mix, beta=beta,
                           num.iterations=0, ll=attr(mix, "ll")[1],
                           max.gradient=max(maxgrad(x, beta, mix, tol=-5)$grad),
                           convergence=1), 
                      class="nspmix") )
  gmax = Inf
  convergence = 1
  for(i in 1:maxit) {
    mix1 = mix
    switch(plot,
           "gradient" = plotgrad(x, mix, beta, pch=1),
           "probability" = plot(x, mix, beta) )
    if(length(mix$pt) < kmax) {
      gp = gridpoints(x, beta, grid)
      g = maxgrad(x, beta, mix, grid=gp, tol=-Inf)
      gmax = g$gmax
      kpt = min(kmax - length(mix$pt), length(g$pt))
      jpt = order(g$grad, decreasing=TRUE)
      mix = disc(c(mix$pt,g$pt[jpt][1:kpt]), c(mix$pr,rep(0,kpt)))
      if(plot=="gradient") points(g$pt, g$grad, col=2, pch=20)
    }
    if(verbose > 0)
      verb(verbose, attr(mix1, "ll")[1], mix1, beta, gmax, i-1)
    lpt = logd(x, beta, mix$pt, which=c(1,0,0))$ld
    S = pmin(exp(lpt - attr(attr(mix1,"ll"),"logd")), 1e100)
    w = weight(x, beta)
    wr = sqrt(w)
    grad.support = colSums(w * (S - 1))
    r = pnnls(wr * S, wr * 2, sum=1)
    sol = r$x / sum(r$x)
    r = lsch(mix, beta, disc(mix$pt,sol,sort=FALSE), beta, x, which=c(1,0,0))
    mix = r$mix
    attr(mix, "ll") = loglik(mix, x, beta, attr=TRUE)
    mix = collapse0(mix, beta, x, tol=0)
    if(max(grad.support) < 1e10) {    # 1e20
      r = switch(model,
        spmle = bfgs(mix, beta, x, which=c(1,1,1)),
        npmle = bfgs(mix, beta, x, which=c(1,1,0))   )
      beta = r$beta
      mix = collapse0(r$mix, beta, x, tol=tol*0.01)
    }
    if((length(mix$pt) == kmax | gmax < 1e-3) &
       attr(mix,"ll")[1] <= attr(mix1,"ll")[1] + tol)
         {convergence = 0; break}
  }
  mix = collapse0(mix, beta, x,
                        tol=max(tol*.1, abs(attr(mix, "ll")[1])*1e-16))
  m = length(mix$pt)
  if(m < length(r$mix$pt)) {
    d = dll0(x, mix, beta, which=c(0,1,1,1))
    grad = c(d$dp, d$dt, d$db)
    names(grad) = c(paste0("pr",1:m), paste0("pt",1:m), 
           if(is.null(d$db)) NULL else paste0("beta",1:length(beta)) )
  }
  else grad = r$grad
  grad[1:m] = grad[1:m] - sum(rep(weight(x, beta), len=length(x)))
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
  init = initial0(x, init)
  beta = init$beta
  mix = init$mix
  nb = length(beta)
  if(is.null(beta) || any(is.na(beta)))
    stop("No proper initial value for beta. Use cnm() for NPMLE problems.")
  switch(plot,
         "gradient" = plotgrad(x, mix, beta, pch=1),
         "probability" = plot(x, mix, beta) )
  r = pll(beta, mix, x, tol=tol.npmle, grid=grid)
  D = diag(-1, nrow=nb)
  convergence = 1
  num.npmle = 1
  if(maxit == 0) {
    r = list(family=class(x)[1], mix=mix, beta=beta, num.iterations=0,
        ll=loglik(mix, x, beta, attr=TRUE),
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
    if(verbose > 0) verb(verbose, r$ll[1], r$mix, r$beta, r$max.grad, i)
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
  init = initial0(x, init)
  beta = init$beta
  if(is.null(beta) || any(is.na(beta)))
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
    gmax = g$gmax
    if(verbose > 0)
      verb(verbose, attr(mix, "ll")[1], mix, beta, gmax, i-1)
    mix = disc(c(mix$pt,g$pt), c(mix$pr,rep(0,length(g$pt))))
    lpt = logd(x, beta, mix$pt, which=c(1,0,0))$ld
    ## ll1 = attr(mix1, "ll")
    # logd.mix1 = log(attr(ll1,"dmix")) + attr(ll1,"ma")
    ## S = pmin(exp(lpt - logd.mix1), 1e100)
    S = pmin(exp(lpt - attr(attr(mix1,"ll"),"logd")), 1e100)
    wr = sqrt(weight(x, beta))
    r = pnnls(wr * S, wr * 2, sum=1)
    sol = r$x / sum(r$x)
    r = lsch(mix, beta, disc(mix$pt,sol,sort=FALSE), beta, x, which=c(1,0,0))
    mix = collapse0(r$mix, beta, x, tol=0)
    r = bfgs(mix, beta, x, which=c(0,0,1))
    beta = r$beta
    mix = collapse0(r$mix, beta, x, tol=0)
    if( attr(mix, "ll")[1] <= attr(mix1, "ll")[1] + tol ) {convergence = 0; break} 
  }
  mix = collapse0(mix, beta, x,
                        tol=max(tol*.1, abs(attr(mix, "ll")[1])*1e-16))
  m = length(mix$pt)
  d = dll0(x, mix, beta, which=c(0,1,1,1))
  grad = c(d$dp, d$dt, d$db)
  names(grad) = c(paste0("pr",1:m), paste0("pt",1:m), 
                  if(is.null(d$db)) NULL else paste0("beta",1:length(beta)) )
  grad[1:m] = grad[1:m] - sum(rep(weight(x, beta), len=length(x)))
  ll = attr(mix,"ll")[1]
  attr(mix,"ll") = NULL
  structure(list(family=class(x)[1], num.iterations=i, grad=grad,
                 max.gradient=gmax, convergence=convergence,
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
  dl = dll0(x, mix, beta, which=c(1,which), ind=TRUE)
  w = weight(x, beta)
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
      lambda = sum(r$grad * rowSums(D1)) / sum(D11)
      d1 = d1 + lambda * rowSums(D1)
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
            lambda = sum(r$grad * rowSums(D1)) / sum(D11)
            d1 = d1 + lambda * rowSums(D1)
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
      if(which[1]) par2[1:k1] else r$mix$pr, sort=FALSE)
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
  dl1 = dll0(x, mix1, beta1, which=c(1,which))
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
  if(d1.norm == 0 || g1s <= 0) {    # wrong direction: should never happan
    ## cat("Warning in lsch(): d1.norm =", d1.norm, " g1s =", g1s, "\n")
    return( list(mix=mix1, beta=beta1, grad=grad1, ll=ll1, convergence=3) )
  }
  
  ## bracketing phase
  a = 0
  b = 1
  if(which[1] && any(mix2$pr == 0)) brkt = FALSE
  for(i in 1:maxit) {
    for(j in 1:1000) {
      m = disc( (1-b) * mix1$pt + b * mix2$pt, (1-b) * mix1$pr + b * mix2$pr)
      m$pt[je] = mix1$pt[je]
      # beta = if(is.null(beta1)) NULL else (1-b) * beta1 + b * beta2
      beta = if(which[3]) (1-b) * beta1 + b * beta2 else beta1
      if( valid0(x, beta, m) ) break
      brkt = FALSE
      b = 0.5 * a + 0.5 * b
    }
    if(j == 1000) stop("Can not produce valid interior point in lsch()")
    dl = dll0(x, m, beta, which=c(1,which))
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
  
  ## sectioning phase
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
    dl = dll0(x, m, beta, which=c(1,which))
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
      if( valid0(x, beta, r$mix) ) break
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


##' @title Log-likelihood value of a mixture
##' 
##' @description Computes the log-likelihood value
##' 
##' \code{x} must belong to a mixture family, as specified by its class.
##' 
##' @param x a data object of a mixture model class.
##' @param mix a discrete distribution, as defined by class
##'   \code{disc}.
##' @param beta the structural parameter, if any.
##' @param attr =FALSE, by default. If TRUE, also returns attributes
##'   "dmix" and "logd"
##'  
##' @return the log-likelihood value.
##'  
##' @author Yong Wang <yongwang@@auckland.ac.nz>
##'  
##' @seealso \code{\link{cnm}}, \code{\link{cnmms}},
##'   \code{\link{npnorm}}, \code{\link{nppois}}, \code{\link{disc}},
##'  
##' @references
##' 
##' Wang, Y. (2007). On fast computation of the non-parametric maximum
##' likelihood estimate of a mixing distribution. \emph{Journal of the
##' Royal Statistical Society, Ser. B}, \bold{69}, 185-198.
##' 
##' Wang, Y. (2010). Maximum likelihood computation for fitting
##' semiparametric mixture models. \emph{Statistics and Computing},
##' \bold{20}, 75-86
##'  
##' @examples
##' 
##' ## Poisson mixture
##' mix0 = disc(c(1,4), c(0.7,0.3))
##' x = rnppois(10, mix0)
##' loglik(mix0, x)
##' 
##' ## Normal mixture
##' x = rnpnorm(10, mix0, sd=2)
##' loglik(mix0, x, 2)
##' 
##' @export 

loglik = function(mix, x, beta=NULL, attr=FALSE) {
  ld = logd(x, beta, mix$pt, which=c(1,0,0))$ld
  ma = matMaxs(ld)
  dmix = drop(exp(ld - ma) %*% mix$pr) + 1e-100
  logd = log(dmix) + ma
  ll = llex(x, beta, mix) + sum(weight(x, beta) * logd)
  if(attr) {
    attr(ll, "dmix") = dmix
    attr(ll, "logd") = logd    # log(mixture density)
  }
  ll
}

## density.snpmle = function(x, beta, mix) {
##   rowSums(exp(logd(x, beta, mix$pt, which=c(1,0,0))$ld +
##       rep(log(mix$pr), rep(length(x), length(mix$pr)))))
## }

# Profile log-likelihood

pll = function(beta, mix0, x, tol=1e-15, grid=100, ...) {
  r = cnm(x, init=list(beta=beta, mix=mix0), tol=tol, grid=grid, maxit=50, ...)
  mix = r$mix
  dl = logd(x, beta, mix$pt, which=c(1,1,0))
  lpt = dl$ld
  ma = matMaxs(lpt)
  dmix = drop(exp(lpt - ma) %*% mix$pr) + 1e-100
  D = pmin(exp(lpt - ma), 1e100)
  p = weight(x, beta) * D/dmix * rep(mix$pr, rep(nrow(D), length(mix$pr)))
  db = apply(sweep(dl$db, c(1,2), p, "*"), c(1,3), sum)
  g = llexdb(x,beta,mix) + colSums(db)
  list(ll=r$ll, beta=r$beta, mix=mix, grad=g, max.gradient=r$max.gradient,
       num.iterations=r$num.iter)
}

## Compute all local maxima of the gradient function using a hybrid
## secant-bisection method. No need for the second derivative of the
## gradient function.

## O(n*m) + #iteration * O(n*mhat)
## n - sample size
## m - grid points
## mhat - number of local maxima

maxgrad = function(x, beta, mix, grid=100, tol=-Inf, maxit=100) {
  if(length(grid) == 1) grid = gridpoints(x, beta, grid)
  dg = grad(grid, x, beta, mix, order=1)$d1
  j = is.na(dg)
  grid = grid[!j]
  tol.pt = 1e-14 * diff(range(grid))
  dg = dg[!j]
  np = length(grid)
  jmax = which(dg[-np] > 0 & dg[-1] < 0)
  # if( length(jmax) < 1 )
  #   return( list(pt=NULL, grad=NULL, num.iterations=0, gmax=NULL) )
  pt = NULL
  if(length(jmax) != 0) {
    left = grid[jmax]
    right = grid[jmax+1]
    pt = (left + right) * .5
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
  gmax = max(g)
  names(pt) = names(g) = NULL
  j = g >= tol
  list(pt=pt[j], grad=g[j], num.iterations=i, gmax=gmax)
}

# gradient function, with O(n*m)

grad = function(pt, x, beta, mix, order=0) {
  w = weight(x, beta)
  if(length(w) != length(x)) w = rep(w, length=length(x))
  if(is.null(attr(mix, "ll")) || is.null(attr(attr(mix, "ll"), "dmix")))
    attr(mix, "ll") = loglik(mix, x, beta, attr=TRUE)
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



##' Plot the Gradient Function
##' 
##' 
##' Function \code{plotgrad} plots the gradient function or its first
##' derivative of a nonparametric mixture.
##' 
##' 
##' \code{data} must belong to a mixture family, as specified by its class.
##' 
##' The support points are shown on the horizontal line of gradient
##' 0. The vertical lines going downwards at the support points are
##' proportional to the mixing proportions at these points.
##' 
##' @param x a data object of a mixture model class.
##' @param mix an object of class 'disc', for a discrete mixing
##'   distribution.
##' @param beta the structural parameter.
##' @param len number of points used to plot the smooth curve.
##' @param order the order of the derivative of the gradient function
##'   to be plotted. If 0, it is the gradient function itself.
##' @param col color for the curve.
##' @param col2 color for the support points.
##' @param add if \code{FALSE}, create a new plot; if \code{TRUE}, add
##'   the curve and points to the current one.
##' @param main,xlab,ylab,cex,pch,lwd,xlim,ylim arguments for
##'   graphical parameters (see \code{par}).
##' @param ... arguments passed on to function \code{plot}.
##'  
##' @author Yong Wang <yongwang@@auckland.ac.nz>
##'  
##' @seealso \code{\link{plot.nspmix}}, \code{\link[lsei]{nnls}},
##'    \code{\link{cnm}}, \code{\link{cnmms}}, \code{\link{npnorm}},
##'    \code{\link{nppois}}.
##'  
##' @references
##' 
##' Wang, Y. (2007). On fast computation of the non-parametric maximum
##' likelihood estimate of a mixing distribution. \emph{Journal of the
##' Royal Statistical Society, Ser. B}, \bold{69}, 185-198.
##' 
##' Wang, Y. (2010). Maximum likelihood computation for fitting
##' semiparametric mixture models. \emph{Statistics and Computing},
##' \bold{20}, 75-86
##'  
##' @examples
##' 
##' ## Poisson mixture
##' x = rnppois(200, disc(c(1,4), c(0.7,0.3)))
##' r = cnm(x)
##' plotgrad(x, r$mix)
##' 
##' ## Normal mixture
##' x = rnpnorm(200, disc(c(0,4), c(0.3,0.7)), sd=1)
##' r = cnm(x, init=list(beta=0.5))   # sd = 0.5
##' plotgrad(x, r$mix, r$beta)
##' 
##' @export

plotgrad = function(x, mix, beta, len=500, order=0, col=4,
                    col2=2, add=FALSE,
                    main=paste0("Class: ",class(x)),
                    xlab=expression(theta),
                    ylab=paste0("Gradient (order = ",order,")"), cex=1,
                    pch=1, lwd=1, xlim, ylim, ...) {
  if(order > 1) stop("order > 1")
  if( missing(xlim) ) xlim = range(gridpoints(x, beta, grid=5))
  pt = seq(xlim[1], xlim[2], len=len)
  if("disc" %in% class(mix))
    pt = sort(c(mix$pt[mix$pt>xlim[1] & mix$pt<xlim[2]], pt))
  g = grad(pt, x, beta, mix, order=order)[[1]]
  if(missing(ylim)) ylim = range(0,g,na.rm=TRUE)
  if(!add) plot(xlim, c(0,0), type="l", col="black", ylim=ylim, main=main,
       xlab=xlab, ylab=ylab, cex = cex, cex.axis = cex, cex.lab = cex, ...)
  lines(pt, g, col=col, lwd=lwd)
  if("disc" %in% class(mix)) {
    j = mix$pr > 0 & mix$pt > xlim[1] & mix$pt < xlim[2]
    points(mix$pt[j], rep(0,length(mix$pt[j])), pch=pch, col=col2)
    lines(mix$pt[j], mix$pr[j]*min(g), type="h", col=col2)
  }
  invisible(list(pt=pt,g=g))
}

##' @title Initialisation
##'
##' @description The functions examines whether the initial values are
##'   proper. If not, proper ones are provided, by employing the
##'   function "initial" provided by the class.
##'
##' @param x an object of a class for data.
##' @param init a list with initial values for beta and mix (as in the
##'   output of \code{initial})
##' @param kmax the maximum allowed number of support points used.
##' 
##' @author Yong Wang <yongwang@@auckland.ac.nz>
##'
##' @return
##'
##' \item{beta}{an initialized value of beta}
##' 
##' \item{mix}{an initialised or updated object of class \code{disc}
##' for the mixing distribution.}
##'
##' @export 

initial0 = function(x, init=NULL, kmax=NULL) {
  if(!is.null(kmax) && kmax == Inf) kmax = NULL
  init = initial(x, init$beta, init$mix, kmax=kmax)
  init$ll = loglik(init$mix, x, init$beta, attr=TRUE)
  if(! valid0(x, init$beta, init$mix)) stop("Invalid initial values!")
  init
}

# Compute the first derivatives of the log-likelihood function of pi,
# theta and beta

# INPUT:
# 
# which     A 4-vector having values either 0 and 1, indicating which
#           derivatives to be computed, in the order of pi, theta and beta
# ind       = TRUE, return the derivatives with respect to individual
#                   observations (which may have weights)
#           = FALSE, return the derivatives
#
# OUTPUT:   a list with ll, dp, dt, db, as specified by which

dll0 = function(x, mix, beta, which=c(1,0,0,0), ind=FALSE) {
  w = weight(x, beta)
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

valid0 = function(x, beta, mix) {
  bs = suppspace(x, beta)
  bs = bs + 1e-14 * c(-1,1) * abs(bs)             # tolerance
  valid(x, beta, mix) && all(mix$pr >= 0, mix$pt >= bs[1], mix$pt <= bs[2])
}

collapse0 = function(mix, beta, x, tol=1e-15) {
  j = mix$pr == 0
  if(any(j)) {
    mix = disc(mix$pt[!j], mix$pr[!j])
    attr(mix, "ll") = loglik(mix, x, beta, attr=TRUE)
  }
  else { if( is.null(attr(mix, "ll")) || is.null(attr(attr(mix, "ll"), "dmix")) )
           attr(mix, "ll") = loglik(mix, x, beta, attr=TRUE) }
  
  if(tol > 0) {
    repeat {
      if(length(mix$pt) <= 1) break
      prec = max(10 * min(diff(mix$pt)), 1e-100)   ## one digit at a time
      ## if(prec > diff(range(gridpoints(x, beta, grid=2))) * .01) break
      mixt = unique(mix, prec=c(prec,0))
      attr(mixt, "ll") = loglik(mixt, x, beta, attr=TRUE)
      # if(attr(mixt,"ll")[1] >= attr(mix,"ll")[1] - tol * abs(attr(mix,"ll")[1]))
      if(attr(mixt,"ll")[1] >= attr(mix,"ll")[1] - tol) mix = mixt  
      else break
    }
  }
  mix
}

verb = function(verbose, ll, mix, beta, maxdg, i) {
  if(verbose > 0)
    cat("## Iteration ", i, ":\n   Log-likelihood = ", as.character(ll), "\n", sep="")
  if(verbose > 1) cat("   Max gradient =", maxdg, "\n")
  if(verbose > 2) cat("   theta =", signif(mix$pt,6), "\n")
  if(verbose > 3) cat("   Proportions =", signif(mix$pr,6), "\n")
  if(verbose > 4)
    cat("   beta =", if(is.null(beta)) "NULL" else signif(beta,6), "\n")
}

##' @title Plots a function for an object of class \code{nspmix}
##' 
##' @description Plots a function for the object of class
##'   \code{nspmix}, currently either using the plot function of the
##'   class or plotting the gradient curve (or its first derivative)
##' 
##' \code{data} must belong to a mixture family, as specified by its class.
##' 
##' @param data a data object of a mixture model class.
##' @param x an object of class \code{nspmix}
##' @param type if \code{probability}, use the plot function of the
##'   data class; if \code{gradient}, plot the gradient curve
##' @param ... arguments passed on to function \code{plot}.
##'  
##' @author Yong Wang <yongwang@@auckland.ac.nz>
##'  
##' @seealso \code{\link{plot.nspmix}}, \code{\link[lsei]{nnls}},
##'   \code{\link{cnm}}, \code{\link{cnmms}}, \code{\link{npnorm}},
##'   \code{\link{nppois}}.
##'  
##' @references
##' 
##' Wang, Y. (2007). On fast computation of the non-parametric maximum
##' likelihood estimate of a mixing distribution. \emph{Journal of the
##' Royal Statistical Society, Ser. B}, \bold{69}, 185-198.
##' 
##' Wang, Y. (2010). Maximum likelihood computation for fitting
##' semiparametric mixture models. \emph{Statistics and Computing},
##' \bold{20}, 75-86
##'  
##' @examples
##' 
##' ## Poisson mixture
##' x = rnppois(200, disc(c(1,4), c(0.7,0.3)))
##' plot(cnm(x), x)
##' 
##' ## Normal mixture
##' x = rnpnorm(200, disc(c(0,4), c(0.3,0.7)), sd=1)
##' r = cnm(x, init=list(beta=0.5))   # sd = 0.5
##' plot(r, x)
##' plot(r, x, type="g")
##' plot(r, x, type="g", order=1)
##' 
##' @export

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

##' @title Density function of a mixture distribution
##' 
##' @description Computes the density or their logarithmic values of a
##'   mixture distribution, where the component family depends on the
##'   class of \code{x}.
##' 
##' \code{x} must belong to a mixture family, as specified by its class.
##' 
##' @param x a data object of a mixture model class.
##' @param mix a discrete distribution, as defined by class
##'   \code{disc}.
##' @param beta the structural parameter, if any.
##' @param log if \code{TRUE}, computes the log-values, or else just
##'   the density values.
##'  
##' @author Yong Wang <yongwang@@auckland.ac.nz>
##'  
##' @seealso \code{\link{cnm}}, \code{\link{cnmms}},
##'   \code{\link{npnorm}}, \code{\link{nppois}}, \code{\link{disc}},
##'  
##' @references
##' 
##' Wang, Y. (2007). On fast computation of the non-parametric maximum
##' likelihood estimate of a mixing distribution. \emph{Journal of the
##' Royal Statistical Society, Ser. B}, \bold{69}, 185-198.
##' 
##' Wang, Y. (2010). Maximum likelihood computation for fitting
##' semiparametric mixture models. \emph{Statistics and Computing},
##' \bold{20}, 75-86
##'  
##' @examples
##' 
##' ## Poisson mixture
##' mix0 = disc(c(1,4), c(0.7,0.3))
##' x = rnppois(10, mix0)
##' dmix(x, mix0)
##' dmix(x, mix0, log=TRUE)
##' 
##' ## Normal mixture
##' x = rnpnorm(10, mix0, sd=1)
##' dmix(x, mix0, 1)
##' dmix(x, mix0, 1, log=TRUE)
##' dmix(x, mix0, 0.5, log=TRUE)
##' 
##' @export

dmix = function(x, mix, beta=NULL, log=FALSE) {
  ld = logd(x, beta, mix$pt, which=c(1,0,0))$ld + 
    rep(log(mix$pr), rep(length(x), length(mix$pr)))
  if(log) {
    ma = matMaxs(ld)
    ma + log(rowSums(exp(ld - ma)))
  }
  else rowSums(exp(ld))
}
