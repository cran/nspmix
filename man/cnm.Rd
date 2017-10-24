\name{cnm}
\alias{cnm}

\title{Maximum Likelihood Estimation of a Nonparametric Mixture Model}

\description{
  
  Function \code{cnm} can be used to compute the maximum likelihood
  estimate of a nonparametric mixing distribution (NPMLE) that has a
  one-dimensional mixing parameter. or simply the mixing proportions
  with support points held fixed.
  
}

\usage{

cnm(x, init=NULL, model=c("npmle","proportions"), maxit=100, tol=1e-6,
    grid=100, plot=c("null", "gradient", "probability"), verbose=0)

}

\arguments{
  
  \item{x}{a data object of some class that is fully defined by the
  user. The user needs to supply certain functions as described below.}
  
  \item{init}{list of user-provided initial values for the mixing
  distribution \code{mix} and the structural parameter \code{beta}.} 
  
  \item{model}{the type of model that is to estimated: the non-parametric MLE
    (if \code{npmle}), or mixing proportions only (if \code{proportions}). }
  
  \item{maxit}{maximum number of iterations.}
  
  \item{tol}{a tolerance value needed to terminate an
    algorithm. Specifically, the algorithm is terminated, if the
    relative increase of the log-likelihood value after an iteration is
    less than \code{tol}. If an algorithm converges rapidly enough, then
    \code{-log10(tol)} is roughly the number of accurate digits in
    log-likelihood.}

  \item{grid}{number of grid points that are used by the algorithm to
    locate all the local maxima of the gradient function. A larger
    number increases the chance of locating all local maxima, at the
    expense of an increased computational cost. The locations of the
    grid points are determined by the function \code{gridpoints}
    provided by each individual mixture family, and they do not have to
    be equally spaced. If needed, a \code{gridpoints} function may
    choose to return a different number of grid points than specified by
    \code{grid}.}

  \item{plot}{whether a plot is produced at each iteration. Useful for
    monitoring the convergence of the algorithm. If \code{="null"}, no plot
    is produced. If \code{="gradient"}, it plots the gradient curves and if
    \code{="probability"}, the \code{plot} function defined by the user for
    the class is used.}

  \item{verbose}{verbosity level for printing intermediate results in
    each iteration, including none (= 0), the log-likelihood value (=
    1), the maximum gradient (= 2), the support points of the mixing
    distribution (= 3), the mixing proportions (= 4), and if available,
    the value of the structural parameter beta (= 5).}
  
}

\details{
  
  A finite mixture model has a density of the form

  \deqn{f(x; \pi, \theta, \beta) = \sum_{j=1}^k \pi_j f(x; \theta_j, \beta).}{f(x; pi, theta, beta) = sum_{j=1}^k pi_j f(x; theta_j, beta),}
  
  where \eqn{pi_j \ge 0}{pi_j >= 0} and \eqn{\sum_{j=1}^k pi_j
  =1}{sum_{j=1}^k pi_j =1}. 

  A nonparametric mixture model has a density of the form

  \deqn{f(x; G) = \int f(x; \theta) d G(\theta),}{f(x; G) = Integral f(x; theta) d G(theta),}
  where \eqn{G} is a mixing distribution that is completely
  unspecified. The maximum likelihood estimate of the
  nonparametric \eqn{G}, or the NPMLE of $\eqn{G}, is known to be a
  discrete distribution function.
  
  Function \code{cnm} implements the CNM algorithm that is proposed in
  Wang (2007) and the hierarchical CNM algorithm of Wang and Taylor
  (2013). The implementation is generic using S3 object-oriented
  programming, in the sense that it works for an arbitrary family of
  mixture models defined by the user. The user, however, needs to supply
  the implementations of the following functions for their self-defined
  family of mixture models, as they are needed internally by function
  \code{cnm}:

  \code{initial(x, beta, mix, kmax)}
  
  \code{valid(x, beta)}

  \code{logd(x, beta, pt, which)}

  \code{gridpoints(x, beta, grid)}
  
  \code{suppspace(x, beta)}

  \code{length(x)}

  \code{print(x, ...)}

  \code{weights(x, ...)}

  While not needed by the algorithm for finding the solution, one may
  also implement
  
  \code{plot(x, mix, beta, ...)}

  so that the fitted model can be shown graphically in a user-defined
  way. Inside \code{cnm}, it is used when \code{plot="probability"} so
  that the convergence of the algorithm can be graphically monitored.

  For creating a new class, the user may consult the implementations of
  these functions for the families of mixture models included in the
  package, e.g., \code{npnorm} and \code{nppois}.

}

\value{
  
  \item{family}{the name of the mixture family that is used to fit to
                the data.}
  
  \item{num.iterations}{Number of iterations required by the algorithm}
  
  \item{max.gradient}{Maximum value of the gradient function, evaluated
    at the beginning of the final iteration}
  
  \item{convergence}{convergence code. \code{=0} means a success, and 
  \code{=1} reaching the maximum number of iterations}
  
  \item{ll}{log-likelihood value at convergence}
  
  \item{mix}{MLE of the mixing distribution, being an object of the
    class \code{disc} for discrete distributions.}
  
  \item{beta}{The structural parameter, the value of which does not
    change during the computation.}
  
}

\author{ Yong Wang <yongwang@auckland.ac.nz>}

\references{

  Wang, Y. (2007). On fast computation of the non-parametric maximum
  likelihood estimate of a mixing distribution. \emph{Journal of the
  Royal Statistical Society, Ser. B}, \bold{69}, 185-198.

  Wang, Y. (2010). Maximum likelihood computation for fitting
  semiparametric mixture models. \emph{Statistics and Computing},
  \bold{20}, 75-86

  Wang, Y. and Taylor, S. M. (2013). Efficient computation of
  nonparametric survival functions via a hierarchical mixture
  formulation. \emph{Statistics and Computing}, \bold{23}, 713-725.

}

\seealso{ \code{\link{nnls}}, \code{\link{npnorm}},
  \code{\link{nppois}}, \code{\link{cnmms}}. }

\examples{
## Simulated data
x = rnppois(1000, c(1,4), c(0.7,0.3)) # Poisson mixture
(r = cnm(x))
plot(r, x)

x = rnpnorm(1000, mu=c(0,4), pr=c(0.3,0.7), sd=1) # Normal mixture
plot(cnm(x), x)                        # sd = 1
plot(cnm(x, init=list(beta=0.5)), x)   # sd = 0.5
mix0 = disc(seq(min(x$v),max(x$v), len=100)) # over a finite grid
plot(cnm(x, init=list(beta=0.5, mix=mix0), model="p"),
     x, add=TRUE, col="blue")          # An approximate NPMLE

## Real-world data
data(thai)
plot(cnm(x <- nppois(thai)), x)     # Poisson mixture

data(brca)
plot(cnm(x <- npnorm(brca)), x)     # Normal mixture

}
\keyword{ function }   % at least one, from doc/KEYWORDS
