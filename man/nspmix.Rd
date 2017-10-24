\name{nspmix}
\alias{nspmix}
\alias{plot.nspmix}

%- Also NEED an '\alias' for EACH other topic documented here.

\title{Class 'nspmix'}
\description{

  Class \code{nspmix} is an object returned by function \code{cnm},
  \code{cnmms}, \code{cnmpl} or \code{cnmap}.
  
  Function \code{plot.nspmix} plots either the mixture model, if the
  family of the mixture provides an implementation of the generic
  \code{plot} function, or the gradient function.
  
}

\usage{

\method{plot}{nspmix}(x, data, type=c("probability","gradient"), ...)

}

\arguments{
  
  \item{x}{an object of a mixture model class}
  
  \item{data}{a data set from the mixture model} 
  
  \item{type}{the type of function to be plotted: the probability model of
    the mixture family (\code{probability}), or the gradient function
    (\code{gradient}). }
  
  \item{...}{arguments passed on to the \code{plot} function called.}

}

\details{
  
  \code{data} must belong to a mixture family, as specified by its
  class.

}

\author{ Yong Wang <yongwang@auckland.ac.nz>}

\references{

  Wang, Y. (2007). On fast computation of the non-parametric maximum
  likelihood estimate of a mixing distribution. \emph{Journal of the
  Royal Statistical Society, Ser. B}, \bold{69}, 185-198.

  Wang, Y. (2010). Maximum likelihood computation for fitting
  semiparametric mixture models. \emph{Statistics and Computing},
  \bold{20}, 75-86

}

\seealso{ \code{\link{nnls}}, \code{\link{cnm}}, \code{\link{cnmms}},
  \code{\link{cnmpl}}, \code{\link{cnmap}}, \code{\link{npnorm}},
  \code{\link{nppois}}. }

\examples{
## Poisson mixture
x = rnppois(200, c(1,4), c(0.7,0.3))
r = cnm(x)
plot(r, x, "p")
plot(r, x, "g")

## Normal mixture
x = rnpnorm(200, mu=c(0,4), pr=c(0.3,0.7), sd=1)
r = cnm(x, init=list(beta=0.5))   # sd = 0.5
plot(r, x, "p")
plot(r, x, "g")

}
\keyword{ function }   % at least one, from doc/KEYWORDS
