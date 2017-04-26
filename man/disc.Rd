\name{disc}
\alias{disc}
\alias{print.disc}

%- Also NEED an '\alias' for EACH other topic documented here.

\title{Class `disc'}

\description{

  Class \code{disc} is used to represent an arbitrary univariate
  discrete distribution with a finite number of support points.
  
  Function \code{disc} creates an object of class \code{disc}, given the
  support points and probability values at these points.
  
  Function \code{print.disc} prints the discrete distribution.

}

\usage{

disc(pt, pr=1)
\method{print}{disc}(x, ...)

}

\arguments{
  
  \item{pt}{a numeric vector for support points.}
  
  \item{pr}{a numeric vector for probability values at the support points.}
  
  \item{x}{an object of class \code{disc}.}
  
  \item{...}{arguments passed on to function \code{print}.}

}

\author{ Yong Wang <yongwang@auckland.ac.nz>}

\seealso{ \code{\link{cnm}}, \code{\link{cnmms}}. }

\examples{

(d = disc(pt=c(0,4), pr=c(0.3,0.7)))

}
\keyword{ class }   % at least one, from doc/KEYWORDS
\keyword{ function }   % at least one, from doc/KEYWORDS
