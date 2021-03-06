% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/disc.R
\name{disc}
\alias{disc}
\alias{print.disc}
\title{Class `disc'}
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
\description{
Class \code{disc} is used to represent an arbitrary univariate discrete
distribution with a finite number of support points.
}
\details{
Function \code{disc} creates an object of class \code{disc}, given the
support points and probability values at these points.

Function \code{print.disc} prints the discrete distribution.
}
\examples{


(d = disc(pt=c(0,4), pr=c(0.3,0.7)))

}
\seealso{
\code{\link{cnm}}, \code{\link{cnmms}}.
}
\author{
Yong Wang <yongwang@auckland.ac.nz>
}
\keyword{class}
\keyword{function}
