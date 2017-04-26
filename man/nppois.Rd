\name{nppois}
\alias{nppois}
\alias{rnppois}
\alias{plot.nppois}

%- Also NEED an '\alias' for EACH other topic documented here.

\title{Class `nppois'}

\description{

  Class \code{nppois} is used to store data that will be processed as
  those of a nonparametric Poisson mixture.
  
  Function \code{nppois} creates an object of class \code{nppois}, given
  values and weights/frequencies.
  
  Function \code{rnppois} generates a random sample from a Poisson mixture and
  saves the data as an object of class \code{nppois}.

  Function \code{plot.nppois} plots the Poisson mixture.

}

\usage{

nppois(v, w=1)
rnppois(n, lambda=1, pr=1)
\method{plot}{nppois}(x, mix, beta, col="red", add=FALSE,
                      components=TRUE, main="Poisson Mixture", lwd=1,
                      lty=1, xlab="Data", ylab="Density", ...)

}

\arguments{
  
  \item{v}{a numeric vector that stores the values of a sample.}
  
  \item{w}{a numeric vector that stores the corresponding
           weights/frequencies of the observations.}
  
  \item{n}{the sample size.}
  
  \item{lambda}{a numeric vector for all the component means.}
  
  \item{pr}{a numeric vector for all the mixing proportions.}
  
  \item{x}{an object of class \code{nppois}.}
  
  \item{mix}{an object of class \code{disc}.}
  
  \item{beta}{the structural parameter, which is not really needed for
              the Poisson mixture.}
  
  \item{col}{the color of the density curve to be plotted.}
  
  \item{add}{if \code{FALSE}, creates a new plot; if \code{TRUE}, adds
            the plot to the existing one.}
  
  \item{components}{if \code{TRUE}, also show the support points and
        mixing proportions.}
  
  \item{main, lwd, lty, xlab, ylab}{arguments for graphical
        parameters (see \code{par}).}
  
  \item{...}{arguments passed on to function \code{plot}.}

}

\details{
  
  When \code{components=TRUE}, the support points are shown on the
  horizontal line of density 0. The component density curves, weighted
  appropriately, are also shown.

}

\author{ Yong Wang <yongwang@auckland.ac.nz>}

\references{

  Wang, Y. (2007). On fast computation of the non-parametric maximum
  likelihood estimate of a mixing distribution. \emph{Journal of the
  Royal Statistical Society, Ser. B}, \bold{69}, 185-198.

}

\seealso{ \code{\link{nnls}}, \code{\link{cnm}}, \code{\link{cnmms}},
  \code{\link{plot.nspmix}}. }

\examples{
mix = disc(pt=c(1,4), pr=c(0.3,0.7))
x = rnppois(200, mix$pt, mix$pr)
plot(x, mix)

}
\keyword{ class }   % at least one, from doc/KEYWORDS
\keyword{ function }   % at least one, from doc/KEYWORDS
