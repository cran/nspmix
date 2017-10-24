\name{plotgrad}
\alias{plotgrad}

%- Also NEED an '\alias' for EACH other topic documented here.

\title{Plot the Gradient Function}
\description{
  
  Function \code{plotgrad} plots the gradient function or its first
  derivative of a nonparametric mixture.
  
}

\usage{

plotgrad(x, mix, beta, len=500, order=0, col="blue", col2="red",
    add=FALSE, main=paste0("Class: ",class(x)), xlab=expression(theta),
    ylab=paste0("Gradient (order = ",order,")"), cex=1, pch=1, lwd=1,
    xlim, ylim, ...)

}

\arguments{
  
  \item{x}{a data object of a mixture model class.}
  
  \item{mix}{an object of class 'disc', for a discrete mixing distribution.}
  
  \item{beta}{the structural parameter.} 
  
  \item{len}{number of points used to plot the smooth curve.}
  
  \item{order}{the order of the derivative of the gradient function to
    be plotted. If 0, it is the gradient function itself.}
  
  \item{col}{color for the curve.}
  
  \item{col2}{color for the support points.}
  
  \item{add}{if \code{FALSE}, create a new plot; if \code{TRUE}, add
    the curve and points to the current one.}
  
  \item{main, xlab, ylab, cex, pch, lwd, xlim, ylim}{arguments for graphical
    parameters (see \code{par}).}
  
  \item{...}{arguments passed on to function \code{plot}.}

}

\details{
  
  \code{data} must belong to a mixture family, as specified by its
  class.

  The support points are shown on the horizontal line of gradient 0. The
  vertical lines going downwards at the support points are proportional
  to the mixing proportions at these points.

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

\seealso{ \code{\link{plot.nspmix}}, \code{\link{nnls}}, \code{\link{cnm}},
  \code{\link{cnmms}}, \code{\link{npnorm}}, \code{\link{nppois}}. }

\examples{
## Poisson mixture
x = rnppois(200, c(1,4), c(0.7,0.3))
r = cnm(x)
plotgrad(x, r$mix)

## Normal mixture
x = rnpnorm(200, mu=c(0,4), pr=c(0.3,0.7), sd=1)
r = cnm(x, init=list(beta=0.5))   # sd = 0.5
plotgrad(x, r$mix, r$beta)

}
\keyword{ function }   % at least one, from doc/KEYWORDS
