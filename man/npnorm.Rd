\name{npnorm}
\alias{npnorm}
\alias{rnpnorm}
\alias{plot.npnorm}

%- Also NEED an '\alias' for EACH other topic documented here.

\title{Class `npnorm'}

\description{

  Class \code{npnorm} is used to store data that will be processed as
  those of a nonparametric normal mixture.
  
  Function \code{npnorm} creates an object of class \code{npnorm}, given
  values and weights/frequencies.
  
  Function \code{rnpnorm} generates a random sample from a normal
  mixture and saves the data as an object of class \code{npnorm}.

  Function \code{plot.npnorm} plots the normal mixture.

}

\usage{

npnorm(v, w=1)
rnpnorm(n, mu=0, pr=1, sd=1)
\method{plot}{npnorm}(x, mix, beta, breaks=NULL, col="red", add=FALSE,
                       border.col=NULL, border.lwd=1, fill="lightgrey",
                       components=TRUE, main, lwd=2,
                       lty=1, xlab="Data", ylab="Density", ...)

}

\arguments{
  
  \item{v}{a numeric vector that stores the values of a sample.}
  
  \item{w}{a numeric vector that stores the corresponding
    weights/frequencies of the observations.}
  
  \item{n}{the sample size.}
  
  \item{mu}{a numeric vector for all the component means.}
  
  \item{pr}{a numeric vector for all the mixing proportions.}
  
  \item{sd}{a scalar for the component standard deviation that is common
    to all components.}
  
  \item{x}{an object of class \code{npnorm}.}
  
  \item{mix}{an object of class \code{disc}.}
  
  \item{beta}{the structural parameter.}
  
  \item{breaks}{the rough number bins used for plotting the histogram.}
  
  \item{col}{the color of the density curve to be plotted.}
  
  \item{add}{if \code{FALSE}, creates a new plot; if \code{TRUE}, adds
    the plot to the existing one.}
  
  \item{border.col}{color for the border of histogram boxes.}
	  
  \item{border.lwd}{line width for the border of histogram boxes.}
  
  \item{fill}{color to fill in the histogram boxes.}
	  
  \item{components}{if \code{TRUE}, also show the support points and
    mixing proportions.}
  
  \item{main, lwd, lty, xlab, ylab}{arguments for graphical
    parameters (see \code{par}).}
  
  \item{...}{arguments passed on to function \code{plot}.}

}

\details{
  
  When \code{components=TRUE}, the component means are shown on the
  horizontal line of density 0. The vertical lines going upwardly at the
  support points are proportional to the mixing proportions at these
  points.

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
mix = disc(pt=c(0,4), pr=c(0.3,0.7))
x = rnpnorm(200, mix$pt, mix$pr, sd=1)
plot(x, mix, beta=1)

}
\keyword{ class }   % at least one, from doc/KEYWORDS
\keyword{ function }   % at least one, from doc/KEYWORDS
