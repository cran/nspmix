\name{brca}
\alias{brca}

\docType{data}

\title{
  Z-values of BRCA Data
}

\description{

  Contains 3226 \eqn{z}-values computed by Efron (2004) from the data
  obtained in a well-known microarray experiment concerning two types of
  genetic mutations causing increased breast cancer risk, BRCA1 and
  BRCA2.

}

\usage{brca}

\format{ A numeric vector containing 3226 \eqn{z}-values. 
  
}

\references{

  Efron, B. (2004). Large-scale simultaneous hypothesis testing: the
  choice of a null hypothesis. \emph{Journal of the American Statistical
  Association}, \bold{99}, 96-104.

  Wang, Y. (2007). On fast computation of the non-parametric maximum
  likelihood estimate of a mixing distribution. \emph{Journal of the
  Royal Statistical Society, Ser. B}, \bold{69}, 185-198.

  Wang, Y. and C.-S. Chee (2012). Density estimation using nonparametric
  and semiparametric mixtures. \emph{Statistical Modelling: An
  International Journal}, \bold{12}, 67-92.
  
}

\seealso{ \code{\link{npnorm}},\code{\link{cnm}}. }

\examples{

data(brca)
x = npnorm(brca)
plot(cnm(x), x)

}

\keyword{datasets}

