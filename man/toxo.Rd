\name{toxo}
\alias{toxo}

\docType{data}

\title{
  Toxoplasmosis Data
}

\description{

  Contains the number of subjects testing positively for toxoplasmosis
  in 34 cities of El Salvador, with various rainfalls.

}

\usage{toxo}

\format{ A numeric matrix with four columns:

  city:     city identification code.

  y:        the number of subjects testing positively for toxoplasmosis.
  
  n:        the number of subjects tested.
  
  rainfall: the annual rainfall of the city, in meters. 
  
}

\references{

  Efron, B. (1986). Double exponential families and their use in
  generalized linear regression. \emph{Journal of the American
  Statistical Association}, \bold{81}, 709-721.

  Aitkin, M. (1996). A general maximum likelihood analysis of
  overdispersion in generalised linear models. \emph{Statistics and
  Computing}, \bold{6}, 251-262.

  Wang, Y. (2010). Maximum likelihood computation for fitting
  semiparametric mixture models. \emph{Statistics and Computing},
  \bold{20}, 75-86.

}

\seealso{ \code{\link{mlogit}},\code{\link{cnmms}}. }

\examples{

data(toxo)
x = mlogit(toxo)
cnmms(x)

}

\keyword{datasets}

