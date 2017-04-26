\name{betablockers}
\alias{betablockers}

\docType{data}

\title{
  Beta-blockers Data
}

\description{

  Contains the data of the 22-center clinical trial of beta-blockers for
  reducing mortality after myocardial infarction.

}

\usage{betablockers}

\format{ A numeric matrix with four columns:

  center:   center identification code.

  deaths:   the number of deaths in the center.
  
  total:    the number of patients taking beta-blockers in the center.
  
  treatment: 0 for control, and 1 for treatment.
  
}

\source{
  
  Aitkin, M. (1999). A general maximum likelihood analysis of variance
  components in generalized linear models. \emph{Biometrics}, \bold{55},
  117-128.

}

\references{

  Wang, Y. (2010). Maximum likelihood computation for fitting
  semiparametric mixture models. \emph{Statistics and Computing},
  \bold{20}, 75-86.

}

\seealso{ \code{\link{mlogit}},\code{\link{cnmms}}. }

\examples{

data(betablockers)
x = mlogit(betablockers)
cnmms(x)

}

\keyword{datasets}

