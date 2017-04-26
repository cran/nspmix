\name{lungcancer}
\alias{lungcancer}

\docType{data}

\title{
  Lung Cancer Data
}

\description{

  Contains the data of 14 studies of the effect of smoking on lung cancer.

}

\usage{lungcancer}

\format{ A numeric matrix with four columns:

  study:    study identification code.

  lungcancer:   the number of people diagnosed with lung cancer.
  
  size:      the number of people in the study.
  
  smoker: 0 for smoker, and 1 for non-smoker.
  
}

\source{
  
  Booth, J. G. and Hobert, J. P. (1999). Maximizing generalized linear
  mixed model likelihoods with an automated Monte Carlo EM
  algorithm. \emph{Journal of the Royal Statistical Society, Ser. B},
  \bold{61}, 265-285.

}

\references{

  Wang, Y. (2010). Maximum likelihood computation for fitting
  semiparametric mixture models. \emph{Statistics and Computing},
  \bold{20}, 75-86.

}

\seealso{ \code{\link{mlogit}},\code{\link{cnmms}}. }

\examples{

data(lungcancer)
x = mlogit(lungcancer)
cnmms(x)

}

\keyword{datasets}

