\name{cvps}
\alias{cvp}
\alias{rcvp}
\alias{cvps}
\alias{rcvps}
\alias{print.cvps}

\title{Class `cvps'}

\description{

  These functions can be used to study a common variance problem (CVP),
  where univariate observations fall in known groups. Observations in
  each group are assumed to have the same mean, but different groups may
  have different means. All observations are assumed to have a common
  variance, despite their different means, hence giving the name of the
  problem. It is a random-effects problem. 

  Class \code{cvps} is used to store the CVP data in a summarized form.
  
  Function \code{cvps} creates an object of class \code{cvps}, given a
  matrix that stores the values (column 2) and their grouping
  information (column 1).
  
  Function \code{rcvp} generates a random sample in the raw form for a
  common variance problem, where the means follow a discrete
  distribution.

  Function \code{rcvps} generates a random sample in the summarized form
  for a common variance problem, where the means follow a discrete
  distribution.

  Function \code{print.cvps} prints the CVP data given in the summarized
  form.

}

\usage{

cvps(x)
rcvp(k, ni=2, mu=0, pr=1, sd=1)
rcvps(k, ni=2, mu=0, pr=1, sd=1)
\method{print}{cvps}(x, ...)

}

\arguments{
  
  \item{x}{CVP data in the raw form as an argument in \code{cvps}, or an
    object of class \code{cvps} in \code{print.cvps}.}
  
  \item{k}{the number of groups.}
  
  \item{ni}{a numeric vector that gives the sample size in each group.}
  
  \item{mu}{a numeric vector for all the theoretical means.}
  
  \item{pr}{a numeric vector for all the probabilities associated with
            the theoretical means.}
  
  \item{sd}{a scalar for the standard deviation that is common to all
            observations.}
  
  \item{...}{arguments passed on to function \code{print}.}

}

\details{
  
  The raw form of the CVP data is a two-column matrix, where each row
  represents an observation. The two columns along each row give,
  respectively, the group membership (\code{group}) and the value
  (\code{x}) of an observation.

  The summarized form of the CVP data is a four-column matrix, where
  each row represents the summarized data for all observations in a
  group. The four columns along each row give, respectively, the group
  number (\code{group}), the number of observations in the group
  (\code{ni}), the sample mean of the observations in the group
  (\code{mi}), and the residual sum of squares of the observations in
  the group (\code{ri}).

}

\author{ Yong Wang <yongwang@auckland.ac.nz>}

\references{

  Neyman, J. and Scott, E. L. (1948). Consistent estimates based on
  partially consistent observations. \emph{Econometrica}, \bold{16},
  1-32.

  Kiefer, J. and Wolfowitz, J. (1956). Consistency of the maximum
  likelihood estimator in the presence of infinitely many incidental
  parameters. \emph{Ann. Math. Stat.}, \bold{27}, 886-906.

  Wang, Y. (2010). Maximum likelihood computation for fitting
  semiparametric mixture models. \emph{Statistics and Computing},
  \bold{20}, 75-86.

}

\seealso{ \code{\link{nnls}}, \code{\link{cnmms}}. }

\examples{

mix = disc(pt=c(0,4), pr=c(0.3,0.7))
x = rcvp(100, ni=3:10, mu=mix$pt, pr=mix$pr, sd=1)
cnmms(cvps(x))

}
\keyword{ class }   % at least one, from doc/KEYWORDS
\keyword{ function }   % at least one, from doc/KEYWORDS
