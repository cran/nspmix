\name{whist}
\alias{whist}

\title{Weighted Histograms}
\description{
  
  Plots or computes the histogram with observations with
  multiplicities/weights.
  
}

\usage{

whist(x, w=1, breaks="Sturges", plot=TRUE, freq=NULL, xlim=NULL,
      ylim=NULL, xlab="Data", ylab=NULL, main=NULL, add=FALSE, col=NULL,
      border=NULL, lwd=1, ...) 

}

\arguments{
  
  \item{x}{a vector of values for which the histogram is desired.}
  
  \item{w}{a vector of multiplicities/weights for the values in
    \code{x}.}
  
  \item{breaks,plot,freq,xlim,ylim,xlab,ylab,main,add,col,border,lwd}{These
    arguments have similar functionalities to their namesakes in
    function \code{hist}.}

  \item{...}{arguments passed on to function \code{plot}.}

}

\details{
  
  Just like \code{hist}, \code{whist} can either plot the histogram or
  compute the values that define the histogram, by setting \code{plot}
  to \code{TRUE} or \code{FALSE}. 

  The histogram can either be the one for frequencies or density, by
  setting \code{freq} to \code{TRUE} or \code{FALSE}.

}

\author{ Yong Wang <yongwang@auckland.ac.nz>}

\seealso{ \code{\link{hist}}. }

\keyword{ function }   % at least one, from doc/KEYWORDS
