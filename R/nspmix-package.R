##'Beta-blockers Data
##'
##'
##'Contains the data of the 22-center clinical trial of beta-blockers for
##'reducing mortality after myocardial infarction.
##'
##'
##'@name betablockers
##'@docType data
##'@format A numeric matrix with four columns:
##'
##'center: center identification code.
##'
##'deaths: the number of deaths in the center.
##'
##'total: the number of patients taking beta-blockers in the center.
##'
##'treatment: 0 for control, and 1 for treatment.
##'@seealso \code{\link{mlogit}},\code{\link{cnmms}}.
##'@references
##'
##'Wang, Y. (2010). Maximum likelihood computation for fitting semiparametric
##'mixture models. \emph{Statistics and Computing}, \bold{20}, 75-86.
##'@source
##'
##'Aitkin, M. (1999). A general maximum likelihood analysis of variance
##'components in generalized linear models. \emph{Biometrics}, \bold{55},
##'117-128.
##'@keywords datasets
##'@examples
##'
##'
##'data(betablockers)
##'x = mlogit(betablockers)
##'cnmms(x)
##'
##'
NULL





##'Z-values of BRCA Data
##'
##'
##'Contains 3226 \eqn{z}-values computed by Efron (2004) from the data obtained
##'in a well-known microarray experiment concerning two types of genetic
##'mutations causing increased breast cancer risk, BRCA1 and BRCA2.
##'
##'
##'@name brca
##'@docType data
##'@format A numeric vector containing 3226 \eqn{z}-values.
##'@seealso \code{\link{npnorm}},\code{\link{cnm}}.
##'@references
##'
##'Efron, B. (2004). Large-scale simultaneous hypothesis testing: the choice of
##'a null hypothesis. \emph{Journal of the American Statistical Association},
##'\bold{99}, 96-104.
##'
##'Wang, Y. (2007). On fast computation of the non-parametric maximum
##'likelihood estimate of a mixing distribution. \emph{Journal of the Royal
##'Statistical Society, Ser. B}, \bold{69}, 185-198.
##'
##'Wang, Y. and C.-S. Chee (2012). Density estimation using nonparametric and
##'semiparametric mixtures. \emph{Statistical Modelling: An International
##'Journal}, \bold{12}, 67-92.
##'@keywords datasets
##'@examples
##'
##'
##'data(brca)
##'x = npnorm(brca)
##'plot(cnm(x), x)
##'
##'
NULL





##'Lung Cancer Data
##'
##'
##'Contains the data of 14 studies of the effect of smoking on lung cancer.
##'
##'
##'@name lungcancer
##'@docType data
##'@format A numeric matrix with four columns:
##'
##'study: study identification code.
##'
##'lungcancer: the number of people diagnosed with lung cancer.
##'
##'size: the number of people in the study.
##'
##'smoker: 0 for smoker, and 1 for non-smoker.
##'@seealso \code{\link{mlogit}},\code{\link{cnmms}}.
##'@references
##'
##'Wang, Y. (2010). Maximum likelihood computation for fitting semiparametric
##'mixture models. \emph{Statistics and Computing}, \bold{20}, 75-86.
##'@source
##'
##'Booth, J. G. and Hobert, J. P. (1999). Maximizing generalized linear mixed
##'model likelihoods with an automated Monte Carlo EM algorithm. \emph{Journal
##'of the Royal Statistical Society, Ser. B}, \bold{61}, 265-285.
##'@keywords datasets
##'@examples
##'
##'
##'data(lungcancer)
##'x = mlogit(lungcancer)
##'cnmms(x)
##'
##'
NULL





##'Class 'nspmix'
##'
##'
##'Class \code{nspmix} is an object returned by function \code{cnm},
##'\code{cnmms}, \code{cnmpl} or \code{cnmap}.
##'
##'Function \code{plot.nspmix} plots either the mixture model, if the family of
##'the mixture provides an implementation of the generic \code{plot} function,
##'or the gradient function.
##'
##'
##'\code{data} must belong to a mixture family, as specified by its class.
##'
##'@name plot.nspmix
##'@aliases nspmix plot.nspmix
##'@param x an object of a mixture model class
##'@param data a data set from the mixture model
##'@param type the type of function to be plotted: the probability model of the
##'mixture family (\code{probability}), or the gradient function
##'(\code{gradient}).
##'@param ... arguments passed on to the \code{plot} function called.
##'@author Yong Wang <yongwang@@auckland.ac.nz>
##'@seealso \code{\link{nnls}}, \code{\link{cnm}}, \code{\link{cnmms}},
##'\code{\link{cnmpl}}, \code{\link{cnmap}}, \code{\link{npnorm}},
##'\code{\link{nppois}}.
##'@references
##'
##'Wang, Y. (2007). On fast computation of the non-parametric maximum
##'likelihood estimate of a mixing distribution. \emph{Journal of the Royal
##'Statistical Society, Ser. B}, \bold{69}, 185-198.
##'
##'Wang, Y. (2010). Maximum likelihood computation for fitting semiparametric
##'mixture models. \emph{Statistics and Computing}, \bold{20}, 75-86
##'@keywords function
##'@examples
##'
##'## Poisson mixture
##'x = rnppois(200, disc(c(1,4), c(0.7,0.3)))
##'r = cnm(x)
##'plot(r, x, "p")
##'plot(r, x, "g")
##'
##'## Normal mixture
##'x = rnpnorm(200, mix=disc(c(0,4), c(0.3,0.7)), sd=1)
##'r = cnm(x, init=list(beta=0.5))   # sd = 0.5
##'plot(r, x, "p")
##'plot(r, x, "g")
##'
##'@usage
##'\method{plot}{nspmix}(x, data, type=c("probability","gradient"), ...)
##'
##'@export plot.nspmix
NULL





##'Illness Spells and Frequencies of Thai Preschool Children
##'
##'
##'Contains the results of a cohort study in north-east Thailand in which 602
##'preschool children participated. For each child, the number of illness
##'spells \eqn{x}, such as fever, cough or running nose, is recorded for all
##'2-week periods from June 1982 to September 1985. The frequency for each
##'value of \eqn{x} is saved in the data set.
##'
##'
##'@name thai
##'@docType data
##'@format A data frame with 24 rows and 2 variables:
##'
##'x: values of \eqn{x}.
##'
##'freq: frequencies for each value of \eqn{x}.
##'@seealso \code{\link{nppois}},\code{\link{cnm}}.
##'@references
##'
##'Wang, Y. (2007). On fast computation of the non-parametric maximum
##'likelihood estimate of a mixing distribution. \emph{Journal of the Royal
##'Statistical Society, Ser. B}, \bold{69}, 185-198.
##'@source
##'
##'Bohning, D. (2000). \emph{Computer-assisted Analysis of Mixtures and
##'Applications: Meta-analysis, Disease Mapping, and Others}. Boca Raton:
##'Chapman and Hall-CRC.
##'@keywords datasets
##'@examples
##'
##'
##'data(thai)
##'x = nppois(thai)
##'plot(cnm(x), x)
##'
##'
NULL





##'Toxoplasmosis Data
##'
##'
##'Contains the number of subjects testing positively for toxoplasmosis in 34
##'cities of El Salvador, with various rainfalls.
##'
##'
##'@name toxo
##'@docType data
##'@format A numeric matrix with four columns:
##'
##'city: city identification code.
##'
##'y: the number of subjects testing positively for toxoplasmosis.
##'
##'n: the number of subjects tested.
##'
##'rainfall: the annual rainfall of the city, in meters.
##'@seealso \code{\link{mlogit}},\code{\link{cnmms}}.
##'@references
##'
##'Efron, B. (1986). Double exponential families and their use in generalized
##'linear regression. \emph{Journal of the American Statistical Association},
##'\bold{81}, 709-721.
##'
##'Aitkin, M. (1996). A general maximum likelihood analysis of overdispersion
##'in generalised linear models. \emph{Statistics and Computing}, \bold{6},
##'251-262.
##'
##'Wang, Y. (2010). Maximum likelihood computation for fitting semiparametric
##'mixture models. \emph{Statistics and Computing}, \bold{20}, 75-86.
##'@keywords datasets
##'@examples
##'
##'
##'data(toxo)
##'x = mlogit(toxo)
##'cnmms(x)
##'
##'
NULL



