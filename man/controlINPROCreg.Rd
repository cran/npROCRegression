\name{controlINPROCreg}
\alias{controlINPROCreg}
\title{
Function used to set several parameters controlling fitting process.
}
\description{
Function used to set several parameters controlling fitting process.
}
\usage{
controlINPROCreg(step.p = 0.02, kbin = 30, p = 1, h = c(-1, -1, -1, -1), 
seed = NULL, nboot = 500, level = 0.95, resample.m = c("coutcome", "ncoutcome"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{step.p}{a numeric value, defaulting to 0.02. ROC curves are calculated at a regular sequence of false positive fractions with \code{step.p} increment.}
  \item{kbin}{an integer value specifying the number of binning knots. By default 30.}
  \item{p}{an integer value specifying the order of the local polynomial kernel estimator for the regression functions. By default 1.}
  \item{h}{a vector of length 4 specifying the bandwidths to be used for the estimation of the regression and variance functions in healthy population and the regression and variance functions in diseased populations (in this order). By default -1 (selected using cross-validation). A value of 0 would indicate a linear fit.}
  \item{seed}{an integer value specifying the seed for the bootstrap resamples. If NULL it is initialized randomly.}
  \item{nboot}{an integer value specifying the number of bootstrap resamples for the construction of the confidence intervals. By default 500.}
  \item{level}{a real value specifying the confidence level for the confidence intervals. By default 0.95.}
  \item{resample.m}{a character string specifying if bootstrap resampling (for the confidence intervals) should be done with or without regard to the disease status (``coutcome'' or ``noutcome''). When the resampling method is done conditionally on the disease status, the resampling is based on the residuals of the regression models in healthy and diseased populations. However, when the bootstrap resampling is done without regard to the disease status, a naive bootstrap is used. By default, the resampling is done conditionally on the disease status.}
}
%\details{}
%\value{}
%\references{}
\author{Maria Xose Rodriguez - Alvarez and Javier Roca-Pardinas
}
%\note{}
\seealso{
See Also \code{\link{INPROCreg}}
}
\examples{
data(endosim)
# Evaluate the effect of age on the accuracy of the body mass index for males
m0.men <- INPROCreg(marker = "bmi", covariate = "age", group = "idf_status", 
						tag.healthy = 0, 
						data = subset(endosim, gender == "Men"), 
						ci.fit = FALSE, test = FALSE, 
						accuracy = c("EQ","TH"),
						accuracy.cal="AROC", 
						control=controlINPROCreg(p=1,kbin=30,step.p=0.01), 
						newdata = data.frame(age = seq(18,85,l=50)))	

summary(m0.men)
plot(m0.men)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
