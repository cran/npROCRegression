\name{controlDNPROCreg}
\alias{controlDNPROCreg}
\title{
Function used to set several parameters controlling the ROC regression fitting process
}
\description{
Function used to set several parameters controlling the ROC regression fitting process
}
\usage{
controlDNPROCreg(step.p = 0.02, card.P = 50, link = c("probit", "logit"), 
kbin = 30, p = 1, seed = NULL, nboot = 500, level = 0.95, 
resample.m = c("coutcome", "ncoutcome"))
}
\arguments{
  \item{step.p}{a numeric value, defaulting to 0.02. ROC curves are calculated at a regular sequence of false positive fractions with \code{step.p} increment.
}
  \item{card.P}{an integer value specifying the cardinality of the set of false positive fractions used in the estimation processs. By default 50.
}
  \item{link}{a character string specifying the link function (``probit'' or ``logit''). By default the link is the probit function.
}
  \item{kbin}{an integer value specifying the number of binning knots. By default 30.
}
  \item{p}{an integer value specifying the order of the local polinomial kernel estimator. By default 1.
}
  \item{seed}{an integer value specifying the seed for the bootstrap resamples. If NULL it is initialized randomly.
}
  \item{nboot}{an integer value specifying the number of bootstrap resamples for the construction of the confidence intervals. By default 500.   
}
  \item{level}{a real value specifying the confidence level for the confidence intervals. By default 0.95  
}
  \item{resample.m}{a character string specifying if bootstrap resampling (for the confidence intervals) should be done with or without regard to the disease status (``coutcome'' or ``noutcome''). In both cases, a naive bootstrap is used. By default, the resampling is done conditionally on the disease status.
}
}
%\details{}
%\value{}
%\references{}
\author{
Maria Xose Rodriguez - Alvarez and Javier Roca-Pardinas
}
%\note{}
\seealso{
See Also \code{\link{DNPROCreg}}
}
\examples{
data(endosim)
# Fit a model including the interaction between age and gender.
m0 <- DNPROCreg(marker = "bmi", formula.h = "~ gender + s(age) + s(age, by = gender)", 
				formula.ROC = "~ gender + s(age) + s(age, by = gender)", 
				group = "idf_status", 
				tag.healthy = 0, 
				data = endosim, 
				control = controlDNPROCreg(card.P=50, kbin=30, step.p=0.02))
summary(m0)				
plot(m0)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
