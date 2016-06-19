\name{INPROCreg}
\alias{INPROCreg}
\title{
Induced nonparametric ROC regression modelling
}
\description{
Estimates the covariate-specific ROC curve (and related measures) in the presence of a one-dimensional continuous covariate based on the induced nonparametric ROC regression approach as presented in Rodriguez-Alvarez et al. (2011).  
}
\usage{
INPROCreg(marker, covariate, group, tag.healthy, data, ci.fit = FALSE, 
test = FALSE, accuracy = NULL, accuracy.cal = c("ROC", "AROC"), 
newdata = NULL, control = controlINPROCreg(), weights = NULL)
}

\arguments{
  \item{marker}{A character string with the name of the diagnostic test variable.
}
  \item{covariate}{A character string with the name of the continuous covariate.
}
  \item{group}{A character string with the name of the variable that distinguishes healthy from diseased individuals.
}
  \item{tag.healthy}{The value codifying the healthy individuals in the variable \code{group}.
}
  \item{data}{Data frame representing the data and containing all needed variables.
}
  \item{ci.fit}{A logical value. If TRUE, confidence intervals are computed.
}
  \item{test}{A logical value. If TRUE, the bootstrap-based test for detecting covariate effect is performed.
}
  \item{accuracy}{A character vector indicating if the Youden index (``YI''), the value for which the TPF and the TNF coincides (``EQ''), and/or optimal threshold (``TH'') based on these two criteria should be computed. 
}
  \item{accuracy.cal}{A character string indicating if the accuracy measures should be calculated based on the covariate-specific ROC curve or on the covariate-adjusted ROC curve (AROC).
}
  \item{newdata}{A data frame containing the values of the covariate at which predictions are required.
}
  \item{control}{Output of the \code{controlINROCreg()} function.
}
  \item{weights}{An optional vector of `prior weights' to be used in the fitting process.
}
}
%\details{}
\value{
As a result, the function \code{INPROCreg()} provides a list with the following components:
	\item{call}{The matched call.}
	\item{X}{The data frame used in the predictions.}
	\item{fpf}{Set of false positive fractions at which the covariate-specific ROC curve has been estimated.}
	\item{h}{Estimated regression and variance functions in healthy population.} 
	\item{d}{Estimated regression and variance functions in diseased population.}
	\item{ROC}{Estimated covariate-specific ROC curve.}	
	\item{AUC}{Estimated covariate-specific AUC, and corresponding confidence intervals if required.}
	\item{AROC}{Estimated covariate-adjusted ROC curve.}
	\item{YI/EQ}{If required, estimated covariate-specific YI (or values at which the TPF and the TNF coincide), and corresponding bootstrap confidence intervals.}
	\item{TH}{If required, estimated optimal threshold values based on either the YI or the criterion of equality of TPF and TNF, and corresponding bootstrap confidence intervals.}
	\item{pvalue}{If required, p-value obtained with the test for checking the effect of the continuous covariate on the ROC curve.}
}
\references{

Gonzalez - Manteiga, W., Pardo-Fernandez, J.C. and van Keilegom, I. (2011). ROC curves in nonparametric location-scale regression models. Scandinavian Journal of Statistics, 38, 169--184.

Rodriguez - Alvarez, M.X., Roca-Pardinas, J. and Cadarso-Suarez, C. (2011). ROC curve and covariates: extending induced methodology to the non-parametric framework. Statistics and Computing, 21(4), 483--499.

Yao, F., Craiu, R.V. and Reiser, B. (2010). Nonparametric covariate adjustment for receiver operating characteristic curves. The Canadian Journal of Statistics, 38, 27--46.
}
\author{
Maria Xose Rodriguez - Alvarez and Javier Roca-Pardinas
}
%\note{}
\seealso{
See Also as \code{\link{DNPROCreg}}, \code{\link{summary.INPROCreg}}, \code{\link{plot.INPROCreg}}, \code{\link{controlINPROCreg}}.
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
# Evaluate the effect of age on the accuracy of the body mass index for females
m0.women <- INPROCreg(marker = "bmi", covariate = "age", group = "idf_status", 
						tag.healthy = 0, 
						data = subset(endosim, gender == "Women"), 
						ci.fit = FALSE, test = FALSE, 
						accuracy = c("EQ","TH"),
						accuracy.cal="ROC", 
						control=controlINPROCreg(p=1,kbin=30,step.p=0.01), 
						newdata = data.frame(age = seq(18,85,l=50)))
						
summary(m0.women)						
plot(m0.women)
\dontrun{
# For computing confidence intervals and testing covariate effect
set.seed(123)
m1.men <- INPROCreg(marker = "bmi", covariate = "age", group = "idf_status", 
						tag.healthy = 0, 
						data = subset(endosim, gender == "Men"), 
						ci.fit = TRUE, test = TRUE, 
						accuracy = c("EQ","TH"),
						accuracy.cal="AROC", 
						control=controlINPROCreg(p=1,kbin=30,step.p=0.01), 
						newdata = data.frame(age = seq(18,85,l=50)))
summary(m1.men)
plot(m1.men)					
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
