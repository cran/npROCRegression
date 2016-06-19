\name{DNPROCreg}
\alias{DNPROCreg}
\title{
Direct nonparametric ROC regression modelling
}
\description{
Estimates the covariate-specific ROC curve in the presence of multidimensional covariates by means of the ROC-GAM regression model presented in Rodriguez- Alvarez et al. (2011)
}
\usage{
DNPROCreg(marker, formula.h = ~1, formula.ROC = ~1, group, tag.healthy, data, 
ci.fit = FALSE, test.partial = NULL, newdata = NULL, control = controlDNPROCreg(), 
weights = NULL)
}
\arguments{
  \item{marker}{A character string with the name of the diagnostic test variable.
}
  \item{formula.h}{Right-hand formula(s) giving the regression and variance model(s) to be fitted in healthy population. Atomic values are also valid, being recycled.
}
  \item{formula.ROC}{Right-hand formula giving the ROC regression to be fitted.
}
  \item{group}{A character string with the name of the variable that distinguishes healthy from diseased individuals.
}
  \item{tag.healthy}{The value codifying the healthy individuals in the variable \code{group}.
}
  \item{data}{Data frame representing the data and containing all needed variables.
}
  \item{ci.fit}{A logical value. If TRUE, confidence intervals are computed.
}
  \item{test.partial}{A numeric vector containing the covariate terms in the ROC regression model to be tested for a possible effect. If NULL, no test is performed.
}
  \item{newdata}{A data frame containing the values of the covariate at which predictions are required.
}
  \item{control}{Output of the \code{controlDNROCreg()} function.
}
  \item{weights}{An optional vector of `prior weights' to be used in the fitting process.
}
}
%\details{}
\value{
As a result, the function \code{DNPROCreg()} provides a list with the following components:
	\item{call}{The matched call.}
	\item{model}{Data frame containing all variables and observations used in the fitting process.}
	\item{fpf}{Set of false positive fractions (FPF) at which the covariate-specific ROC curve has been estimated.}
	\item{newdata}{Data frame containing the values of the covariates at which estimates has been obtained.}
	\item{pfunctions}{Matrices containing the estimates of each component of the additive predictor of the ROC regression model. One matrix contains the effects of the covariates, the other the effect of the FPF. Confidence intervals are returned if required).}	
	\item{coefficients}{Vector of parametric coefficient of the fitted ROC regression model.}
	\item{ROC}{Estimated covariate-specific ROC curve.}
	\item{AUC}{Estimated covariate-specific AUC, and corresponding confidence intervals if required.}
	\item{pvalue}{If required, p-values are obtained  - with two different bootstrap-based tests - for each model term indicated in argument \code{test.partial} (T2: \eqn{L_{2}}-based test; and T1: \eqn{L_{1}}-based test). See Rodriguez-Alvarez et al. (2016).}
}
\references{

Rodriguez- Alvarez, M.X., Roca-Pardinas, J. and Cadarso-Suarez, C. (2011). A new flexible direct ROC regression model - Application to the detection of cardiovascular risk factors by anthropometric measures. Computational Statistics and Data Analysis, 55(12), 3257--3270.

Rodriguez- Alvarez, M.X., Roca-Pardinas, J. and Cadarso-Suarez, C. (2016). Boostrap-based procedures for inference in nonparametric ROC regression analysis. Technical report.
}
\author{
Maria Xose Rodriguez-Alvarez and Javier Roca-Pardinas
}
%\note{}
\seealso{
See Also as \code{\link{INPROCreg}}, \code{\link{summary.DNPROCreg}}, \code{\link{plot.DNPROCreg}}, \code{\link{controlDNPROCreg}}, \code{\link{DNPROCregData}}.
}
\examples{
data(endosim)
# Fit a model including the interaction between age and gender.
m0 <- DNPROCreg(marker = "bmi", formula.h = "~ gender + s(age) + s(age, by = gender)", 
				formula.ROC = "~ gender + s(age) + s(age, by = gender)", 
				group = "idf_status", 
				tag.healthy = 0, 
				data = endosim, 
				control = list(card.P=50, kbin=30, step.p=0.02))
summary(m0)				
plot(m0)

\dontrun{
# For confidence intervals
set.seed(123)
m1 <- DNPROCreg(marker = "bmi", formula.h = "~ gender + s(age) + s(age, by = gender)", 
				formula.ROC = "~ gender + s(age) + s(age, by = gender)", 
				group = "idf_status", 
				tag.healthy = 0, 
				data = endosim, 
				control = list(card.P=50, kbin=30, step.p=0.02),
				ci.fit = TRUE)
summary(m1)
plot(m1)

# For testing the presence of interaccion between age and gender
set.seed(123)
m2 <- DNPROCreg(marker = "bmi", formula.h = "~ gender + s(age) + s(age, by = gender)", 
				formula.ROC = "~ gender + s(age) + s(age, by = gender)", 
				group = "idf_status", 
				tag.healthy = 0, 
				data = endosim, 
				control = list(card.P=50, kbin=30, step.p=0.02),
				test.partial = 3)
summary(m2)
plot(m2)			
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
