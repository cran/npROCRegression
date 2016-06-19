\name{DNPROCregData}
\alias{DNPROCregData}
\title{
Selects an adequate set of points from a data set for obtaining predictions or plots.
}
\description{
Selects an adequate set of points from a data set to be used as a default dataset for obtaining predictions or plots.
}
\usage{
DNPROCregData(data, names.cov, group)
}
\arguments{
  \item{data}{Data set from which the new set of covariate values is obtained.
}
  \item{names.cov}{Character vector with the names of the covariates to be included in the new data set.
}
  \item{group}{A character string with the name of the variable in the original data set that distinguishes healthy from diseased individuals.
}
}
%\details{}
\value{a data frame containing selected values of all needed covariates. For those that are continuous, 30 different values are selected.}
%\references{}
\author{
Maria Xose Rodriguez - Alvarez and Javier Roca-Pardinas
}
%\note{}
\seealso{
See Also \code{\link{DNPROCreg}}.
}
\examples{
data(endosim)
# Fit a model including the interaction between age and gender.
m0 <- DNPROCreg(marker = "bmi", formula.h = "~ gender + s(age) + s(age, by = gender)", 
				formula.ROC = "~ gender + s(age) + s(age, by = gender)", 
				group = "idf_status", 
				tag.healthy = 0, 
				data = endosim, 
				control = list(card.P=50, kbin=30, step.p=0.02), 
				ci.fit = FALSE, 
				test.partial = NULL,
				newdata = NULL)
summary(m0)				
plot(m0)			 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
