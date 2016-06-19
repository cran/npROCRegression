\name{print.DNPROCreg}
\alias{print.DNPROCreg}
\title{
Print method for DNPROCreg objects
}
\description{
Print method for DNPROCreg objects
}
\usage{
\method{print}{DNPROCreg}(x, ...)
}
\arguments{
  \item{x}{an object of class \code{DNPROCreg} as produced by \code{DNPROCreg()}}
  \item{\dots}{further arguments passed to or from other methods. Not yet implemented}
}
%\details{}
%\value{}
%\references{}
\author{
Maria Xose Rodriguez-Alvarez and Javier Roca-Pardinas
}
%\note{}
\seealso{
See Also as \code{\link{DNPROCreg}}.
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
m0				
summary(m0)				
plot(m0)
}
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
