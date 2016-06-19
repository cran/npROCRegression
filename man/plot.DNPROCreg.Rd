\name{plot.DNPROCreg}
\alias{plot.DNPROCreg}
\title{
Default DNPROCreg plotting
}
\description{
Takes a fitted \code{DNPROCreg} object produced by \code{DNPROCreg()} and plots the covariate-specific ROC curve and associated AUC.
}
\usage{
\method{plot}{DNPROCreg}(x, ask = TRUE, ...)
}
\arguments{
  \item{x}{an object of class \code{DNPROCreg} as produced by \code{DNPROCreg()}}
  \item{ask}{a logical value. If TRUE, the default, the user is asked for confirmation, before a new figure is drawn}
  \item{\dots}{further arguments passed to or from other methods.}
}
%\details{}
%\value{}
%\references{}
\author{
Maria Xose Rodriguez-Alvarez and Javier Roca-Pardinas
}
%\note{}
\seealso{
See Also \code{\link{DNPROCreg}}, \code{\link{summary.DNPROCreg}}.
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
}
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
