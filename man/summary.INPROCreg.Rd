\name{summary.INPROCreg}
\alias{summary.INPROCreg}
\title{
Summary method for INPROCreg objects.
}
\description{
Summary method for INPROCreg objects.
}
\usage{
\method{summary}{INPROCreg}(object, ...)
}
\arguments{
  \item{object}{an object of class \code{INPROCreg} as produced by \code{INPROCreg()}}.
  \item{\dots}{further arguments passed to or from other methods. Not yet implemented.}
}
%\details{}
%\value{}
%\references{}
\author{
Maria Xose Rodriguez - Alvarez and Javier Roca-Pardinas
}
%\note{}
\seealso{
See Also \code{\link{INPROCreg}}.
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
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }
