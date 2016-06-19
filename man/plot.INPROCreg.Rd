\name{plot.INPROCreg}
\alias{plot.INPROCreg}
\title{
Default INPROCreg plotting
}
\description{
Default INPROCreg plotting (see details)
}
\usage{
\method{plot}{INPROCreg}(x, ask = TRUE, ...)
}
\arguments{
  \item{x}{an object of class \code{INPROCreg} as produced by \code{INPROCreg()}.}
  \item{ask}{a logical value. If TRUE, the default, the user is asked for confirmation, before a new figure is drawn.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
The function produces the following plots:
\describe{
	\item{(a)}{the estimated regression and variance functions in both the healthy and diseased populations.
	}
	\item{(b)}{the covariate-specific ROC curve and AUC.
	}
	\item{(c)}{the covariate-adjusted ROC curve (AROC).
	}
	\item{(d)}{(optionally) the Youden Index (YI) or the value for which the TPF and the TNF coincides (EQ).
	}
	\item{(e)}{(optionally) the optimal thresholds based on these criteria (TH).
	}
}
}
%\value{}
%\references{}
\author{
Maria Xose Rodriguez-Alvarez and Javier Roca-Pardinas
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
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
