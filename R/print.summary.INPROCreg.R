print.summary.INPROCreg <-
function(x,...) {
	print.INPROCreg(x)
	if(x$pvalue != -1) {
		cat("\n----------------------------------------------")
		cat(paste("\nTest for continuous covariate effect (p-value): ",round(x$pvalue,4), "\n ",sep=""))
	}
	invisible(x)
}
