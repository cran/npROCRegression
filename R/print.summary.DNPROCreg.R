print.summary.DNPROCreg <-
function(x,...) {
	print.DNPROCreg(x)
	cat("\n----------------------------------------------\n")
	cat("Parametric coefficients (ROC curve):\n")
		print(x$coefficients)		 
	if(ncol(x$pvalue)!=0) {
		cat("\n----------------------------------------------\n")	   
		cat("Tests for effects (p-values)\n")
		print(x$pvalue)		 
	}
	invisible(x)	
}
