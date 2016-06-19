print.DNPROCreg <-
function(x,...) {
	cat("\nCall:\n") 
	print(x$call)	
	cat("\n*************************************************\n")
	cat("Direct non-parametric ROC regression")
	cat("\n*************************************************\n")
	cat(paste("Marker: ",x$marker,"\n"))
	cat(paste("Group (status): ",x$group,"\n"))
	cat(paste("Healthy regression model (mean): ", deparse(x$formula.h[[1]]),"\n"))
	cat(paste("Healthy regression model (variance): ", deparse(x$formula.h[[2]]),"\n"))
	cat(paste("ROC regression model: ",x$formula.ROC,"\n"))		
	invisible(x)   
}
