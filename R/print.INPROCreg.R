print.INPROCreg <-
function(x,...) {
	cat("\nCall:\n") 
	print(x$call)    
	cat("\n*************************************************\n")
	cat("Induced non-parametric ROC regression")
	cat("\n*************************************************\n")
	cat(paste("MARKER: ",x$marker,"\n"))
	cat(paste("COVARIATE: ",x$covariate,"\n"))
	cat(paste("STATUS: ",x$group,"\n"))
	invisible(x)   
}
