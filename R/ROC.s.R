ROC.s <-
function(x1=NULL,x2=NULL, by=NULL) {   
	args <- match.call()
	if(!is.null(args$x1) & !is.null(args$x2) & !is.null(args$by)) 
		stop("Factor interaction terms (by argument) are not allowed in surface interactions")
	if(is.null(args$x1) & is.null(args$x2))
		stop("x1 or x2 must be indicated")
	if(!is.null(args$x1) & is.null(args$x2) & is.null(args$by)) {			   
		cov = c("", deparse(args$x1, backtick = TRUE, width.cutoff = 500))
	} else if (!is.null(args$x1) & !is.null(args$x2) & is.null(args$by)) {	  
		cov = c(deparse(args$x1, backtick = TRUE, width.cutoff = 500),deparse(args$x2, backtick = TRUE, width.cutoff = 500))		
	} else if (!is.null(args$x1) & is.null(args$x2) & !is.null(args$by)) {	  
		cov = c(deparse(args$by, backtick = TRUE, width.cutoff = 500), deparse(args$x1, backtick = TRUE, width.cutoff = 500))
	} else {
		stop("Invalid expression")
	}
	res <- list(cov=cov)
	res
}
