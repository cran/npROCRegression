interpret.ROCformula <-
function(formula,data) {
	env <- environment(formula) 
	if(inherits(formula, "character"))		  
		formula <- as.formula(formula)
	tf <- terms.formula(formula, specials = c("s"))
	terms <- attr(tf, "term.labels")
	if(length(grep(":",terms))!=0)  stop("Symbol '*' is not allowed")
	
	nt <- length(terms)	 
	ns <- attr(tf,"specials")$s
	II <- list()
	h  <- list()
	partial=vector()
	k <- 0
	if(nt) {
		for (i in 1:nt) {
			if (i %in% ns) {
				k=k+1				   
				st <- eval(parse(text = paste("ROC.",terms[i],sep="")))
				II[[k]] <- match(st$cov,names(data),-1)
				h[[k]] <- -1
				partial[k] <- terms[i]														  
			} else {
				pos <- match(terms[i],names(data))				  
				if(!is.na(pos)) {
					k=k+1
					II[[k]]<- c(-1,pos)
					h[[k]] <- 0
					partial[k] <- terms[i]
				}
			}
		}		   
	}	   
	II <- if(length(II)) {
		matrix(unlist(II), nrow=2)
	} else {
		matrix(0,nrow=2)
	}	   
	res<-list(II=II, h=unlist(h), npartial=k, partial = partial)
	res
}
