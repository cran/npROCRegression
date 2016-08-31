create.design.matrices <- function(iROCf, names.cov, mode, data, data.f) {
	np <- ncol(iROCf$II)
	dm <- list()
	npar <- 0
	nparl <- NULL
	ncoeff <- 0
	for(i in 1:np){
		if(iROCf$II[1,i] == -1) {
			if(mode[iROCf$II[2,i]] == 6) {
				npar <- npar + 1
				levels <- levels(data[ ,names.cov[iROCf$II[2,i]], drop=TRUE])
				levels <- levels[na.omit(unique(data.f[ ,names.cov[iROCf$II[2,i]], drop=TRUE]))]
				dm[[npar]] <- contr.sum(levels)
				attr(dm[[npar]], "varnames") <- paste(names.cov[iROCf$II[2,i]],"_", levels, sep = "")
				ncoeff <- ncoeff + length(levels) - 1
				nparl <- c(nparl, length(levels) - 1)
			}
			if (mode[iROCf$II[2,i]] == 5 & iROCf$h[i] == 0) {
				npar <- npar + 1
				dm[[npar]] <- matrix(1, ncol = 1, nrow = 1)
				attr(dm[[npar]], "varnames") <- names.cov[iROCf$II[2,i]]
				ncoeff <- ncoeff + 1
				nparl <- c(nparl, 1)
			}
		}			 
	}
	res <- list(dm = dm, ncoeff = ncoeff, nparl = nparl)
	res
}