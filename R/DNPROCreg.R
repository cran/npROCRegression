DNPROCreg <-
function(marker, formula.h=~1, formula.ROC=~1, group, tag.healthy, data, ci.fit=FALSE, test.partial = NULL, newdata=NULL, control = controlDNPROCreg(), weights=NULL) {
	control <- do.call("controlDNPROCreg", control)	
	if(inherits(formula.h, "formula")) formula.h <- c(formula.h)
	if(length(formula.h) == 1) formula.h[[2]] <- formula.h[[1]]
	if(inherits(formula.h, "character")) {
		m <- list()
		m[[1]] <- as.formula(formula.h[[1]])
		m[[2]] <- as.formula(formula.h[[2]])
		formula.h = m
	}
	if(inherits(formula.ROC, "formula")) formula.ROC <- c(formula.ROC)  
	if(inherits(formula.ROC, "character")) formula.ROC <- c(as.formula(formula.ROC))
		
	names.cov.hm <- all.vars(formula.h[[1]])
	names.cov.hv <- all.vars(formula.h[[2]])	
	names.cov.ROC <- all.vars(formula.ROC[[1]])		 
	names.cov <- c(names.cov.hm, names.cov.hv[is.na(match(names.cov.hv,c(names.cov.hm,names.cov.ROC)))], names.cov.ROC[is.na(match(names.cov.ROC,c(names.cov.hm,names.cov.hv)))])		   
	
	
	if(!is.null(newdata) && !inherits(newdata, "data.frame"))
		stop("Newdata must be a data frame")
	if(sum(is.na(match(c(marker,names.cov,group), names(data)))))
		stop("Not all needed variables are supplied in data")	
	if(sum(!is.null(newdata) && is.na(match(names.cov, names(newdata)))))
		stop("Not all needed variables are supplied in newdata")		 
	if(length(unique(data[,group]))!=2)
		stop(paste(group," variable must have only two different values (for healthy and diseased individuals)"), sep="")
		
	data.new <- na.omit(data[,c(marker,group,names.cov)])
	data.new[,group]<- data.new[,group]!=tag.healthy
	
	data.h <- data.new[data.new[,group]==0,]
	data.d <- data.new[data.new[,group]!=0,]
		
	n <- nrow(data.new)
	n0<- nrow(data.h)
	n1<- nrow(data.d)
	
	mode <- lapply(names.cov, function(x,data) switch(class(data[,x,drop=TRUE]),"integer" = 5, "numeric" = 5,"factor" = 6,"character" = 6), data=data.new)
				
	extract.fhm <- interpret.ROCformula(formula.h[[1]], data.new[,names.cov,drop=FALSE])	
	extract.fhv <- interpret.ROCformula(formula.h[[2]], data.new[,names.cov,drop=FALSE])			
	extract.fROC<- interpret.ROCformula(formula.ROC[[1]], data.new[,names.cov,drop=FALSE])
	
	if(any((mode[extract.fROC$II[2,]])[test.partial] == 6))
		stop("Testing for the effect of a factor variable not yet implemented")	 
	
	if(is.null(newdata)) {
		newdata <- DNPROCregData(data.new, names.cov, group)
	} else {
		newdata <- na.omit(newdata[,names.cov,drop=FALSE])
	}
	# Data frames to Fortran
	data.f <- data.new
	newdata.f <- newdata
	names.cat <- names.cov[mode==6]
	if (length(names.cat))  
		for (i in 1:length(names.cat)){		 
			levels <- levels(data.new[,names.cat[i], drop=TRUE])				
			data.f[,names.cat[i]] <- match(data.new[,names.cat[i], drop=TRUE],levels)
			newdata.f[,names.cat[i]] <- match(newdata[,names.cat[i], drop=TRUE],levels)
		}
	
	set.t <- seq(0, 1, by = control$step.p)
	l.set.t <- length(set.t)
	l.set.cont <- nrow(newdata)
	
	ROC <- rep(0.0,l.set.t*l.set.cont)	
	AUC <- rep(0.0,l.set.cont*3)
	pfunctions <- rep(0.0, l.set.t*l.set.cont*(extract.fROC$npartial+1)*3)
	pvalue <- rep(-1, 2*length(test.partial))
	
	# Create design matrices for the coefficients
	dm <- create.design.matrices(extract.fROC, names.cov, mode, data.new, data.f)
	coeff <- rep(-1, dm$ncoeff + 1)
	
	if(is.null(weights))
		weights=rep(1.0,n)
		
	fit=.Fortran("DLLROCDirect",									  
				 Z=matrix(as.double(as.matrix(data.f[,names.cov])),ncol=length(names.cov)),
				 X=as.double(data.f[,marker]),
				 W=as.double(weights),			   
				 Status = as.integer(data.f[,group]),
				 n=as.integer(c(n,n0,n1)),			   
				 tagh=as.integer(0),
				 nt = as.integer(control$card.P),				
				 nvar = as.integer(length(names.cov)),
				 mode = as.integer(mode),
				 nparm = as.integer(extract.fhm$npartial),			   
				 IIm = matrix(as.integer(extract.fhm$II),nrow=2),
				 hm = as.double(extract.fhm$h),
				 nparv = as.integer(extract.fhv$npartial),			   
				 IIv = matrix(as.integer(extract.fhv$II),nrow=2),
				 hv = as.double(extract.fhv$h),
				 nparr = as.integer(extract.fROC$npartial),			  
				 IIr = matrix(as.integer(extract.fROC$II),nrow=2),
				 hr = as.double(extract.fROC$h),
				 family = as.integer(ifelse(control$link=="probit",7,1)),
				 Zb=matrix(as.double(as.matrix(newdata.f[,names.cov])),ncol=length(names.cov)),
				 nb=as.integer(l.set.cont),
				 ntb=as.integer(l.set.t),
				 p=as.integer(control$p),				   
				 kbin=as.integer(control$kbin),
				 coutcome= is.element("coutcome",control$resample.m),
				 nboot=as.integer(control$nboot),
				 seed=as.integer(control$seed),
				 cifit=as.logical(ci.fit),
				 level=as.double(control$level),
				 tpartial = as.integer(test.partial),
				 npartial = length(test.partial),
				 pfunctions=array(as.double(pfunctions),dim=c(l.set.t*l.set.cont, extract.fROC$npartial+1,3)),
				 coeff=as.double(coeff),	
				 ROC=as.double(ROC),
				 AUC=as.double(AUC),
				 pvalue=matrix(as.double(as.matrix(pvalue)), nrow = 2), PACKAGE = "npROCRegression")
		
	columns  <- switch(as.character(ci.fit),"TRUE" = 1:3, "FALSE" = 1)					 
	col.names <- c("","ll","ul")
	
	# Partial functions 
	pftemp <- array(fit$pfunctions, dim = c(l.set.t*l.set.cont, extract.fROC$npartial+1, 3))
	colnames.cov <- NULL
	m <- matrix(ncol = (extract.fROC$npartial+1)*length(columns), nrow = l.set.t*l.set.cont)	
	for (i in 1:extract.fROC$npartial) {		
		colnames.cov <- c(colnames.cov, paste(extract.fROC$partial[i], col.names,sep="")[columns])
		m[,((i-1)*(length(columns))+1):(i*(length(columns)))] <- matrix(pftemp[,i,], ncol = 3)[, columns, drop=FALSE]
	}
	m[,((extract.fROC$npartial)*(length(columns))+1):((extract.fROC$npartial+1)*(length(columns)))] <- matrix(pftemp[,extract.fROC$npartial+1,], ncol=3)[,columns,drop=FALSE]   
	m <- data.frame(m)
	colnames.fpf <- paste("s(fpf)",col.names, sep="")[columns]
	colnames(m)<- c(colnames.cov, colnames.fpf)
	
	# Return two "matrices": one for the covariates the other for the FPF
	p.functions.cov <- m[seq(1,l.set.t*l.set.cont, by = l.set.t), colnames.cov, drop = FALSE]
	p.functions.fpf <- m[1:l.set.t, colnames.fpf, drop = FALSE]
	
	# Named the coefficients
	e <- cumsum(dm$nparl)
	s <- e - dm$nparl + 1
	ncoeff <- fit$coeff[-1]
	nncoeff <- fit$coeff[1]
	names(nncoeff) <- "(Intercept)"
	if(!is.null(dm$nparl) && dm$nparl != 0) 
		for(i in 1:length(dm$nparl)) {
			aux <- dm$dm[[i]]%*%ncoeff[s[i]:e[i]]
			names(aux) <- attr(dm$dm[[i]], "varnames")
			nncoeff <- c(nncoeff, aux)
		}  
	
	res <- list()
	res$call  <- match.call()
	res$marker <- marker
	res$group <- group
	res$formula.h <- formula.h
	res$formula.ROC <- formula.ROC	
	res$ci.fit <- ci.fit	
	res$model <- data.new						 
	res$fpf <- set.t
	res$newdata <- newdata
	res$pfunctions <- list(covariates = p.functions.cov, fpf = p.functions.fpf)
	res$coefficients <- nncoeff
	res$ROC <- t(array(fit$ROC, dim=c(l.set.t,l.set.cont), dimnames=list(set.t,1:l.set.cont)))
	res$AUC <- array(fit$AUC, dim=c(l.set.cont,3),dimnames=list(1:l.set.cont, paste("AUC",col.names,sep="")))[,columns,drop=FALSE]		  
	res$pvalue <- matrix(round(fit$pvalue,4), nrow = 2, dimnames = list(c("T2","T1"), extract.fROC$partial[test.partial]))
	class(res)<-"DNPROCreg"
	res
}
