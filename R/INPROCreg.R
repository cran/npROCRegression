INPROCreg <-
function(marker, covariate, group, tag.healthy, data, ci.fit=FALSE, test=FALSE, accuracy=NULL, accuracy.cal = c("ROC","AROC"), newdata=NULL, control = controlINPROCreg(), weights=NULL) {
	control <- do.call("controlINPROCreg", control)    
	accuracy.cal = match.arg(accuracy.cal)
	set.accuracy <- c("EQ","YI", "TH")
	if(!is.null(newdata) && !inherits(newdata, "data.frame"))
		stop("Newdata must be a data frame")   
	if(sum(is.element(c("EQ","YI"), accuracy))==2)
		stop("Threshold values must be calculated by the Youden Index (YI) criterion or Se=Sp (EQ) criterion")   
	if(sum(is.na(match(c(marker,covariate,group), names(data)))))
		stop("Not all needed variables are supplied in data")    
	if(sum(!is.null(newdata) && is.na(match(covariate, names(newdata)))))
		stop("Not all needed variables are supplied in newdata")
	if(length(unique(data[,group]))!=2)
		stop(paste(group," variable must have only two different values (for healthy and diseased individuals)"), sep="")
	if(is.null(newdata)) {
		newdata <- unique(na.omit(data[order(data[,covariate]),covariate]))
	} else {
		newdata <- newdata[,covariate]
	}
	data <-	na.omit(data[,c(marker,covariate,group)])
	data[,group]<- data[,group]!=tag.healthy
	
	data.h <- data[data[,group]==0,]
	data.d <- data[data[,group]!=0,]

	n <- nrow(data)
	n0<- nrow(data.h)
	n1<- nrow(data.d)

	set.t <- seq(0, 1, by = control$step.p)
	l.set.t <- length(set.t)
	l.set.cont <- length(newdata)
	accuracy.call <- c(is.element("ROC",accuracy.cal), is.element("YI",accuracy), is.element("YI",accuracy) | is.element("EQ",accuracy), is.element("TH",accuracy))
    
	ROC<-rep(0.0,l.set.t*l.set.cont)
	AROC<-rep(0.0,l.set.t*3)
	AUC<-YI<-TH<-M0<-V0<-M1<-V1<-rep(0.0,l.set.cont*3)
	pvalor = -1.0
	if(is.null(weights))
		weights=rep(1.0,n)
		
	fit=.Fortran("DLLROCInduced",					 				  
				Z=as.double(data[,covariate]), 				      
				X=as.double(data[,marker]),
				W=as.double(weights), 
				Status = as.integer(data[,group]),
				n=as.integer(c(n,n0,n1)),
				tagh=as.integer(0),
				Zb=as.double(newdata), 
				nb=as.integer(l.set.cont),
				nt=as.integer(l.set.t),
				p=as.integer(control$p),
				h=matrix(as.double(as.matrix(control$h, ncol = 2)), ncol = 2),
				kbin=as.integer(control$kbin),
				coutcome= is.element("coutcome",control$resample.m),
				nboot=as.integer(control$nboot),
				seed=as.integer(control$seed),
				cifit=as.logical(ci.fit),
				level=as.double(control$level),   
				test=as.logical(test),
				accuracy=as.logical(accuracy.call),
				M0=as.double(M0),
				V0=as.double(V0),
				M1=as.double(M1),
				V1=as.double(V1),
				ROC=as.double(ROC),
				AUC=as.double(AUC),
				AROC=as.double(AROC),
				YI=as.double(YI),
				TH=as.double(TH),
				pvalor=as.double(pvalor), PACKAGE = "npROCRegression")
		
	columns <-switch(as.character(ci.fit),"TRUE"=1:3,"FALSE"=1)			    		
	col.names<-c("fit","fitll","fitul")
	
	res <- list()
	res$call <- match.call()
	res$marker <- marker
	res$covariate <- covariate
	res$group <- group
	res$ci.fit <- ci.fit
	res$X <- matrix(newdata, ncol=1, dimnames=list(1:l.set.cont,covariate))
	res$fpf <- set.t
	res$h <- list(M0=array(fit$M0, dim=c(l.set.cont,3), dimnames=list(1:l.set.cont,col.names))[,columns,drop=FALSE] , V0=array(fit$V0, dim=c(l.set.cont,3), dimnames=list(1:l.set.cont,col.names))[,columns,drop=FALSE])	
	res$d <- list(M1=array(fit$M1, dim=c(l.set.cont,3), dimnames=list(1:l.set.cont,col.names))[,columns,drop=FALSE] , V1=array(fit$V1, dim=c(l.set.cont,3), dimnames=list(1:l.set.cont,col.names))[,columns,drop=FALSE])
	res$ROC <- t(array(fit$ROC, dim=c(l.set.t,l.set.cont), dimnames=list(set.t,1:l.set.cont)))			
	res$AUC <- array(fit$AUC, dim=c(l.set.cont,3),dimnames=list(1:l.set.cont,col.names))[,columns,drop=FALSE]		
	res$AROC <- array(fit$AROC, dim=c(l.set.t,3), dimnames=list(set.t,col.names))[,columns,drop=FALSE]	
	if(is.element("YI",accuracy)|is.element("EQ",accuracy))	 		
		res[[ifelse(is.element("YI",accuracy),"YI","EQ")]] <- array(fit$YI, dim=c(l.set.cont,3), dimnames=list(1:l.set.cont,col.names))[,columns,drop=FALSE] 
	if(is.element("TH",accuracy))
		res$TH <- array(fit$TH, dim=c(l.set.cont,3), dimnames=list(1:l.set.cont,col.names))[,columns,drop=FALSE]
	res$pvalue <- fit$pvalor
	class(res)<-"INPROCreg"
	res
}
