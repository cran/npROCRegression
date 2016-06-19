plot.INPROCreg <-
function(x, ask = TRUE, ...) {
	plot.accuracy <- function(x, accuracy, dots, ask, ci.fit=FALSE) {
		# Order the covariate value
		if(accuracy != "AROC") {
			ord <- order(x$X)
			x$X <- x$X[ord]
			x[[accuracy]] <- x[[accuracy]][ord,,drop = FALSE]
		}
		if(ask)
			readline("Press return for next page....")
		plot(0, 0, xlab = ifelse(accuracy=="AROC","FPF",x$covariate), ylab = accuracy, xlim = if(accuracy=="AROC") c(0,1) else range(x$X), ylim = if(accuracy == "TH") range(x[[accuracy]]) else c(0,1), type = "n")
		lines(if(accuracy=="AROC") x$fpf else x$X, x[[accuracy]][,1])
		if(ci.fit) {
			lines(if(accuracy=="AROC") x$fpf else x$X, x[[accuracy]][,2], lty=2)
			lines(if(accuracy=="AROC") x$fpf else x$X, x[[accuracy]][,3], lty=2)
		}			  	
		if (accuracy == "AUC") abline(h = 0.5, col = "grey")
		if (accuracy == "AROC") abline(0,1, col = "grey")
	}
	ci.fit <- ifelse(is.null(x$ci.fit),FALSE,x$ci.fit)
	# Means
	range.x <- c(min(x$X),max(x$X))
	range.y <- c(min(min(x$h$M0),min(x$d$M1)),max(max(x$h$M0),max(x$d$M1)))
	ord <- order(x$X)   
	plot(x$X[ord], x$d$M1[ord,1], xlim=range.x, ylim=range.y, type="l", lwd=2,main="Mean", xlab=x$covariate, ylab=x$marker)
	lines(x$X[ord], x$h$M0[ord,1], lwd=2, lty=2)
	if(ci.fit) {
		lines(x$X[ord],x$d$M1[ord,2], lwd=1)
		lines(x$X[ord],x$d$M1[ord,3], lwd=1)
		lines(x$X[ord],x$h$M0[ord,2], lty=2)
		lines(x$X[ord],x$h$M0[ord,3], lty=2)
	}
	legend("bottomright", legend=c("Diseased","Healthy"), lty=1:2, lwd=2)
	if(ask)
			readline("Press return for next page....")   
	# Variances	
	range.y <- sqrt(c(max(min(min(x$h$V0),min(x$d$V1)),0),max(max(x$h$V0),max(x$d$V1))))  
	plot(x$X[ord],sqrt(x$d$V1[ord,1]), xlim=range.x, ylim=range.y, type="l", lwd=2, main="Standard deviation", xlab=x$covariate, ylab=x$marker)
	lines(x$X[ord],sqrt(x$h$V0[ord,1]), lwd=2, lty=2)
	if(ci.fit) {
		lines(x$X[ord],sqrt(x$d$V1[ord,2]), lwd=1)
		lines(x$X[ord],sqrt(x$d$V1[ord,3]), lwd=1)
		lines(x$X[ord],sqrt(x$h$V0[ord,2]), lty=2)
		lines(x$X[ord],sqrt(x$h$V0[ord,3]), lty=2)
	}
	legend("bottomright", legend=c("Diseased","Healthy"), lty=1:2, lwd=2)
	
	#### Accuracy measures	
	dots <- list(...)
	set.accuracy <- c("EQ", "YI", "TH")
	set.plot=c("AUC","AROC",set.accuracy[is.element(set.accuracy, names(x))])
		
	if(ask)
		readline("Press return for next page....")
	delete.obs <- duplicated(x$X)	
	X.tmp <- x$X[!delete.obs]
	ROC.tmp <- x$ROC[!delete.obs,,drop = FALSE]
	ord.tmp <- order(X.tmp) 	
	persp(x$fpf, X.tmp[ord.tmp], t(ROC.tmp[ord.tmp,,drop = FALSE]),
	xlab = "FPF", ylab = x$covariate, zlab = "TPF",	
  	theta = if (!is.null(dots$theta))dots$theta else 20,
	phi = if (!is.null(dots$phi))dots$phi else 30,
	col = if(!is.null(dots$col))dots$col else "white",
	shade = if(!is.null(dots$shade))dots$shade else 0.5, ticktype = "detailed",
	cex.axis = dots$cex.axis, cex.lab = dots$cex.lab, cex.sub = dots$cex.sub,cex = dots$cex)

	for(i in 1:length(set.plot))
		plot.accuracy(x, set.plot[i], dots, ask, ci.fit) 
}
