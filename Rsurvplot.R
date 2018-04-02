# Survial plot of selected prognostic features
# Usage: Rsurvplot(file, clinical, lassotable)
# file: p by n matrix of data, where p is the number
#       of features, n is the number of samples; 
#       rownames are feature names; colnames are sample names
# clinical: n by q matrix of clinical data, where n is 
#           the number of samples, q is the clinical informtion,
#           including at least survival time and survival event
# lasso: r by 2 lasso coeficients matrix, r is the number of 
#         selected prognostic feature values,the first column is 
#         the feature names, the second column is the 
#         coefficient of each feature 
# Created by Yidi Sun (sunyidi@picb.ac.cn)

# load required packages
load_packages=function(){
  a = function(){
    if ( require (survival))
      return()
    else
      install.packages("survival")
    }
    
  a()
}
load_packages()

# LASSO regression
require(survival)

Rsurvplot = function(file, clinical, lassotable){
	data=read.table(file, header=TRUE, row.names=1)
	clin=read.table(clinical, header=TRUE, row.names=1)
	lasso=read.table(lassotable, header=TRUE)
	data=t(data)
	dbnew=data[,colnames(data) %in% lasso[,1]]
	coef = as.numeric(as.character(lasso[,2]))
	lassovalue = apply(dbnew,1,function(x){sum(x*coef)})
	lassovalue = as.data.frame(lassovalue)
	lassovalue[lassovalue$lassovalue>median(lassovalue$lassovalue),"class"] = 1
	lassovalue[lassovalue$lassovalue<=median(lassovalue$lassovalue),"class"] = 0
	survlasso = survfit(Surv(clin$surv_year,clin$vital_status)~lassovalue$class)
	pdf("survplot.pdf",width=5,height=5)
	survplot(survlasso,group.names=c("Low risk","High risk"),extra.left.margin=4,
  			loc.legend='bottomleft',col.surv = c("blue","red"),mark="",
  			xlab="Time(Years)",ylab="Overall Survival Probability",las=1)
	dev.off()
}

survplot <- function(survfit, mark=3, simple=FALSE,
	xaxis.at=pretty(survfit$time), xaxis.lab=xaxis.at, 
	lty.surv=1, lwd.surv=1, col.surv=1,
	lty.ci=0, lwd.ci=1, col.ci=col.surv, #By default (lty.ci=0), confidence intervals are not plotted.
	group.names=NULL, group.order=seq(length(survfit$n)), extra.left.margin=4, 
	label.n.at.risk=TRUE, draw.lines=FALSE, cex.axis=1,
	xlab='', ylab='', main='', xlim=c(0,max(survfit$time)), ylim=c(0,1),
	grid=FALSE, lty.grid=1, lwd.grid=1, col.grid=grey(.9),
	legend=!is.null(survfit$strata), loc.legend='topright', add=FALSE,
	... # ... is passed to par()
	) {

	# xaxis.at specifies where 'n at risk' will be computed and printed.
	# xaxis.lab specifies what will be printed at xaixs.at.  (see example)

	# If group names are long, add extra left margin by setting extra.left.margin to something greater than 0.

	# line specifications (lty.surv, lwd.surv, col.surv) will be recycled.
	# Set lty.ci to 1 if confidence intervals are needed.
	# group.names will overwrite whatever is specified in survfit() output.
	# group.order specifies the order of groups from top in 'n at risk'.  1 is at top, 2 next, and so on.

        # if add=TRUE, then par() is not refreshed.  allows multiple panels by
        # using, e.g., par(mfrow=c(2,2)).

	# op <- par(no.readonly = TRUE)

	ng0 <- length( survfit$strata )
	ng <- max(ng0,1) 
	# When only one group...
	if(ng0==0){
		survfit$strata <- length(survfit$time)
		names(survfit$strata) <- 'All'
		legend <- draw.lines <- FALSE
	} 

	lty.surv <- rep(lty.surv, ng)
	lwd.surv <- rep(lwd.surv, ng)
	col.surv <- rep(col.surv, ng)
	lty.ci <- rep(lty.ci, ng)
	lwd.ci <- rep(lwd.ci, ng)
	col.ci <- rep(col.ci, ng)

	## group names and error checking	
	gr <- c(survfit$strata)
	if( is.null(group.names) ) {
		group.names <- names(survfit$strata)
	}
	if( length(unique(group.names)) != ng ) {
		stop('\n','length(unique(group.names)) != number of groups.')
	}
	if( suppressWarnings(any( sort(group.order) != 1:ng)) ) {
		stop('\n', 'Something wrong with group.order.','\n','sort(group.order) must equal 1:', ng, '.')
	}
	group.names <- gsub(' *$', '', group.names)  #to remove unwanted white spaces in group.names.
	if(ng==1 & (group.names[1]=='group.names') ) {
		group.names <- 'N at risk'
		label.n.at.risk = FALSE
	}

	## graphic parameters
    if(!add) {
		par(list(oma=c(1,1,1,1), mar=c(4+ng,2+extra.left.margin,2,0)+.1))
		if(simple) par( mar=c(3,4,2,1)+.1 )
		par( list(...) )
    }

	## reformat survival estimates
	dat <- data.frame(time=survfit$time, n.risk=survfit$n.risk, n.event=survfit$n.event, survival=survfit$surv, std.err=survfit$std.err, 
						lower=survfit$lower, upper=survfit$upper, group=rep( group.names, gr) )
	dat.list <- split(dat, f=dat$group)

	## plot (but not survival curves) 
	plot(0, type='n', xlim=xlim, ylim=ylim, xaxt='n', yaxt='n', xaxs='i', yaxs="i",xlab='', ylab='' )
	if(grid){
		par('xpd'=FALSE)
		abline(v=xaxis.at, lty=lty.grid, lwd=lwd.grid, col=col.grid )
		abline(h=pretty(c(0,1)), lty=lty.grid, lwd=lwd.grid, col=col.grid )
	}
	axis( side=2, at=pretty(c(0,1)), cex.axis=cex.axis )	
	axis( side=1, at=xaxis.at, label=xaxis.lab, line=-0.5, tick=FALSE, cex.axis=cex.axis )
	axis( side=1, at=xaxis.at, label=rep('',length(xaxis.at)), line=0, tick=TRUE )
	title(xlab=xlab, line=1.5, adj=.5, ...) ; title(ylab=ylab, ... )

	if(!simple){
	## write group names
		group.name.pos <- (par()$usr[2]-par()$usr[1]) / -8 ; padding <- abs( group.name.pos / 8 )
		line.pos <- (1:ng)[order(group.order)] + 2
		mtext( group.names, side=1, line=line.pos, at=group.name.pos, adj=1, col=1, las=1, cex=cex.axis)
		## numbers at risk
		kms <- summary(survfit, times=xaxis.at) ; if(is.null(kms$strata)) kms$strata <- rep(1,length(kms$time) )
		d1 <- data.frame(time = kms$time, n.risk = kms$n.risk, strata = c(kms$strata))
		d2 <- split(d1, f=d1$strata)

		## Right-justifying the numbers 
		ndigits <- lapply(d2, function(x) nchar(x[,2]) )
		max.len <- max( sapply(ndigits, length) )
		L <- do.call('rbind', lapply(ndigits, function(z){ length(z) <- max.len ; z} ) )
		nd <- apply( L, 2, max, na.rm=T )
		for( i in seq(ng) ){
			this <- d2[[i]] 
			w.adj <- strwidth('0', cex=cex.axis, font=par('font')) / 2 * nd[1:nrow(this)]
			mtext( side=1, at=this$time+w.adj, text=this$n.risk, line=line.pos[i], cex=cex.axis, adj=1, col=1, las=1)
		}
		if(label.n.at.risk) mtext( side=1, text='No. at risk', at=group.name.pos, line=1.5, adj=1, col=1, las=1, cex=cex.axis )
	} ## End of if(!simple)	

	## pvalue
	sdiff <- survdiff(eval(survfit$call$formula), data = eval(survfit$call$data))
    pval <- pchisq(sdiff$chisq, length(sdiff$n)-1, lower.tail = FALSE)
	pvaltxt <- ifelse(pval < 0.001, "Log-rank P < 0.001", paste("Log-rank P =", round(pval, 3)))
	# Legend
	rlp <- group.order
	if(legend){
		bgc <- ifelse( par('bg')=='transparent', 'white', par('bg') )
		legend(x=loc.legend, legend=c(group.names[rlp],pvaltxt), col=c(col.surv[rlp],0), lty=lty.surv[rlp], lwd=lwd.surv[rlp],
			bty='o', cex=cex.axis, bg=bgc, box.col='transparent', inset=.01 ) #,text.font = 2
	}

	## draw confidence intervals
	##for(i in 1:ng){
		##this <- dat.list[[i]] 
		##x <- this$time
		##L <- this$lower
		##U <- this$upper
		##S <- this$survival
		##naL <- which( is.na(L) )
		##L[naL] <- L[naL-1]
		##U[naL] <- U[naL-1]
		##lines( x, L, type='s', col=col.ci[i], lty=lty.ci[i], lwd=lwd.ci[i] )
		##lines( x, U, type='s', col=col.ci[i], lty=lty.ci[i], lwd=lwd.ci[i] )
	##}
	# draw curves
	lines(survfit, conf.int=FALSE, col=col.surv, lty=lty.surv, lwd=lwd.surv, mark=mark, xmax=xlim[2], ymin=ylim[1])

	box(bty=par('bty'))
	
	# par(op) 
}
