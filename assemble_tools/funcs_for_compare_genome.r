fit_val <- function ( x=c(1:10), low=0, high=150 ) {
	x[ x>=low & x<=high ]
}

bwa.stats <- function ( x ) {
  back      <- list()
  back$qt   <- quantile( x, probs=c(0.25,0.75,0.99,0.97,0.95) ) 
  back$low  <- back$qt[1] - 2 * ( back$qt[2] - back$qt[1] )
  back$high <- back$qt[2] + 2 * ( back$qt[2] - back$qt[1] )
  back$used <- !(is.na(x)) & x >= back$low & x <= back$high
  back$mean <- mean( x[back$used] )
  back$sd   <- sd( x[back$used] )
  back$mean.norm <- mean( x )
  back$sd.norm <- sd( x )
  back
}# bwa.stats

plot_pair <- function ( bn6=bn6, id1="scaffold116_cov122", id2="scaffold12", xrange=NULL, yrange=NULL ) {
	kk <- bn6$V1 == id1 & bn6$V2 == id2
	# if ( is.null(xrange) ) { xrange <- range(bn6$V7[kk], bn6$V8[kk]) } 
	# if ( is.null(yrange) ) { yrange <- range(bn6$V9[kk], bn6$V10[kk]) } 
	if ( is.null(xrange) ) { xrange <- c(1, bn6$V13[kk][1]) } 
	if ( is.null(yrange) ) { yrange <- c(1, bn6$V14[kk][1]) } 
	plot(x=xrange, y=yrange, type='n', xlab=paste0(id1, " (", format(bn6$V13[kk][1], digits=0, big.mark=",")," bp)"), ylab=paste0(id2, " (", format(bn6$V14[kk][1], digits=0, big.mark=","), " bp)"))
	seg.col <- ifelse( bn6$V9[kk] > bn6$V10[kk], 'blue', 'darkred' )
	segments( x0=bn6$V7[kk], y0=bn6$V9[kk], x1=bn6$V8[kk], y1=bn6$V10[kk], col=seg.col )
	abline(v=c(1, bn6$V13[kk][1]), col="red")
	abline(h=c(1, bn6$V14[kk][1]), col="red")
}

plot_psl_pair <- function ( psl=psl, id1="scaffold116_cov122", id2="scaffold12", xrange=NULL, yrange=NULL ) {
	kk <- psl$V10 == id1 & psl$V14 == id2
	if ( is.null(xrange) ) { xrange <- c(1, psl$V11[kk][1]) } 
	if ( is.null(yrange) ) { yrange <- c(1, psl$V15[kk][1]) } 
	plot(x=xrange, y=yrange, type='n', xlab=paste0(id1, " (", format(psl$V11[kk][1], digits=0, big.mark=",")," bp)"), ylab=paste0(id2, " (", format(psl$V15[kk][1], digits=0, big.mark=","), " bp)"))
	psl$V12[kk] <- psl$V12[kk]+1; 
	psl$V16[kk] <- psl$V16[kk]+1; 
	mm <- psl$V9[kk] == '-' 
	mm.1 <- psl$V16[kk][mm]
	psl$V16[kk][mm] <- psl$V17[kk][mm]
	psl$V17[kk][mm] <- mm.1 
	seg.col <- ifelse( mm, 'purple', 'darkred' )
	segments( x0=psl$V12[kk], y0=psl$V16[kk], x1=psl$V13[kk], y1=psl$V17[kk], col=seg.col )
	abline(v=c(1, psl$V11[kk][1]), col="red")
	abline(h=c(1, psl$V15[kk][1]), col="red")
}

plot_psl_1_to_many <- function ( psl=psl, id1="scaffold116_cov122", id2=c("scaffold12"), xrange=NULL, yrange=NULL, bpPerInch=2e6, max_plot_inch = 6.9 ) {
	id2 <- as.character(id2) 
	id1 <- as.character(id1)
	id2.list <- list()
	id2.cumS <- list() 
	id2.cumS[[ id2[1] ]] <- 1 
	id2.len <- list() 
	id2.cumMax <- 0 
	for (id2.i in 1:length(id2)) {
		if ( id2[id2.i] == '' ) { next }
		id2.list[[ id2[id2.i] ]] <- id2.i
		id2.len[[ id2[id2.i] ]] <- psl$V15[ psl$V14 == id2[id2.i] ][1]
		id2.cumMax <- ifelse( id2.cumMax > id2.len[[ id2[id2.i] ]], id2.cumMax, id2.len[[ id2[id2.i] ]] )
		if ( id2.i > 1 ) {
			id2.cumS[[ id2[id2.i] ]] <- id2.cumS[[ id2[id2.i-1] ]] + id2.len[[ id2[id2.i-1] ]]
			id2.cumMax <- ifelse( id2.cumMax > id2.cumS[[ id2[id2.i] ]] + id2.len[[ id2[id2.i] ]] -1 , id2.cumMax, id2.cumS[[ id2[id2.i] ]] + id2.len[[ id2[id2.i] ]] -1 )
		}
	}
	rm(id2.i)
	
	kk <- psl$V10 == id1 & apply( as.matrix( psl$V14 ), MARGIN=1, FUN=function(x) { !is.null( id2.list[[ x[1] ]] ) } )
	if ( is.null(xrange) ) { xrange <- c(0, psl$V11[kk][1]) } 
	if ( is.null(yrange) ) { yrange <- c(0, id2.cumMax) } 
	par.xaxs <- par('xaxs')
	par.yaxs <- par('yaxs')
	par(xaxs= 'i')
	par(yaxs= 'i')
	par.pin <- par('pin') # 5.747554 5.151251
	par.mai <- par('mai')
	# par('mai' = c(par.mai[1:2], 0.2, 0.2))
	par.r2 <- par.mai[2] / par.pin[1]
	par.r4 <- par.mai[4] / par.pin[1]
	par.r1 <- par.mai[1] / par.pin[2]
	par.r3 <- par.mai[3] / par.pin[2]
	par.width  <- xrange[2] / ( bpPerInch * ( 1 - par.r2 - par.r4 ) )
	par.height <- yrange[2] / ( bpPerInch * ( 1 - par.r1 - par.r3 ) )
	par.hw_weight <- ifelse( max(par.width, par.height) <= max_plot_inch , 1, max_plot_inch/max(par.width, par.height) )
	par.width <- par.width * par.hw_weight + 1
	par.height <- par.height * par.hw_weight
	# par('pin' = c(par.width, par.height), 'mai' = c(par.mai[1:2], 0.2, 0.2))
	
	plot(x=xrange, y=yrange, type='n', xlab=paste0(id1, " (", format(psl$V11[kk][1], digits=0, big.mark=",")," bp)"), ylab=paste0(id2, " (", format(psl$V15[kk][1], digits=0, big.mark=","), " bp)", collapse="; "), 	)
	psl$V12[kk] <- psl$V12[kk]+1; 
	psl$V16[kk] <- psl$V16[kk]+1; 
	mm <- psl$V9[kk] == '-' 
	mm.1 <- psl$V16[kk][mm]
	psl$V16[kk][mm] <- psl$V17[kk][mm]
	psl$V17[kk][mm] <- mm.1 
	seg.col <- ifelse( mm, 'purple', 'darkred' )
	
	plot.y0 <- apply( X=cbind(psl$V14[kk], psl$V16[kk]), MARGIN=1, FUN=function (x) { as.numeric(x[2]) + id2.cumS[[ x[1] ]] - 1 }  )
	plot.y1 <- apply( X=cbind(psl$V14[kk], psl$V17[kk]), MARGIN=1, FUN=function (x) { as.numeric(x[2]) + id2.cumS[[ x[1] ]] - 1 }  )
	segments( x0=psl$V12[kk], y0=plot.y0, x1=psl$V13[kk], y1=plot.y1, col=seg.col )
	#segments( x0=psl$V12[kk], y0=psl$V16[kk], x1=psl$V13[kk], y1=psl$V17[kk], col=seg.col )
	abline(v=c(1, psl$V11[kk][1]), col="red")
	abline(h=c(unlist(id2.cumS), id2.cumMax), col="red")
	#abline(h=c(1, psl$V15[kk][1]), col="red")
	
	final.par.pin <- par('pin') 
	par('xaxs' = par.xaxs)
	par('yaxs' = par.yaxs)
	par('pin' = par.pin ) 
	par('mai' = par.mai)
	
	return( c(psl$V11[kk][1], id2.cumMax, par.width, par.height, final.par.pin) ) 
}# plot_psl_1_to_many


