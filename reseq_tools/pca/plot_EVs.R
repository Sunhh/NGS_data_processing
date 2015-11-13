argvs <- commandArgs( trailingOnly=TRUE ) ;

tblFn <- argvs[1]
legendTag <- argvs[2]
grpLis <- argvs[-c(1,2)]
# tblFn <- 'in_pca.set01_all.pca.evec.tbl'

aa <- read.table( file=tblFn, header=T, stringsAsFactors=F )
if (length(grpLis) == 0) {
        grpLis <- as.character( levels( as.factor( aa$Grp ) ) ) 
}
if (is.null(legendTag) || is.na(legendTag)) {
	legendTag <- 'topleft'
}
grpNum <- length( grpLis ) 
grpColor <- rainbow( grpNum )
grpShape <- sample(1:20, size=grpNum)

kk <- aa$Grp %in% grpLis
indvColor <- rep('while', times=length(aa$Indv))
indvShape <- rep(20,      times=length(aa$Indv))
for ( i in 1:length(grpLis) ) {
	indvColor[ grpLis[i] == aa$Grp ] <- grpColor[i]
	indvShape[ grpLis[i] == aa$Grp ] <- grpShape[i]
}
pdf( file = "plot_EVs.pdf" )
plot( aa$EV1[kk] , aa$EV2[kk], col=indvColor[kk], pch=indvShape[kk], xlab="EV1", ylab="EV2" )
legend( legendTag , legend=grpLis, col=grpColor, pch=grpShape)
write.table( cbind("Grp"=grpLis, "GrpColor"=grpColor, "GrpShape"=grpShape), file="plot_EVs.para_plot", quote=FALSE, sep="\t", row.names=FALSE )
dev.off() 

