#!/usr/bin/Rscript

fixNegLog10P <- function (x) {
	x[ is.infinite(x) & !is.na(x) ] <- 1 + max( x[ !is.infinite(x) ],na.rm=TRUE )
	x[ is.na(x) ] <- 0
	x
}

argvs <- commandArgs( trailingOnly=TRUE ) ;

inF <- as.character( argvs[1] ) ; # This is the input file for Gprime, with only one colum and without header. 

x_Gp <- scan( inF )

# 1. Let W_Gp = log(X_Gp), with value name 'lnGp'; 
lnGp   <- log(x_Gp)

# 2. Let s_W = MAD_l(W_Gp), MAD_l(Y) = median(abs(y_i - median(Y))) , where y_i <= median(Y)
med_lnGp <- median(lnGp, na.rm=TRUE) 
s_W      <- median( abs( lnGp[ lnGp <= med_lnGp ] - med_lnGp ), na.rm=TRUE )

# 3. Use Hampel’s rule to identify outliers, O_W , as all w_i in W_Gp that satisfy : w_i - median(W_Gp) > g(N, alpha_N) * MAD_l(W_Gp)
#   MAD_l(W_Gp) is s_W ; 
#   O_W is a data set of outliers, which satisfies the above formula; 
#   w_i is one of W_Gp, which is value lnGp; 
k_outliers <- ( lnGp - med_lnGp > 5.2 * s_W ) | x_Gp == 0 | is.na( lnGp ) | is.infinite( lnGp )

# 4. Construct a trimmed data set X_T that doesn't have data from O_W; 
in_Gp <- x_Gp[ !k_outliers ]

# 5. Let theta_mu_bar = log(median(X_T)) , and theta_sd_bar = theta_mu_bar - log(mode_r(X_T)) ; 
#    mode_r(X_T) is a robust estimator of the mode for continuous variables; 
#    I saw the QTLseqr package used half-sample mode as this estimator. 
mode_r <- modeest::mlv( x= in_Gp, bw= 0.5, method= 'hsm' )$M
# mode_r <- modeest::mlv( x= in_Gp, bw= 0.3, method= 'hrm' )$M
mu_bar <- log( median(in_Gp) )  ;               # This is the mean of log scaled distribution of Gprime; 
sd_bar <- sqrt( abs(mu_bar - log( mode_r )) ) ; # This is the standard deviation of log scaled distribution of Gprime; 

# 6. Get p-values; 
pval   <- 1 - plnorm( q= x_Gp, meanlog= mu_bar, sdlog= sd_bar, lower.tail=TRUE ); # lower.tail=TRUE, because the outlier_Gprime should be bigger than null_in_Gprime; 
negLog10P <- -log10(pval)
# 7. Simply adjusted P value; 
qval   <- p.adjust( p= pval, method= 'BH' ) # BH - Benjamini & Hochberg, alias 'fdr'
negLog10AdjP <- -log10(qval)

# 8. Write to file. 
negLog10P    <- fixNegLog10P(negLog10P)
negLog10AdjP <- fixNegLog10P(negLog10AdjP) 
ooo    <- cbind( pval, negLog10P, qval, negLog10AdjP )
write.table( ooo, file= paste0(inF, '.pval', sep=''), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE )

# 9. Plot distribution; 
rrr    <- range( x_Gp[ !is.na(x_Gp) ] )
rrr[1] = round(rrr[1]-1, digits= 0)
rrr[2] = round(rrr[2]+1, digits= 0)
simu_xx <- seq(0, rrr[2], 0.001)
simu_yy <- dlnorm( simu_xx, meanlog= mu_bar, sdlog= sd_bar )
plot_xmax <- qlnorm(0.999,  meanlog= mu_bar, sdlog= sd_bar )[1]
simu_q01  <- qlnorm(0.99,   meanlog= mu_bar, sdlog= sd_bar )[1]
pdf( file=paste0(inF, '.dist.pdf', sep='') )
hist( x_Gp, breaks=seq(rrr[1], rrr[2], 0.01), xlim=c(0, plot_xmax), freq=FALSE )
points( simu_xx, simu_yy, type='l', col='red' )
legend( 'topright', legend=c(paste0('mean=',round(exp(mu_bar), digits=4)), paste0('mode=', round(mode_r, digits=4)), paste0('top 1%=', round(simu_q01, digits=4))), col=c('red', 'green', 'blue'), pch=20 )
abline( v=c(exp(mu_bar), mode_r, simu_q01), col=c('red', 'green', 'blue') )
dev.off()



