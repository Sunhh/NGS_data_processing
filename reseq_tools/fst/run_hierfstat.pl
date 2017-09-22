#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use Getopt::Long; 
my %opts;
GetOptions(\%opts, 
	"help!", 
	"fst_in:s", 
	"mrk_info:s", 
	"inList!", 
	"maxNmissR:f", # 1 
	"minGrp1N:i",  # 1
	"minGrp2N:i",  # 1 
	"rmNegNeiFst!", 
	"rmNegWcFst!", 
	"exe_Rscript:s", # ~/bin/Rscript 
	"generateR:s", 
); 
$opts{'exe_Rscript'} //= '~/bin/Rscript'; 
$opts{'maxNmissR'} //= 1; 
$opts{'minGrp1N'} //= 1; 
$opts{'minGrp2N'} //= 1; 
my $help_txt = <<HH; 

perl $0 -fst_in in_fst_prefix -mrk_info mrk_info -exe_Rscript '~/bin/Rscript' 

-inList      [Bool] This option mask -mrk_info . 
-maxNmissR   [$opts{'maxNmissR'}] [0-1] Maximum N missing ratio allowed in each group. 
-minGrp1N    [$opts{'minGrp1N'}] [1-...] minimum genotyped individuals accepted in group 1. 
-minGrp2N    [$opts{'minGrp2N'}] [1-...] minimum genotyped individuals accepted in group 2.
-rmNegNeiFst [Bool] Set NA to sites with negative Nei_Fst if given. 
-rmNegWcFst  [Bool] Set NA to sites with negative Wc_Fst if given. 

-generateR   [''] Highest. 

-help

HH

$opts{'help'} and &LogInforSunhh::usage($help_txt); 
(defined $opts{'fst_in'}) or defined $opts{'generateR'} or &LogInforSunhh::usage($help_txt); 
( $opts{'maxNmissR'} >= 0 and $opts{'maxNmissR'} <= 1 ) or &LogInforSunhh::usage($help_txt); 
( $opts{'minGrp1N'}  >= 1 ) or &LogInforSunhh::usage($help_txt); 
( $opts{'minGrp2N'}  >= 1 ) or &LogInforSunhh::usage($help_txt); 

if (defined $opts{'generateR'} and $opts{'generateR'} ne '') {
	&_write_fst_R_byList( $opts{'generateR'} ); 
	exit(0); 
}

if ($opts{'inList'}) {
	&run_fst_R_byList(); 
	my $fst_in_aref = &load_list_info( $opts{'fst_in'} ); 
	for my $a1 (@$fst_in_aref) {
		my $mrk2loc_href = &load_mrk_info( $a1->[2] ); 
		&print_site_fst( "$a1->[1].fst.perSite", $mrk2loc_href,   "$a1->[1].fst.perSiteChrPos" ); 
		&print_wind_fst( "$a1->[1].fst.perWind", $a1->[0],        "$a1->[1].fst.perWindLine" ); 
	}
} else {
	&run_fst_R(); 
	my $mrk2loc_href = &load_mrk_info( $opts{'mrk_info'} ); 
	&print_site_fst( "$opts{'fst_in'}.fst.perSite", $mrk2loc_href,   "$opts{'fst_in'}.fst.perSiteChrPos" ); 
	&print_wind_fst( "$opts{'fst_in'}.fst.perWind", $opts{'fst_in'}, "$opts{'fst_in'}.fst.perWindLine" ); 
}

&tsmsg("[Rec] done.\n"); 

sub load_list_info {
	# $_[0] : list of fst_in 
	my $fh = &openFH($_[0], '<'); 
	my @dd; 
	while (<$fh>) {
		chomp; 
		my @ta = split(/\t/, $_); 
		push(@dd, [@ta[0,1,2]]); # [fst_in.file, fst_in.o_prefix, fst_in.mrk_info]
	}
	close($fh); 
	return \@dd; 
}

sub print_wind_fst {
	# $_[0] : fst_in.perWind
	# $_[1] : fst_in
	# $_[2] : fst_in.perWind.line
#  ''      x
#  Ho      0.0053
#  Hs      0.0935
	my $fh  = &openFH($_[0], '<'); 
	my $ofh = &openFH($_[2], '>'); 
	my (@h, @v); 
	while (<$fh>) {
		chomp; 
		my @ta = split(/\t/, $_); 
		$. == 1 and next; 
		push(@h, $ta[0]); 
		push(@v, $ta[1]); 
	}
	print {$ofh} join("\t", 'WindID',    @h)."\n"; 
	print {$ofh} join("\t", $_[1], @v)."\n"; 
	close($ofh); 
	close($fh); 
}
sub print_site_fst {
	# $_[0] : fst_in.perSite 
	# \%mrk2loc
	# $_[1] : fst_in.perSite.chrPos
#   ''      Ho      Hs      Ht      Dst     Htp     Dstp    Fst     Fstp    Fis     Dest     WCFst
#   L1      0       0.386   0.3813  -0.0047 0.3766  -0.0094 -0.0123 -0.025  1       -0.0153  -0.02
#   L2      0       0.0302  0.0303  0       0.0303  1e-04   0.001   0.0021  1       1e-04    0.0015
	my $fh  = &openFH($_[0], '<'); 
	my $ofh = &openFH($_[2], '>'); 
	while (<$fh>) {
		chomp; 
		my @ta = split(/\t/, $_); 
		if ($. == 1) {
			print {$ofh} join("\t", qw/chr pos/, @ta[1..$#ta])."\n"; 
			next; 
		}
		print {$ofh} join("\t", $_[1]->{$ta[0]}[0], $_[1]->{$ta[0]}[1], @ta[1..$#ta])."\n"; 
	}
	close($fh); 
	close($ofh); 
}

sub load_mrk_info {
# [Sunhh@Penguin test]$ head mrk_info
# L1      WM97_Chr01      11
# L2      WM97_Chr01      13
	my %mrk2loc; 
	my $fh = &openFH($_[0], '<'); 
	while (<$fh>) {
		chomp; 
		my @ta = split(/\t/, $_); 
		$mrk2loc{$ta[0]} = [@ta[1,2]]; 
	}
	close($fh); 
	return \%mrk2loc; 
}

sub run_fst_R_old {
	my $tmp_R = &fileSunhh::new_tmp_file(); 
	&_write_fst_R( $tmp_R ); 
	&exeCmd_1cmd("$opts{'exe_Rscript'} $tmp_R $opts{'fst_in'}") and &stopErr("[Err] Failed to run Rscript.\n"); 
	unlink($tmp_R); 
	return $tmp_R; 
}

sub run_fst_R {
	my $tmp_R = &fileSunhh::new_tmp_file(); 
	&_write_fst_R_byList( $tmp_R ); 
	my $tmp_in_list = &fileSunhh::new_tmp_file(); 
	&_write_fst_in_list( $tmp_in_list, $opts{'fst_in'} ); 
	&exeCmd_1cmd("$opts{'exe_Rscript'} $tmp_R $tmp_in_list") and &stopErr("[Err] Failed to run Rscript $tmp_R $tmp_in_list.\n"); 
	unlink($tmp_R); 
	unlink($tmp_in_list); 
	return $tmp_R; 
}

sub _write_fst_in_list {
	# $_[0] : out file name
	# $_[1] : fst_in.fmt file. 
	my $fh = &openFH($_[0], '>'); 
	print {$fh} join("\t", $_[1], $_[1], 'NA')."\n"; 
	close($fh); 
}# _write_fst_in_list() 

sub run_fst_R_byList {
	my $tmp_R = &fileSunhh::new_tmp_file(); 
	&_write_fst_R_byList( $tmp_R ); 
	&exeCmd_1cmd("$opts{'exe_Rscript'} $tmp_R $opts{'fst_in'}") and &stopErr("[Err] Failed in run_fst_R_byList() for $tmp_R $opts{'fst_in'}\n"); ; 
	unlink($tmp_R); 
	return($tmp_R); 
}

=head1 _write_fst_R( $outR_filename ) 
=cut
sub _write_fst_R {
	my $fh = &openFH($_[0], '>'); 
	print {$fh} <<'HH';

argvs <- commandArgs( trailingOnly=TRUE ) ; 
#flist <- read.table( file = argvs[1], stringsAsFactors=F, colClasses=c('character'), header=F ) ; 
# The first column is input file name, the second column is output file prefix. 
library(hierfstat); 

#for ( i in 1:nrow(flist) ) {
	# input file    : flist[i,1]
	# output prefix : flist[i,2]
	aa <- read.table( file= argvs[1], header=T, colClasses="numeric", stringsAsFactors=F ) 
	aa.stats <- basic.stats( aa )
	write.table( aa.stats$perloc,  file=paste0(argvs[1], ".fst.perSite", sep=""), append=F, row.names=T, col.names=NA, quote=F, sep="\t" )
	write.table( aa.stats$overall, file=paste0(argvs[1], ".fst.perWind", sep=""), append=F, row.names=T, col.names=NA, quote=F, sep="\t" )
#}


HH
	close($fh); 
	return ($_[0]); 
}# _write_fst_R() 

sub _write_fst_R_byList {

	my $fh = &openFH($_[0], '>'); 
	print {$fh} <<'L0'; 
maxNR <- 1  # maximum N missing ratio [0-1]
rmNegativeNeiFst <- FALSE # Tell if I should remove sites with nei_fst < 0 
rmNegativeWCFst  <- FALSE # Tell if I should remove sites with wc_fst < 0 

L0
	print {$fh} "maxNR  <- $opts{'maxNmissR'}\n"; 
	print {$fh} "minGrp1N <- $opts{'minGrp1N'}\n"; 
	print {$fh} "minGrp2N <- $opts{'minGrp2N'}\n"; 
	$opts{'rmNegNeiFst'} and print {$fh} "rmNegativeNeiFst <- TRUE\n"; 
	$opts{'rmNegWcFst'}  and print {$fh} "rmNegativeWCFst  <- TRUE\n"; 

	print {$fh} <<'L1'; 
.tsmsg <- function(...) {
	message("[", date(), "]: ", ...)
}# End .tsmsg

library(hierfstat);
argvs <- commandArgs( trailingOnly=TRUE ) ;
if ( length(argvs) == 0 ) {
	.tsmsg("[Err] Should be : ~/bin/Rscript   aa.R  file_list [AllowMaxNmissingR[0-1]   0_1_ifRmNegNeiFst   0_1_ifRmNegWcFst]\n")
	q() 
}
flist <- read.table( file = argvs[1], stringsAsFactors=F, colClasses=c('character'), header=F ) ;
# The first column is input file name, the second column is output file prefix.
maxNR <- ifelse( !is.na(argvs[2]) & !is.na(as.numeric(argvs[2])), as.numeric(argvs[2]), maxNR )
minGrp1N <- minGrp1N # I don't bother setting this parameter. 
minGrp2N <- minGrp2N # I don't bother setting this parameter. 
rmNegativeNeiFst <- switch( as.character(argvs[3]), '0'=FALSE, '1'=TRUE, 'TRUE'=TRUE, 'FALSE'=FALSE, 'T'=TRUE, 'F'=FALSE, rmNegativeNeiFst )
rmNegativeWCFst  <- switch( as.character(argvs[4]), '0'=FALSE, '1'=TRUE, 'TRUE'=TRUE, 'FALSE'=FALSE, 'T'=TRUE, 'F'=FALSE, rmNegativeWCFst  )

# Output : list( qw/ nei_perloc wc_perloc_fst nei_overall wc_overall_fst / )
cnt_fst <- function( x=NULL, stats=NULL, maxNR=1, minGrp1N=1, minGrp2N=1, x.kk=NULL, rmNegativeNeiFst=FALSE, rmNegativeWCFst=FALSE ) {
	# Note that length(x.kk) is always 1 less than nrow(x)
	#   and x.kk should be a vector of TRUE/FALSE ; 
	# Output : list( qw/ nei_perloc wc_perloc_fst nei_overall wc_overall_fst / )
	back <- list() 
	if (is.null(x)) {
		.tsmsg("No input for x!\n")
		q()
	}
	samp.ttl   <- nrow(x)
	site.ttl   <- ncol(x)-1
	samp.num   <- as.numeric( table( x[,1] ) )
	site.name  <- colnames(x) ; site.name  <- site.name[-1] ; 
	if (is.null(x.kk)) { x.kk <- rep(TRUE, times=site.ttl) }
	if ( length(x.kk) != site.ttl ) {
		.tsmsg("[Err] Bad length of x.kk:",length(x.kk), " VS. site.total:", site.ttl, "\n"); 
		q(); 
	}
	
	nei_perloc_colName  <- c('Ho', 'Hs', 'Ht', 'Dst', 'Htp', 'Dstp', 'Fst', 'Fstp', 'Fis', 'Dest', 'Hs_p1', 'Hs_p2')
	back$nei_perloc     <- matrix( rep(NaN, times=length(nei_perloc_colName)*site.ttl), byrow=TRUE, nrow=site.ttl )
	rownames(back$nei_perloc) <- site.name
	colnames(back$nei_perloc) <- nei_perloc_colName
	back$nei_overall    <- rep( NaN, times=length(nei_perloc_colName))
	names( back$nei_overall ) <- nei_perloc_colName
	back$wc_perloc_fst  <- rep( NaN, times=site.ttl)
	back$wc_overall_fst <- NaN
	back$error          <- NULL 
	
	# Remove sites with only NA and set it as FALSE in x.kk 
	if ( sum(x.kk) == 0 ) {
		return(back); 
	} else if ( sum(x.kk) == 1 ) {
		na.inSite <- sum( is.na( x[, c(which(x.kk)+1)] ) )
		x.kk[x.kk][na.inSite == samp.ttl] <- FALSE
	} else {
		na.inSite <- apply( x[, c(which(x.kk)+1)], MARGIN=2, FUN=function(y) { sum(is.na(y)) } )
		x.kk[x.kk][na.inSite == samp.ttl] <- FALSE
	}
	
	# Check how many available sites in x.kk and determine what to do for basic.stats() 
	if ( sum(x.kk) == 0 ) {
		; 
	} else if ( sum(x.kk) == 1 ) {
		x1 <- x[, c(1, which(x.kk)+1, which(x.kk)+1)]
		tmp.back <- cnt_fst( x=x1, stats=NULL, maxNR=maxNR, minGrp1N=minGrp1N, minGrp2N=minGrp2N, x.kk=c(TRUE,TRUE), rmNegativeNeiFst=rmNegativeNeiFst, rmNegativeWCFst=rmNegativeWCFst )
		back$nei_overall          <- tmp.back$nei_overall 
		back$wc_overall_fst       <- tmp.back$wc_overall_fst
		back$wc_perloc_fst[x.kk]  <- tmp.back$wc_perloc_fst[1]
		back$error                <- tmp.back$error 
		
		tmp <- matrix(tmp.back$nei_perloc[1,], nrow=1, byrow=TRUE)
		colnames(tmp) <- colnames(tmp.back$nei_perloc)
		rownames(tmp) <- rownames(tmp.back$nei_perloc)[1]
		back$nei_perloc[x.kk, ] <- tmp 
	} else if ( sum(x.kk) >= 2 ) {
		if (is.null(stats)) {
			x1 <- x[, c(1, which(x.kk)+1)]
			x1.stats <- basic.stats(x1)
			tmp.back <- cnt_fst( x=x1, stats=x1.stats, maxNR=maxNR, minGrp1N=minGrp1N, minGrp2N=minGrp2N, rmNegativeNeiFst=rmNegativeNeiFst, rmNegativeWCFst=rmNegativeWCFst )
			back$wc_overall_fst         <- tmp.back$wc_overall_fst
			back$nei_overall            <- tmp.back$nei_overall
			back$wc_perloc_fst[x.kk]    <- tmp.back$wc_perloc_fst
			back$nei_perloc[x.kk, ]     <- tmp.back$nei_perloc
			back$error                  <- tmp.back$error 
		} else {
			good1 <- sum(x.kk) # Must have ( good2 >= 2 ) 
			need_num1 <- samp.num[1]-floor(maxNR*samp.num[1])
			need_num2 <- samp.num[2]-floor(maxNR*samp.num[2])
			enough_num1 <- sapply( stats$n.ind.samp[,1], FUN=function(y) {isTRUE(y >= need_num1 & y >= minGrp1N)} )
			enough_num2 <- sapply( stats$n.ind.samp[,2], FUN=function(y) {isTRUE(y >= need_num2 & y >= minGrp2N)} )
			# Remove sites without variations, because it will cause wc() fail if there isn't any variation in dataset. 
			var_site <- apply( x[,-1], MARGIN=2, FUN= function( x ) { dim(table(x)) > 1 } )
			x.kk <- x.kk & enough_num1 & enough_num2 & var_site
			if (isTRUE( rmNegativeNeiFst )) {
				x.kk <- x.kk & sapply( stats$perloc[,7], FUN=function(y) {isTRUE(y >= 0)} )
			}
			good2 <- sum(x.kk) # Maybe 0,1,2 .. good2 
			x1 <- x[, c(1, which(x.kk)+1)]
			if (good2 == 0) {
				; 
			} else if (good1 == good2) {
				x1.wc <- wc(x1) 
				if ( length(x1.wc$per.loc$FST) != sum(x.kk) ) {
					.tsmsg("[Err] The wc() return different length of fst. Need to check input data\n")
					back$error <- paste0("wc() return different length rawN=[", sum(x.kk),'] wcN=[',length(x1.wc$per.loc$FST),']', sep="")
					return(back)
				}
				if ( isTRUE(rmNegativeWCFst) ) {
					x1.kk <- sapply( x1.wc$per.loc$FST, FUN=function(y) {isTRUE(y >= 0)} )
					x1.cnt_fst <- cnt_fst( x=x1, stats=NULL, x.kk=x1.kk, maxNR=maxNR, minGrp1N=minGrp1N, minGrp2N=minGrp2N, rmNegativeNeiFst=rmNegativeNeiFst, rmNegativeWCFst=FALSE )
					back$wc_overall_fst      <- x1.cnt_fst$wc_overall_fst
					back$nei_overall         <- x1.cnt_fst$nei_overall
					back$wc_perloc_fst[x.kk] <- x1.cnt_fst$wc_perloc_fst
					back$nei_perloc[x.kk,]   <- x1.cnt_fst$nei_perloc
					back$error               <- x1.cnt_fst$error
				} else {
					back$wc_overall_fst      <- x1.wc$FST
					back$nei_overall         <- c( stats$overall, Hs_p1=mean( stats$Hs[,1] , na.rm=TRUE), Hs_p2=mean( stats$Hs[,2], na.rm=TRUE ) )
					back$wc_perloc_fst[x.kk] <- x1.wc$per.loc$FST
					back$nei_perloc          <- cbind( as.matrix( stats$perloc ), Hs_p1=stats$Hs[,1], Hs_p2=stats$Hs[,2] )
				}
			} else if (good2 < good1) {
				tmp.back <- cnt_fst( x=x1, stats=NULL, maxNR=maxNR, minGrp1N=minGrp1N, minGrp2N=minGrp2N, rmNegativeNeiFst=rmNegativeNeiFst, rmNegativeWCFst=rmNegativeWCFst )
				back$wc_overall_fst      <- tmp.back$wc_overall_fst
				back$nei_overall         <- tmp.back$nei_overall
				back$wc_perloc_fst[x.kk] <- tmp.back$wc_perloc_fst
				back$nei_perloc[x.kk, ]  <- tmp.back$nei_perloc
				back$error               <- tmp.back$error
			} else {
				.tsmsg(paste0("Why come here? good2=", good2, " good1=", good1, sep=""))
				q()
			}
		}
	} else {
		.tsmsg(paste0("Why come here? good2=", good2, " good1=", good1, sep=""))
		q()
	}
	
	back ; 
}# cnt_fst()

#### The main body #### 
for ( i in 1:nrow(flist) ) {
		# input file    : flist[i,1]
		# output prefix : flist[i,2]
		aa <- read.table( file= flist[i,1], header=T, colClasses="numeric", stringsAsFactors=F )
		aa.cnt_fst <- cnt_fst( x=aa, stats=NULL, maxNR=maxNR, minGrp1N=minGrp1N, minGrp2N=minGrp2N, rmNegativeNeiFst=rmNegativeNeiFst, rmNegativeWCFst=rmNegativeWCFst )
		if (!is.null(aa.cnt_fst$error)) {
			.tsmsg("[Err] Bad return of cnt_fst() for file [", flist[i,1],"]; error is \"", aa.cnt_fst$error, "\"\n")
			q()
		}
		
		out.perloc  <- cbind( aa.cnt_fst$nei_perloc, WCFst=aa.cnt_fst$wc_perloc_fst )
		out.overall <- c( aa.cnt_fst$nei_overall, WCFst=aa.cnt_fst$wc_overall_fst )

		if ( ncol(out.perloc) > 11 ) {
			out.perloc <- out.perloc[ ,c( 1:10, ncol(out.perloc) , 11:(ncol(out.perloc)-1) ) ]
		}
		if ( length(out.overall) > 11 ) {
			out.overall <- out.overall[ c( 1:10, length(out.overall), 11:(length(out.overall)-1) ) ] 
		}

		write.table( out.perloc,  file=paste0(flist[i, 2], ".fst.perSite", sep=""), append=F, row.names=T, col.names=NA, quote=F, sep="\t" )
		write.table( out.overall, file=paste0(flist[i, 2], ".fst.perWind", sep=""), append=F, row.names=T, col.names=NA, quote=F, sep="\t" )
}

#########################################################
##### Some information
## About the wc() : If all of the samples' alleles in one site are NA, wc() will skip this stie. 
##                  If all of a sub-pop's alleles in one site are NA, wc() may produce bad estimate on that site. 
##                  So I want to remove any site in the above two cases. 
## About basic.stats() : When there are too many sites in the above two case, basic.stats() will produce warnings(). But the matrix is full size. 
## 
#########################################################

#########################################################
##### This is a method using apply to do iteration, but it appears slower than using for
#cnt_fst <- function ( xx ) {
#       aa <- read.table( file=xx[1], header=T, colClasses="numeric", stringsAsFactors=F )
#       aa.stats <- basic.stats( aa )
#       write.table( aa.stats$perloc,  file=paste0(xx[2], ".fst.perSite", sep=""), append=F, row.names=T, col.names=NA, quote=F, sep="\t" )
#       write.table( aa.stats$overall, file=paste0(xx[2], ".fst.perWind", sep=""), append=F, row.names=T, col.names=NA, quote=F, sep="\t" )
#}
#apply( X=flist, MARGIN=1, FUN=cnt_fst )
#####
#########################################################


L1
	close($fh); 
	return ($_[0]); 
}

