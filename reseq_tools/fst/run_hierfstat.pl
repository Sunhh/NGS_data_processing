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
	"exe_Rscript:s", # ~/bin/Rscript 
); 
$opts{'exe_Rscript'} //= '~/bin/Rscript'; 
$opts{'maxNmissR'} //= 1; 

my $help_txt = <<HH; 

perl $0 -fst_in in_fst_prefix -mrk_info mrk_info -exe_Rscript '~/bin/Rscript' 

-inList      [Bool] This option mask -mrk_info . 
-maxNmissR   [$opts{'maxNmissR'}] [0-1] Maximum N missing ratio allowed in each group. 

-help

HH

$opts{'help'} and &LogInforSunhh::usage($help_txt); 
(defined $opts{'fst_in'}) or &LogInforSunhh::usage($help_txt); 

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
#   ''      Ho      Hs      Ht      Dst     Htp     Dstp    Fst     Fstp    Fis     Dest
#   L1      0       0.386   0.3813  -0.0047 0.3766  -0.0094 -0.0123 -0.025  1       -0.0153
#   L2      0       0.0302  0.0303  0       0.0303  1e-04   0.001   0.0021  1       1e-04
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

sub run_fst_R {
	my $tmp_R = &fileSunhh::new_tmp_file(); 
	&_write_fst_R( $tmp_R ); 
	&exeCmd_1cmd("$opts{'exe_Rscript'} $tmp_R $opts{'fst_in'}"); 
	unlink($tmp_R); 
	return $tmp_R; 
}

sub run_fst_R_byList {
	my $tmp_R = &fileSunhh::new_tmp_file(); 
	&_write_fst_R_byList( $tmp_R ); 
	&exeCmd_1cmd("$opts{'exe_Rscript'} $tmp_R $opts{'fst_in'}"); 
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
	print {$fh} <<'LL'; 
argvs <- commandArgs( trailingOnly=TRUE ) ; 
flist <- read.table( file = argvs[1], stringsAsFactors=F, colClasses=c('character'), header=F ) ; 
# The first column is input file name, the second column is output file prefix. 
library(hierfstat); 

#########################################################
##### This is a method using apply to do iteration, but it appears slower than using for
#cnt_fst <- function ( xx ) {
#	aa <- read.table( file=xx[1], header=T, colClasses="numeric", stringsAsFactors=F )
#	aa.stats <- basic.stats( aa )
#	write.table( aa.stats$perloc,  file=paste0(xx[2], ".fst.perSite", sep=""), append=F, row.names=T, col.names=NA, quote=F, sep="\t" )
#	write.table( aa.stats$overall, file=paste0(xx[2], ".fst.perWind", sep=""), append=F, row.names=T, col.names=NA, quote=F, sep="\t" )
#}
#apply( X=flist, MARGIN=1, FUN=cnt_fst )
#####
#########################################################


for ( i in 1:nrow(flist) ) {
	# input file    : flist[i,1]
	# output prefix : flist[i,2]
	aa <- read.table( file= flist[i,1], header=T, colClasses="numeric", stringsAsFactors=F ) 
	np.ttl <- nrow(aa); 
	samp.num <- as.numeric( table( aa[,1] ) ) 
	aa.stats <- basic.stats( aa )
	aa.stats.perloc  <- as.matrix( aa.stats$perloc ) 
	aa.stats.overall <- aa.stats$overall 
	aa.wc <- wc( aa )
	aa.wc.perloc.fst <- aa.wc$per.loc$FST 
	aa.wc.overall.fst <- aa.wc$FST 

LL

print {$fh} <<"L1";
	aa.kk <- aa.stats\$n.ind.samp[,1] >= ceiling( (1-$opts{'maxNmissR'})*samp.num[1] ) \& aa.stats\$n.ind.samp[,2] >= ceiling( (1-$opts{'maxNmissR'})*samp.num[2] ) \& !is.na(aa.stats\$n.ind.samp[,1]) \& !is.na(aa.stats\$n.ind.samp[,2])
L1
print {$fh} <<'L2'; 
	if ( sum(aa.kk) == 0 ) {
		aa.stats.perloc[ !is.nan( aa.stats.perloc ) ] <- NaN
		aa.stats.overall[ !is.na( aa.stats.overall ) ] <- NaN 
		aa.wc.perloc.fst[ !is.na( aa.wc.perloc.fst ) ] <- NaN 
		aa.wc.overall.fst <- NaN 
	} else {
		if ( sum(aa.kk) == 1 ) {
			bb <- aa[ ,c(1, which(aa.kk)+1, which(aa.kk)+1) ] 
		} else {
			bb <- aa[ ,c(1, which(aa.kk)+1) ] 
		}
		bb.stats <- basic.stats( bb ) 
		bb.wc <- wc( bb )
		aa.stats.perloc.toNA <- aa.stats.perloc[ !aa.kk, ]
		aa.stats.perloc.toNA[ !is.nan( aa.stats.perloc.toNA ) ] <- NaN
		aa.stats.perloc[ !aa.kk, ] <- aa.stats.perloc.toNA
		aa.stats.overall <- bb.stats$overall
		aa.wc.perloc.fst[ !aa.kk ] <- NaN
		aa.wc.overall.fst <- bb.wc$FST
	}
	out.perloc <- cbind( aa.stats.perloc, WCFst=aa.wc.perloc.fst )
	out.overall <- c( aa.stats.overall, WCFst=aa.wc.overall.fst )

	write.table( out.perloc,  file=paste0(flist[i, 2], ".fst.perSite", sep=""), append=F, row.names=T, col.names=NA, quote=F, sep="\t" )
	write.table( out.overall, file=paste0(flist[i, 2], ".fst.perWind", sep=""), append=F, row.names=T, col.names=NA, quote=F, sep="\t" )
}
L2
	close($fh); 
	return ($_[0]); 
}

