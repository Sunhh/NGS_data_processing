#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 

!@ARGV and die "perl $0 background_pwy subset_geneList [padj_cutoff_0.05]\n"; 

my $fn_bg = shift; 
my $fn_sub = shift; 
my $padj_cutoff = shift; 
$padj_cutoff //= 0.05; 

# ($bg_total, $bg_pwy, $sub_total, $sub__pwy)
my %bg = %{&_readInBg($fn_bg)}; 

my %sublis; 
my $fh_sub = &openFH($fn_sub); 
while (&wantLineC($fh_sub)) {
	my @ta = &splitL("\t", $_); 
	$sublis{'geneNum'} ++; 
	defined $bg{'gene2pwy'}{$ta[0]} or next; 
	for my $pwy (keys %{$bg{'gene2pwy'}{$ta[0]}}) {
		$sublis{'chk_pwy'}{$pwy} = $.; 
		$sublis{'pwy2gene'}{$pwy}{$ta[0]} = $.; 
	}
}
close ($fh_sub); 

my $wrk_dir = &fileSunhh::new_tmp_dir(); 
mkdir($wrk_dir); 
for my $pwy (sort keys %{$sublis{'chk_pwy'}}) {
	my $bg_total = $bg{'geneNum'}; 
	my $bg_pwy = scalar(keys %{$bg{'pwy2gene'}{$pwy}}); 
	my $sub_total = $sublis{'geneNum'}; 
	my $sub_pwy = scalar(keys %{$sublis{'pwy2gene'}{$pwy}}); 
	&fileSunhh::write2file( "$wrk_dir/input", join("\t", $pwy, $bg_total, $bg_pwy, $sub_total, $sub_pwy, join(";", sort keys %{$sublis{'pwy2gene'}{$pwy}}))."\n", '>>' ); 
}
my $fh_r = &openFH("$wrk_dir/c.r", '>'); 
print {$fh_r} "library(stats)\n"; 
print {$fh_r} "tbl <- read.table( \"$wrk_dir/input\", header=FALSE, sep=\"\\t\", stringsAsFactors=FALSE, quote=\"\" )\n"; 
print {$fh_r} "pval <- apply(tbl, MARGIN=1, FUN=function (x) { q=as.numeric(x[5])-1 ; m=as.numeric(x[3]); n=as.numeric(x[2])-as.numeric(x[3]); k=as.numeric(x[4]); phyper( q=q, m, n, k, lower.tail = FALSE ) })\n";
print {$fh_r} "padj <- p.adjust( pval, method= \'BH\' )\n";
print {$fh_r} "oo <- cbind(padj, pval, tbl)[padj < $padj_cutoff, ]\n";
print {$fh_r} "oo <- oo[ order(oo[,1], decreasing=FALSE), ]\n";
print {$fh_r} "write.table( oo, file=\"$wrk_dir/output\", append=FALSE, quote=FALSE, sep=\"\\t\", row.names=FALSE, col.names=FALSE )\n";
close($fh_r); 
system("Rscript $wrk_dir/c.r"); 

my $fh_o = &openFH("$wrk_dir/output", '<'); 
while (<$fh_o>) {
	print STDOUT $_; 
}
close($fh_o); 

&fileSunhh::_rmtree($wrk_dir); 



sub _readInBg {
	my $fn = shift; 
	my $fh = &openFH($fn); 
	my %back; 
	while (&wantLineC($fh)) {
		my @ta = &splitL("\t", $_); # gene_ID, pwy_ID!!KO_ID!!pwy_desc;;pwy_ID!!KO_ID!!pwy_desc...
		$back{'geneNum'} ++; 
		( defined $ta[1] and $ta[1] !~ m/^\s*$/ ) or next; 
		for my $tb ( &splitL(";;", $ta[1]) ) {
			my @tc = split(/!!/, $tb); 
			$tb = "$tc[0] ($tc[2])"; 
			$back{'gene2pwy'}{$ta[0]}{$tb} = $.; 
			$back{'pwy2gene'}{$tb}{$ta[0]} = $.; 
		}
	}
	close($fh); 

	return(\%back); 
}

