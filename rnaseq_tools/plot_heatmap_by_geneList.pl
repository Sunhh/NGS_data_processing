#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"pdf_width:i", # 10
	"pdf_height:i", # 10
); 
$opts{'pdf_width'}  //= 10; 
$opts{'pdf_height'} //= 10; 
!@ARGV and die "perl $0 <exp.txt> <gene_list> <out_pref>\n"; 
my $expF = shift; 
my $genF = shift; 
my $pref = shift; 

my $wrkDir = &fileSunhh::new_tmp_dir( 'create' => 1 ); 

my $cmd = ''; 
my $i; 
my $hh = ''; 

my %eleID; 
$i=0; 
for my $lr (&fileSunhh::load_tabFile($genF, 0)) {
	$i++; 
	$eleID{$lr->[0]} //= $i; 
}

my %expData; 
$i=0; 
for my $lr (&fileSunhh::load_tabFile($expF, 1)) {
	$i++; 
	if ($i == 1) {
		$hh = join("\t", @$lr); 
		next; 
	}
	defined $eleID{$lr->[0]} or next; 
	$expData{$lr->[0]} = join("\t", @$lr); 
}

my $ofh = &openFH("$wrkDir/exp", '>'); 
print {$ofh} "$hh\n"; 
for my $id (sort { $eleID{$a} <=> $eleID{$b} } keys %expData) {
	print {$ofh} "$expData{$id}\n"; 
}
close($ofh); 

my $ofh2 = &openFH("$wrkDir/hmap1.r", '>'); 
print {$ofh2} <<"RRRR"; 
library(RColorBrewer)
library(pheatmap)
cc = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))

dat1 <- read.table( \"$wrkDir/exp\", header=T, sep=\"\\t\", row.names=1, stringsAsFactors=F )
mat1 <- log2( dat1 + 0.01 )
mat2 <- dat1
mat1.noRowSD <- 0 != apply(as.matrix(mat1), 1, sd, na.rm=T)
mat2.noRowSD <- 0 != apply(as.matrix(mat2), 1, sd, na.rm=T)
dist1_log2_noScale  <- dist( mat1, method = "euclidean" )
dist3_raw_noScale   <- dist( mat2, method = "euclidean" )
scale_rows <- function(x) {
  # From pheatmap.r 
  m = apply(x, 1, mean, na.rm = T); 
  s = apply(x, 1, sd, na.rm = T); 
  return((x-m)/s); 
}
dist2_log2_rowScale <- dist( scale_rows(mat1), method = "euclidean" )
dist4_raw_rowScale  <- dist( scale_rows(mat2), method = "euclidean" )

pdf( file="${pref}_hmaps.pdf", height=$opts{'pdf_height'}, width=$opts{'pdf_width'} )
# clustering_distance_rows 
pheatmap( mat1, main="heatmap with log2(EXP+0.01) [no scale]", 
	na_col = "grey", color=cc(100), 
	scale  = "none", 
	cluster_rows= T, cluster_cols=  F, 
	show_colnames=T, show_rownames= T, 
	# cellwidth = 20, cellheight = 20, 
	display_numbers = F, legend=T
)
pheatmap( mat1[mat1.noRowSD,], main="heatmap with log2(EXP+0.01) [row scaled]", 
	na_col = "grey", color=cc(100), 
	scale  = "row", 
	cluster_rows = T, cluster_cols = F, 
	show_colnames=T, show_rownames = T, 
	# cellwidth = 20, cellheight = 20, 
	display_numbers = F, legend=T
)
pheatmap( mat2, main="heatmap with EXP [no scale]", 
	na_col = "grey", color=cc(100), 
	scale  = "none", 
	cluster_rows = T, cluster_cols = F, 
	show_colnames=T, show_rownames = T, 
	# cellwidth = 20, cellheight = 20, 
	display_numbers = F, legend=T
)
pheatmap( mat2[mat2.noRowSD,], main="heatmap with EXP [row scaled]", 
	na_col = "grey", color=cc(100), 
	scale  = "row", 
	cluster_rows = T, cluster_cols = F, 
	show_colnames=T, show_rownames = T, 
	# cellwidth = 20, cellheight = 20, 
	display_numbers = F, legend=T
)
dev.off()

RRRR
close($ofh2); 
$cmd = "/home/Sunhh/bin/Rscript $wrkDir/hmap1.r"; 
&exeCmd_1cmd( $cmd ) and &stopErr("[Err] Failed at cmd: $cmd\n"); 

&fileSunhh::_rmtree( $wrkDir ); 

