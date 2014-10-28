#!/usr/bin/perl -w
# Design: 
#  Input files : 
#    gff      - PG1All_v2.scf.fa.gff85  
#    full_LTR - PG1All_v2.scf.fa.out85
#    ScfFa    - PG1All_v2.scf.fa
#    refGff   - PG1All_v2.scf.fa.gff99
#    full_LTR - PG1All_v2.scf.fa.out99
#  Tools : 
#    ch_gff_to_tab.pl : 
#    ch_seqID.pl      : 
#    blastn           : 
#    makeblastdb      : 

use strict;
use warnings; 
use LogInforSunhh;
use Cwd 'abs_path';
use File::Basename;


!@ARGV and die "perl $0 in_gff in_gff.LTR.fa ScfFa ref_gff ref_gff.LTR.fa\n"; 
my ($inGff, $inGffSeq, $inScfFa, $refGff, $refGffSeq) = @ARGV; 

## Setting tools being used. 
my %tool; 
{
$tool{pathCfg_dir} = dirname( abs_path($0) );
$tool{pathCfg_file} = "$tool{pathCfg_dir}/path.conf";

&getPath(\%tool, $tool{pathCfg_file});


#$tool{exe_makeblastdb} = '/usr/local/bin/makeblastdb'; 
#$tool{exe_blastn} = '/usr/local/bin/blastn'; 
#$tool{exe_mv} = '/usr/bin/mv'; 
#my $dd = '/workdir/laopopo/spinach/genome/Repeat/05.LTR/tools'; 
#$tool{pl_ch_gff_to_tab} = "$dd/ch_gff_to_tab.pl"; 
#$tool{pl_ch_seqID} = "$dd/ch_seqID.pl"; 
#my $d2 = '/home/laopopo/tools/github/NGS_data_processing'; 
#$tool{pl_deal_fasta} = "$d2/deal_fasta.pl"; 
#$tool{pl_deal_table} = "$d2/deal_table.pl"; 
}

my %para; 
$para{minCov} = 90; 

## Making files for analysis 
&exeCmd("perl $tool{pl_ch_gff_to_tab} $refGff 1> $refGff.tmp.tab"); 
&exeCmd("perl $tool{pl_ch_seqID} $refGffSeq $refGff.tmp.tab 1>$refGffSeq.tmp.full.fa 2>$refGff.tmp.tab1"); 
&exeCmd("$tool{exe_mv} $refGff.tmp.tab1 $refGff.tmp.tab"); 
&exeCmd("perl $tool{pl_deal_fasta} $inScfFa -drawByList -drawList $refGff.tmp.tab -drawLcol 15,5,6,3 > $refGff.tmp.ltr.fa"); 
&exeCmd("perl $tool{pl_deal_fasta} $inScfFa -drawByList -drawList $refGff.tmp.tab -drawLcol 15,7,8,3 >> $refGff.tmp.ltr.fa"); 

&exeCmd("perl $tool{pl_ch_gff_to_tab} $inGff 1> $inGff.tmp.tab"); 
&exeCmd("perl $tool{pl_ch_seqID} $inGffSeq $inGff.tmp.tab 1>$inGffSeq.tmp.full.fa 2>$inGff.tmp.tab1"); 
&exeCmd("$tool{exe_mv} $inGff.tmp.tab1 $inGff.tmp.tab"); 
&exeCmd("perl $tool{pl_deal_fasta} $inScfFa -drawByList -drawList $inGff.tmp.tab -drawLcol 15,5,6,3,0 > $inGff.tmp.ltr.raw.fa"); 
&exeCmd("perl $tool{pl_deal_fasta} $inScfFa -drawByList -drawList $inGff.tmp.tab -drawLcol 15,7,8,3,0 >> $inGff.tmp.ltr.raw.fa"); 

{
&tsmsg("[Rec] Change names in $inGff.tmp.ltr.fa\n"); 
# Change names in "$inGff.tmp.ltr.fa"
open F,'<',"$inGff.tmp.ltr.raw.fa" or &stopErr("[Err] Failed to open $inGff.tmp.ltr.raw.fa $!\n") ;
open O,'>',"$inGff.tmp.ltr.fa" or &stopErr("[Err] Failed open $!\n"); 
my %h; 
while (<F>) {
	chomp; s/\s+$//; 
	if (m/^\s*>/) {
		s/>(\S+)// or &stopErr("[Err] $_\n"); 
		my $tk = $1; 
		if (defined $h{$tk}) {
			$_ = ">${tk}_1"; 
		}else{
			$_ = ">${tk}_2"; 
			$h{$tk} = 1; 
		}
	}
	print O "$_\n"; 
}
close O; 
close F; 
}

&tsmsg("[Rec] Run blastn from $inGff.tmp.ltr.fa to $refGff.tmp.ltr.fa\n"); 
&exeCmd("$tool{exe_makeblastdb} -dbtype nucl -in $refGff.tmp.ltr.fa"); 
&exeCmd("$tool{exe_blastn} -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand\" -num_threads 20 -task blastn -dust no -db $refGff.tmp.ltr.fa -query $inGff.tmp.ltr.fa -out $inGff.tmp.ltr.fa.toTRs.bn6"); 

&tsmsg("[Rec] Searching good LTR with reported TR patterns in $refGff\n"); 
my %keptLTR; 
{
open BN6,'<',"$inGff.tmp.ltr.fa.toTRs.bn6" or &stopErr("[Err] open $inGff.tmp.ltr.fa.toTRs.bn6 . $!\n"); 
my %keptTR; 
while (<BN6>) {
	chomp; s/[^\S\t]+$//; 
	my @ta = split(/\t/, $_); 
	$ta[3] >= $para{minCov}/100 * $ta[12] or $ta[3] >= $para{minCov}/100 * $ta[13] or next; 
	if ( $ta[0] =~ m/^RR(\d+)_(1|2)$/ ) {
		$keptTR{$1}{$2}++; 
	} else {
		&stopErr("[Err] Unkonwn ta[0]=$ta[0]\n"); 
	}
}
close BN6; 
open O,'>',"$inGff.keptLTR.ID" or &stopErr("[Err] Failed to open $inGff.keptLTR.ID . $!\n"); 
for my $rid (sort { $a<=>$b } keys %keptTR) {
	(defined $keptTR{$rid}{1} and defined $keptTR{$rid}{2}) or next; 
	$keptLTR{$rid} ++; 
	print O "RR$rid\tRR${rid}_\n"; 
}
close O; 
my $cnt = scalar(keys %keptLTR); 
&tsmsg("[Rec] Total $cnt LTR kept in $inGff.keptLTR.ID .\n"); 
&exeCmd("perl $tool{pl_deal_table} $inGff.tmp.tab -kSrch_idx $inGff.keptLTR.ID -kSrch_idxCol 0 -kSrch_srcCol 0 > $inGff.keptLTR.tab"); 
}

&tsmsg("[Rec] Filter $inGff to generate $inGff.keptLTR.gff\n"); 
{
# Filter $inGff ("PG1All_v2.scf.fa.gff85") by $inGff.keptLTR.tab 
open IG,'<',"$inGff" or &stopErr("[Err] Failed to open $inGff . $!\n"); 
open OG,'>',"$inGff.keptLTR.gff" or &stopErr("[Err] Failed to open $inGff.keptLTR.gff . $!\n"); 
my $curID = -1; 
while (<IG>) {
	chomp; s/[^\S\t]+$//; 
	if ( m/^\s*#|^\s*$/ ) {
		print OG "$_\n"; 
		next; 
	}
	my @ta = split(/\t/, $_); 
	if ($ta[8] =~ m/^ID=repeat_region(\d+)$/) {
		$curID = $1; 
	} elsif ( $ta[8] =~ m/^ID=LTR_retrotransposon(\d+);/ ) {
		$curID == $1 or &stopErr("[Err] Format strange : $_\n"); 
	} elsif ( $ta[8] =~ m/^Parent=(?:repeat_region|LTR_retrotransposon)(\d+)/ ) {
		$curID == $1 or &stopErr("[Err] Format strange : $_\n"); 
	} else {
		&stopErr("[Err] Unparsed : $_\n"); 
	}
	if ( defined $keptLTR{$curID} ) {
		print OG "$_\n"; 
	}
}
close OG; 
close IG; 
}

&tsmsg("[Rec] All done.\n"); 

### Sub routines.
sub getPath {
	my ($toolR, $cfg_file) = @_;
	open (CF,'<',"$cfg_file") or &stopErr("[Err] file [$cfg_file] $!\n");
	while (<CF>) {
		m/^\s*$/ and next;
		s/[^\S\t]+$//;
		my ($tk, $tv) = split(/\t/, $_);
		while ($tv =~ m/__([^\s_]+)__/) {
			my $pk = $1;
			defined $toolR->{$pk} or &stopErr("[Err] Unknown key [$pk]\n");
			my $pv = $toolR->{$pk};
			$tv =~ s/__${pk}__/$pv/;
		}
		$toolR->{$tk} = $tv;
		&tsmsg("[Msg] Setting $tk=$tv\n");
	}
	close CF;
	return 0;
}#End sub getPath

