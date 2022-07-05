#!/usr/bin/perl
# Input file : 
use strict;
use warnings; 
use LogInforSunhh; 
use Cwd 'abs_path';
use File::Basename;
use Getopt::Long;
my %opts;
GetOptions(\%opts,
	'help!',
	"path_conf:s", # $tool{pathCfg_dir}/path.conf
);
my $help_txt = <<HH;
perl $0 fasfd\n

-help

-path_conf       [filename] check format of path.conf.

HH

!@ARGV and &LogInforSunhh::usage($help_txt);
$opts{'help'} and &LogInforSunhh::usage($help_txt);

# 2.2.	Collection of relatively old LTR retrotransposons

# tools
my %tool; 
{
$tool{pathCfg_dir} = dirname( abs_path($0) );
$tool{pathCfg_file} = "$tool{pathCfg_dir}/path.conf";
defined $opts{'path_conf'} and $tool{'pathCfg_file'} = $opts{'path_conf'};

&getPath(\%tool, $tool{pathCfg_file});

#my $dir1 = "/data/Sunhh/P1_repeat/02.LTR/repAnno_tools"; 
#$tool{pl_ch_gff_to_tab} = "$dir1/ch_gff_to_tab.pl"; 
#$tool{pl_ch_seqID} = "$dir1/ch_seqID.pl"; 
#$tool{pl_filter_tab_byPBSPPT} = "$dir1/filter_tab_byPBSPPT.pl"; 
#$tool{pl_name_from_tab} = "$dir1/name_from_tab.pl"; 
#$tool{pl_filter_flank} = "$dir1/filter_flank.pl"; 
#$tool{pl_filter_RepMsk_out} = "$dir1/filter_RepMsk_out.pl"; 
#$tool{pl_build_Examplar_byFa} = "$dir1/build_Examplar_byFa.pl"; 
#$tool{pl_lis_masked_RepMsk_out} = "$dir1/lis_masked_RepMsk_out.pl"; 
#$tool{pl_get_LTR_wi_Termi} = "$dir1/get_LTR_wi_Termi.pl"; 
#
#$tool{pl_deal_fasta} = "/home/Sunhh/tools/github/NGS_data_processing/deal_fasta.pl"; 
#$tool{pl_deal_table} = "/home/Sunhh/tools/github/NGS_data_processing/deal_table.pl"; 
#
#$tool{exe_RepeatMasker} = "/data/Sunhh/src/Annot/repeatmasker/RepeatMasker/RepeatMasker"; 
#$tool{exe_gt} = '/data/Sunhh/src/Annot/genometools/gt-1.5.3-complete/bin/gt';
#
}

my %input; 
{
$input{eu_tRNA} = $tool{'eu_tRNA'} // '/data/Sunhh/P1_repeat/db/eukaryotic-tRNAs.fa';

$input{refFa} = $tool{'refFa'} // "P1Genom_Gt5h.scf.fa";
$input{refIdx} = $tool{'refIdx'} // 'P1GenomeGt5hScf';

$input{hvt_gff} = "$input{refFa}.gffT85"; 
$input{hvt_outFa} = "$input{refFa}.outT85"; 
$input{hvt_innFa} = "$input{refFa}.outinnerT85"; 
$input{hvt_res} = "$input{refFa}.resultT85"; 
$input{dgt_gff} = "$input{refFa}.gffT85.dgt"; 

$input{ref_hvt_gff} = "$input{refFa}.gffT99"; 
$input{ref_hvt_outFa} = "$input{refFa}.outT99"; 
$input{ref_dbLTR} = $tool{'ref_dbLTR_trim99'} // $tool{'ref_dbLTR'} // "TRIM99_named.lib"; 
}


# Step 2.2.1 Collection of candidate elements with LTRs that are 85% or more in similarity using LTRharvest 
&exeCmd("$tool{exe_gt} suffixerator -db $input{refFa} -indexname $input{refIdx} -tis -suf -lcp -des -ssp -dna"); 
&exeCmd("$tool{exe_gt} ltrharvest -index $input{refIdx} -out $input{hvt_outFa} -outinner $input{hvt_innFa} -gff3 $input{hvt_gff} -minlenltr 70 -maxlenltr 500 -mindistltr 280 -maxdistltr 1500 -mintsd 5 -maxtsd 5 -vic 10  > $input{hvt_res}"); 

# Step 2.2.1.2 Since the terminal sequence motif is not specified, only elements with terminal sequences with patterns that are previously reported are retained.
&exeCmd("perl $tool{pl_get_LTR_wi_Termi} $input{hvt_gff} $input{hvt_outFa} $input{refFa} $input{ref_hvt_gff} $input{ref_hvt_outFa} -path_conf $tool{'pathCfg_file'}"); 

# Step 2.2.2. Using LTRdigest to find elements with PPT (poly purine tract) or PBS (primer binding site)
&exeCmd("$tool{exe_gt} gff3 -sort $input{hvt_gff}.keptLTR.gff > $input{hvt_gff}.use"); 
&exeCmd("$tool{exe_gt} ltrdigest -trnas $input{eu_tRNA} $input{hvt_gff}.use $input{refIdx} > $input{dgt_gff}"); 
## Format results and filter outFa/innFa with PP. 
&exeCmd("perl $tool{pl_ch_gff_to_tab} $input{dgt_gff} 1>dgt.tab"); 
&exeCmd("perl $tool{pl_ch_seqID} $input{hvt_outFa} dgt.tab 1>full_LTR.fa 2>dgt.tab1"); 
&exeCmd("perl $tool{pl_ch_seqID} $input{hvt_innFa} dgt.tab 1>inner.fa 2>dgt.tab2"); 
&exeCmd("mv dgt.tab1 dgt.tab"); 
&exeCmd("perl $tool{pl_filter_tab_byPBSPPT} dgt.tab > dgt.tab.wPP"); # table with PBS / PPT information. 
&exeCmd("perl $tool{pl_name_from_tab} dgt.tab.wPP > dgt.tab.wPPID"); 
## Run repeatmasker to remove sequence masked by LTR99.lib 
&exeCmd("perl $tool{pl_deal_fasta} full_LTR.fa -drawByList -drawWhole -drawList dgt.tab.wPPID -drawLcol 0 > dgt.tab.wPP.full_LTR.fa"); 
&exeCmd("$tool{exe_RepeatMasker} -lib $input{ref_dbLTR} dgt.tab.wPP.full_LTR.fa -nolow -norna -no_is -pa 20 -a 1>stdout.RepMsk_refLTR 2>stderr.RepMsk_refLTR"); 
&exeCmd("perl $tool{pl_deal_fasta} dgt.tab.wPP.full_LTR.fa -attr key:len > dgt.tab.wPP.full_LTR.fa.kl"); 
&exeCmd("perl $tool{pl_lis_masked_RepMsk_out} dgt.tab.wPP.full_LTR.fa.kl dgt.tab.wPP.full_LTR.fa.out > dgt.tab.wPP.full_LTR.fa.kl.toRM"); 
&exeCmd("perl $tool{pl_deal_fasta} dgt.tab.wPP.full_LTR.fa -drawByList -drawWhole -dropMatch -drawLcol 0 -drawList dgt.tab.wPP.full_LTR.fa.kl.toRM > notMaskedByrefLTR.full_LTR.fa"); 

# Step 2.1.3. Further filtering of the candidate elements
&exeCmd("perl $tool{pl_deal_fasta} -listSite \'[nN]+\' notMaskedByrefLTR.full_LTR.fa | awk \' \$5 >= 50 \' | perl -e \' while (<>) { m/^Key/ and next; m/^RR(\\d+)/ or die \"\$_\\n\"; print \"RR\$1_\tRR\$1\\n\"; }  ' > notMaskedByrefLTR.full_LTR.fa.badNlis"); 
&exeCmd("perl $tool{pl_deal_table} dgt.tab.wPP -kSrch_idx dgt.tab.wPP.full_LTR.fa.kl.toRM -kSrch_idxCol 2 -kSrch_srcCol 0 -kSrch_drop | perl $tool{pl_deal_table} -kSrch_idx notMaskedByrefLTR.full_LTR.fa.badNlis -kSrch_idxCol 1 -kSrch_srcCol 0 -kSrch_drop > dgt.tab.wPP.rmrefLTR.filtN"); 
&exeCmd("perl $tool{pl_filter_flank} $input{refFa} dgt.tab.wPP.rmrefLTR.filtN 1>dgt.tab.wPP.rmrefLTR.filtN.chkFlank 2>dgt.tab.wPP.rmrefLTR.filtN.chkFlank.err"); 
&exeCmd("awk '\$1 != \"eleID\" \&\& ( (\$17 >= 0.5 \&\& \$18 >= 0.6) || (\$19 >= 0.5 \&\& \$20 >= 0.6) ) ' dgt.tab.wPP.rmrefLTR.filtN.chkFlank > dgt.tab.wPP.rmrefLTR.filtN.badFlank"); 
&exeCmd("perl $tool{pl_deal_table} dgt.tab.wPP.rmrefLTR.filtN.chkFlank -kSrch_idx dgt.tab.wPP.rmrefLTR.filtN.badFlank -kSrch_idxCol '0-11' -kSrch_srcCol '0-11' -kSrch_drop > dgt.tab.wPP.rmrefLTR.filtN.filtFlank"); 
&exeCmd("perl $tool{pl_name_from_tab} dgt.tab.wPP.rmrefLTR.filtN.filtFlank > dgt.tab.wPP.rmrefLTR.filtN.filtFlankID"); 
&exeCmd("perl $tool{pl_deal_fasta} inner.fa    -drawByList -drawWhole -drawLcol 0 -drawList dgt.tab.wPP.rmrefLTR.filtN.filtFlankID > dgt.tab.wPP.rmrefLTR.filtN.filtFlank.inner.fa"); 
&exeCmd("perl $tool{pl_deal_fasta} full_LTR.fa -drawByList -drawWhole -drawLcol 0 -drawList dgt.tab.wPP.rmrefLTR.filtN.filtFlankID > dgt.tab.wPP.rmrefLTR.filtN.filtFlank.full_LTR.fa"); 
## Here we get all candidate elements. 

# Step 2.1.4. Identify elements with nested insertions
&exeCmd("perl $tool{pl_deal_fasta} $input{refFa} -drawByList -drawList dgt.tab.wPP.rmrefLTR.filtN.filtFlank -drawLcol 15,5,6,3 >  dgt.tab.wPP.rmrefLTR.filtN.filtFlank.ltr.fa"); 
&exeCmd("perl $tool{pl_deal_fasta} $input{refFa} -drawByList -drawList dgt.tab.wPP.rmrefLTR.filtN.filtFlank -drawLcol 15,7,8,3 >>  dgt.tab.wPP.rmrefLTR.filtN.filtFlank.ltr.fa"); 
&exeCmd("$tool{exe_RepeatMasker} -lib dgt.tab.wPP.rmrefLTR.filtN.filtFlank.ltr.fa dgt.tab.wPP.rmrefLTR.filtN.filtFlank.inner.fa -nolow -norna -no_is -pa 20 -a 1>stdout.RepMsk_ltr2Inn 2>stderr.RepMsk_ltr2Inn"); 
&exeCmd("perl $tool{pl_filter_RepMsk_out} dgt.tab.wPP.rmrefLTR.filtN.filtFlank.inner.fa.out > dgt.tab.wPP.rmrefLTR.filtN.filtFlank.inner.fa.out.use"); 
&exeCmd("perl $tool{pl_deal_table} -symbol \'\\s+\' -column 4 dgt.tab.wPP.rmrefLTR.filtN.filtFlank.inner.fa.out.use | perl $tool{pl_deal_table} -UniqColLine 0 | perl -e \' while (<>) { m/^(RR\\d+)_/ or die \"\$_\\n\"; chomp; print \"\$_\\t\$1\\t\$1_\\n\";  }  ' > nested_LTR_list"); 
&exeCmd("perl $tool{pl_deal_fasta} dgt.tab.wPP.rmrefLTR.filtN.filtFlank.inner.fa -drawByList -drawWhole -dropMatch -drawLcol 0 -drawList nested_LTR_list > dgt.tab.wPP.rmrefLTR.filtN.filtFlank.woNest.inner.fa"); 


# Step 2.1.5 Building examplars
## get examplar according to inner without nested
&exeCmd("perl $tool{pl_build_Examplar_byFa} dgt.tab.wPP.rmrefLTR.filtN.filtFlank.woNest.inner.fa 1>stdout.build_examplar_woN 2>stderr.build_examplar_woN"); 
&exeCmd("grep \\> dgt.tab.wPP.rmrefLTR.filtN.filtFlank.woNest.inner.fa.examplars | perl -e 'while (<>) { m/^>(RR\\d+)_/ or die \"\$_\\n\"; print \"\$1_\\t\$1\\n\"; }' > dgt.tab.wPP.rmrefLTR.filtN.filtFlank.woNest.inner.fa.examplars.LTRID"); 
&exeCmd("perl $tool{pl_deal_fasta} full_LTR.fa -drawByList -drawList dgt.tab.wPP.rmrefLTR.filtN.filtFlank.woNest.inner.fa.examplars.LTRID -drawWhole -drawIDmatch -drawLcol 0 > woNest.fullLTR.examplars"); 

## get examplar after removing inner-examplars
&exeCmd("$tool{exe_RepeatMasker} -lib woNest.fullLTR.examplars dgt.tab.wPP.rmrefLTR.filtN.filtFlank.full_LTR.fa -nolow -norna -no_is -pa 20 -a 1>stdout.RepMsk_woN2all 2>stderr.RepMsk_woN2all"); 
&exeCmd("perl $tool{pl_deal_fasta} dgt.tab.wPP.rmrefLTR.filtN.filtFlank.full_LTR.fa -attr key:len > dgt.tab.wPP.rmrefLTR.filtN.filtFlank.full_LTR.fa.kl"); 
&exeCmd("perl $tool{pl_lis_masked_RepMsk_out} dgt.tab.wPP.rmrefLTR.filtN.filtFlank.full_LTR.fa.kl dgt.tab.wPP.rmrefLTR.filtN.filtFlank.full_LTR.fa.out > dgt.tab.wPP.rmrefLTR.filtN.filtFlank.full_LTR.fa.kl.toRM"); 
&exeCmd("perl $tool{pl_deal_fasta} dgt.tab.wPP.rmrefLTR.filtN.filtFlank.full_LTR.fa -drawByList -drawWhole -dropMatch -drawLcol 0 -drawList dgt.tab.wPP.rmrefLTR.filtN.filtFlank.full_LTR.fa.kl.toRM > notMasked.full_LTR.fa"); 
&exeCmd("perl $tool{pl_build_Examplar_byFa} notMasked.full_LTR.fa 1>stdout.build_examplar_notMsk 2>stderr.build_examplar_notMsk"); 

## Combine both examplars. 
&exeCmd("cat woNest.fullLTR.examplars notMasked.full_LTR.fa.examplars | perl $tool{pl_deal_fasta} -frag_head -frag_width 80 -frag 0-0 | perl $tool{pl_deal_fasta} -chopKey \':\\d+\-\\d+\$' > TRIM85.lib"); 

### Sub routines.
sub getPath {
	my ($toolR, $cfg_file) = @_;
	open (CF,'<',"$cfg_file") or &stopErr("[Err] file [$cfg_file] $!\n");
	while (<CF>) {
		m/^\s*($|#)/ and next;
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

