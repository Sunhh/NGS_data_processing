#!/usr/bin/perl
# Reference : http://bioinf.uni-greifswald.de/bioinf/wiki/pmwiki.php?n=IncorporatingRNAseq.Tophat
# Ref.02 : http://bioinf.uni-greifswald.de/bioinf/wiki/pmwiki.php?n=Augustus.ParallelPred
# Ref.03 : http://bioinf.uni-greifswald.de/bioinf/wiki/pmwiki.php?n=Augustus.UTRTraining

use strict; 
use warnings; 
use fileSunhh; 
use LogInforSunhh; 
use File::Path qw(make_path remove_tree); 
use IPC::Open3; 
use Symbol; 
use Parallel::ForkManager; 
use Cwd 'abs_path'; 

use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"stepLis:s", "continue!", 
	"cpuN:i", "path_samtools:s", "path_bwt2:s", "path_bwt2build:s", "path_bamtools:s", 
	"samtools_cpu:i", "samtools_mem:s", 
	"inDir:s", "inRdLis:s", "faRd!", 
	"genomFasMsk:s", "genomFasRaw:s", 
	"oPref:s", 
	"path_tophat2:s", "mode:s", "dbBwt2:s", "para_th:s", "para_cl:s", "dir_cuff:s", "exex_flank:i", 
	"dir_aug:s", "species:s", "chunk_overlap:i", "chunk_size:i", "min_ctg_len:i", "aug_utr!"
); 
$opts{'stepLis'} = $opts{'stepLis'} // ''; 
$opts{'cpuN'} = $opts{'cpuN'} // 30; 
$opts{'path_bwt2'} = $opts{'path_bwt2'} // 'bowtie2'; 
$opts{'path_bwt2build'} = $opts{'path_bwt2build'} // 'bowtie2-build'; 
$opts{'path_bamtools'} = $opts{'path_bamtools'} // 'bamtools'; 
$opts{'inDir'} = $opts{'inDir'} // ''; 
$opts{'inRdLis'} = $opts{'inRdLis'} // 'cult_inLis' ; 
$opts{'oPref'} = $opts{'oPref'} // 'oPref'; 
$opts{'path_tophat2'} = $opts{'path_tophat2'} // 'tophat2'; 
$opts{'dir_cuff'} = $opts{'dir_cuff'} // "/data/Sunhh/src/Annot/Cufflinks/cufflinks-2.2.1.Linux_x86_64"; 
$opts{'mode'} = $opts{'mode'} // 'single'; 
$opts{'dbBwt2'} = $opts{'dbBwt2'} // 'db/bwt2_db'; 
$opts{'para_th'} = $opts{'para_th'} // "--library-type=fr-firststrand --read-mismatches 1 --splice-mismatches 0 --min-intron-length 30"; 
$opts{'para_cl'} = $opts{'para_cl'} // "--min-intron-length 30 --min-frags-per-transfrag 5"; 
$opts{'path_samtools'} = $opts{'path_samtools'} // 'samtools'; 
$opts{'dir_aug'} = $opts{'dir_aug'} // "/data/Sunhh/src/Annot/maker/maker/exe/augustus"; 
$opts{'species'} = $opts{'species'} // 'arabidopsis'; 
$opts{'genomFasMsk'} = $opts{'genomFasMsk'} // 'genome.msk.fa'; 
$opts{'genomFasRaw'} = $opts{'genomFasRaw'} // 'genome.raw.fa'; 
$opts{'chunk_overlap'} = $opts{'chunk_overlap'} // 50000; 
$opts{'chunk_size'} = $opts{'chunk_size'} // 1000000; 
$opts{'min_ctg_len'} = $opts{'min_ctg_len'} // 200; 
$opts{'samtools_cpu'} = $opts{'samtools_cpu'} // 1; 
$opts{'samtools_mem'} = $opts{'samtools_mem'} // ''; 
$opts{'exex_flank'} = $opts{'exex_flank'} // 150; 

sub usage {
	print STDOUT <<HH;
################################################################################
#   perl $0 -stepLis 123 [-continue]
#
#   -stepLis      [123] 1-Run tophat (step1/). 
#                       2-Run the first run of rna_augustus (step1/)
#                       3-Run the second run of rna_augustus (step2/)
#   -continue     [boolean] Will skip finished steps. 
#
#
### Parameters for basic settings. 
#   -cpuN             [30] threads number used. 
#
#   -path_samtools    [samtools] path to samtools bin. 
#   -samtools_cpu     [1] number of \-\@ input for samtools. 
#   -samtools_mem     [\'\'] No memory size assigned. Use this carefully because some samtools version doesn't support it. 
#
#   -path_bwt2        [bowtie2] path to bowtie2 bin. 
#   -path_bwt2build   [bowtie2-build] path to bowtie2-build bin. 
#
#   -path_bamtools    [bamtools] path to bamtools bin. 
#
#   -path_tophat2     [tophat2]
#   -dir_cuff         [/data/Sunhh/src/Annot/Cufflinks/cufflinks-2.2.1.Linux_x86_64]
#
#   -dir_aug          [/data/Sunhh/src/Annot/maker/maker/exe/augustus]
#   -species          [arabidopsis]
#
#
### Input personal parameters. 
#   -inDir            [\'\'] path of directory added to file name from rnaseq list. 
#   -inRdLis          [cult_inLis] rnaseq list. 
#   -mode             [single] rnaseq data construction. 
#   -faRd             [boolean] rnaseq seq format is fasta. 
#   -genomFasMsk      [genome.msk.fa]
#   -genomFasRaw      [genome.raw.fa]
#   -oPref            [oPref]
#   -dbBwt2           [db/bwt2_db] 
#
#   -para_th          [--library-type=fr-firststrand --read-mismatches 1 --splice-mismatches 0 --min-intron-length 30]
#   -para_cl          [--min-intron-length 30 --min-frags-per-transfrag 5]
#   -exex_flank       [150] Should be no less than read length. 
#   -chunk_overlap    [50000]
#   -chunk_size       [1000000]
#   -min_ctg_len      [200]
#   -aug_utr          [boolean] predict with \'UTR=on\' in augustus. 
################################################################################
HH
	exit 1; 
}

$opts{'help'} and &usage(); 
( defined $opts{'stepLis'} and $opts{'stepLis'} ne '' ) or &usage(); 



my %need_step; 
for my $tk (split(//, $opts{'stepLis'})) {
	$need_step{$tk} = 1; 
}


my @rdFiles; 
{
	my $inFh = &openFH($opts{'inRdLis'}, "<"); 
	while (<$inFh>) {
		chomp; 
		my @ta = split(/\s+/, $_); 
		if ( $ta[0] =~ m!^/! ) {
			push(@rdFiles, $ta[0]); 
		} elsif ( defined $opts{'inDir'} and $opts{'inDir'} ne '' ) {
			# push(@rdFiles, Cwd::abs_path("$opts{'inDir'}/$ta[0]")); 
			push(@rdFiles, "$opts{'inDir'}/$ta[0]"); 
		} else {
			push(@rdFiles, $ta[0]); 
		}
	}
	close ($inFh); 
}

my $step1_oDir = "step1/$opts{'oPref'}_thout"; 
my $step1_oBam = "${step1_oDir}/accepted_hits.bam"; 
if ( defined $need_step{1} ) {
	# Step1 : Run tophat2 to align all reads to reference genome. 
	#  Out information from step1 : 
	#   $step1_oBam = "${step1_oDir}/accepted_hits.bam"; 
	&tsmsg("[Rec] Start stepLis=1\n"); 
	my $tmp_rdLis = join(',', @rdFiles); 
	&if_redo("step1/$opts{'oPref'}_thout/", $step1_oBam) and &exeCmd("$opts{'path_tophat2'} $opts{'para_th'} -o step1/$opts{'oPref'}_thout $opts{'dbBwt2'} $tmp_rdLis "); 
	undef($tmp_rdLis); 
	&tsmsg("[Rec] Finish stepLis=1\n"); 
}# End if ( need_step 1 ) 

my $step2_oAug = 'step1/aug1.out'; 
my $step2_hint = 'step1/hints.gff'; 
if ( defined $need_step{2} ) {
	# Step2 : The first run of augustus. 
	#  
	## samtools sort -@ 10 -m 10G AllToPG1_thout/accepted_hits.bam step1/both.ssf
	## $opts{'dir_aug'}/auxprogs/bam2hints/bam2hints
	## 
	&tsmsg("[Rec] Start stepLis=2\n"); 
	-d 'step1' or &stopErr("[Err] No step1 directory found!\n"); 
	# samtools sort output_directory/accepted_hits.sf.bam both.ssf
	my $smT_para = ''; 
	$opts{'samtools_cpu'} > 1 and $smT_para .= " -\@ $opts{'samtools_cpu'}"; 
	$opts{'samtools_mem'} ne '' and $smT_para .= " -m $opts{'samtools_mem'}"; 
	&if_redo( "step1/both.ssf.bam" ) and &exeCmd("$opts{'path_samtools'} sort $smT_para $step1_oBam step1/both.ssf"); 
	$smT_para = ''; 
	# bam2hints --intronsonly --in=both.ssf.bam --out=hints.gff
	&if_redo( $step2_hint ) and &exeCmd("$opts{'dir_aug'}/auxprogs/bam2hints/bam2hints --intronsonly --in=step1/both.ssf.bam --out=$step2_hint"); 
	# Edit extrinsic.cfg file. 
	#### The exonpart malus of .992 means a weak penalty factor for every predicted exonic base that is not supported by any exonpart hints. The exonpart bonus for hints of source W of 1.005 mean that gene structures get this bonus factor for every exonpart hint of multiplicity 1 that is completely included in an exon. Introns that are not supported by any intron hint are penalized by .34, and introns that are supported by RNA-Seq hints are rewarded by a factor of 100,000. The 0.985 are local malus factors. The concept of a local malus was recently introduced to account for different levels of missing information, that becomes more important with with higher general coverage of RNA-Seq. Like the (normal) malus, the local malus also applies to exonic bases that are unsupported by hints. In contract to the normal malus, the local malus only applies to exons, that are well-supported at some region and not supported at another region of the same exon. The following picture illustrates a typical case where the local malus applies.
	## $opts{'dir_aug'}/config/extrinsic/extrinsic.M.RM.E.W.cfg
	#### The optimal values will depend on various parameters, including 
	####  1) Amount of RNA-Seq and percentage of genes expressed in the library. 
	####     The more transcripts are expressed and covered by RNA-Seq the stronger the malus can be. 
	####     In the extreme (hypothetical) case, all exonic bases are covered with RNA-Seq and 
	####     unsupported predicted exonic regions can be punished hard.
	####  2) The stringency of filtering. 
	####     If, for some other reason, you have chosen to include more than one hit of a read or 
	####     you allow best alignments also when they are ambiguous, then you may get a significant 
	####     false positive rate of hints and you may have to reduce bonus factors significantly.
	####  3) Alignment method. 
	$opts{'cfgEx_aug'} = $opts{'cfg_aug'} // "$opts{'dir_aug'}/config/extrinsic/extrinsic.M.RM.E.W.cfg"; 
	&if_redo( "step1/extrinsic.cfg" ) and &exeCmd("cp $opts{'cfgEx_aug'} step1/extrinsic.cfg"); 
	# augustus --species=yourSpecies 
	#  --extrinsicCfgFile=extrinsic.cfg 
	#  --alternatives-from-evidence=true 
	#  --hintsfile=hints.gff 
	#  --allow_hinted_splicesites=atac 
	#  --introns=on --genemodel=complete 
	#  genome.fa 
	#  > aug1.out
	$step2_oAug = &multiRunAug (
	  'cpuN'      => $opts{'cpuN'} , 
	  'dir_aug' => $opts{'dir_aug'}, 
	  'path_aug'  => "$opts{'dir_aug'}/bin/augustus",
	  'species'   => $opts{'species'}, 
	  'genomFas'  => $opts{'genomFasRaw'}, 
	  'outAug'    => ( ( defined $need_step{3} ) ? 'step1/aug1.out' : 'step1/aug.out' ), 
	  'para'      => " --extrinsicCfgFile=step1/extrinsic.cfg --alternatives-from-evidence=true "
	    . ( ( defined $opts{'aug_utr'} ) ? ' --UTR=on ' : '' ) 
	    . " --hintsfile=$step2_hint --allow_hinted_splicesites=atac "
		. ( ( defined $need_step{3} ) ? '' : " --introns=on --genemodel=complete " ) , 
	  'min_ctg_len' => $opts{'min_ctg_len'}, 
	  'file_chr_list' => "step1/chr.lst", 
	  'chunk_overlap' => $opts{'chunk_overlap'}, 
	  'chunk_size' => $opts{'chunk_size'}, 
	  'dir_SoutAug' => 'step1/SoutAugDir', 
	  'joblist' => 'step1/jobs.lst'
	); 
	&tsmsg("[Rec] Finish stepLis=2\n"); 
}# End Step2 : aug_step1 

if ( defined $need_step{3} ) {
	# Step3 : The second run of augustus. 
	#  Input : $step2_oAug , $step2_hint 
	&tsmsg("[Rec] Start stepLis=3\n"); 
	-d 'step2' or mkdir('step2', 0755); 
	# cat aug1.out | tee aug.prelim.gff | grep -P "\tintron\t" > aug1.introns.gff
	&if_redo("step1/aug1.introns.gff") and &exeCmd("cat $step2_oAug | grep -P \"\\tintron\\t\" > step1/aug1.introns.gff"); 
	# cat hints.gff aug1.introns.gff | 
	#  perl -ne '@array = split(/\t/, $_);print "$array[0]:$array[3]-$array[4]\n";'| sort -u > introns.lst
	&if_redo("step1/introns.lst") and &exeCmd("cat $step2_hint step1/aug1.introns.gff " 
	  . " | perl -ne '\@array=split(/\\t/, \$_); print \"\$array[0]:\$array[3]-\$array[4]\\n\"; ' " 
	  . " | sort -u > step1/introns.lst "
	); 
	# intron2exex.pl --introns=introns.lst --seq=genome.masked.fa --exex=exex.fa --map=map.psl
	&if_redo("step1/exex.fa", "step1/map.psl") and &exeCmd("perl $opts{'dir_aug'}/scripts/intron2exex.pl --flank=$opts{'exex_flank'} --introns=step1/introns.lst --seq=$opts{'genomFasMsk'} --exex=step1/exex.fa --map=step1/map.psl"); 
	# bowtie2-build exex.fa your_species_exex1
	my @chk_bt2_files; 
	for my $suff (qw/.1.bt2 .2.bt2 .3.bt2 .4.bt2 .rev.1.bt2 .rev.2.bt2/) {
		push(@chk_bt2_files, "step1/exex.fa$suff"); 
	}
	&if_redo(@chk_bt2_files) and &exeCmd("$opts{'path_bwt2build'} step1/exex.fa step1/exex.fa"); 
	# bowtie2 -x your_species_exex1 -U rnaseq.fastq -S bowtie.sam
	my $smT_para = ''; 
	$opts{'samtools_cpu'} > 1 and $smT_para .= " -\@ $opts{'samtools_cpu'}"; 
	my $tmp_rdLis = join(',', @rdFiles); 
	&if_redo("step1/bowtie.F.sam") and &exeCmd("$opts{'path_bwt2'} -p $opts{'cpuN'} " . ( ($opts{'faRd'}) ? ' -f ' : ' -q ' ) . " -x step1/exex.fa -U $tmp_rdLis " . " | $opts{path_samtools} view $smT_para -S -F 4 -o step1/bowtie.F.sam -"); 
	undef($tmp_rdLis); 
	# samMap.pl bowtie.F.sam map.psl > bowtie.global.sam
	&if_redo("step1/bowtie.global.sam") and &exeCmd("$opts{dir_aug}/scripts/samMap.pl step1/bowtie.F.sam step1/map.psl > step1/bowtie.global.sam"); 
	# discard intron containing alignments from the original bam file
	# bamtools filter -in accepted_hits.bam -out output_directory/accepted_hits.noN.bam 
	#  -script operation_N_filter.txt
	&if_redo("${step1_oDir}/accepted_hits.noN.bam") and &exeCmd("$opts{'path_bamtools'} filter -in $step1_oBam -out ${step1_oDir}/accepted_hits.noN.bam -script $opts{'dir_aug'}/auxprogs/auxBamFilters/operation_N_filter.txt"); 
	# create a bam file with header from the bowtie.global.sam file
	# cat header.txt bowtie.global.sam > bowtie.global.h.sam
	# samtools view -bS -o bowtie.global.h.bam bowtie.global.h.sam
	&if_redo("step1/header.txt") and &exeCmd("$opts{'path_samtools'} view $smT_para -H $step1_oBam > step1/header.txt"); 
	&if_redo("step1/bowtie.global.h.bam") and &exeCmd("cat step1/header.txt step1/bowtie.global.sam " 
	  ." | $opts{'path_samtools'} view $smT_para -bS -o step1/bowtie.global.h.bam -"
	); 
	# join bam files
	# bamtools merge -in bowtie.global.h.bam -in output_directory/accepted_hits.noN.bam -out both.bam
	&if_redo("step1/both.bam") and &exeCmd("$opts{'path_bamtools'} merge -in step1/bowtie.global.h.bam -in ${step1_oDir}/accepted_hits.noN.bam -out step1/both.bam"); 
	# samtools sort -n both.bam both.s
	$smT_para = ''; 
	$opts{'samtools_cpu'} > 1 and $smT_para .= " -\@ $opts{'samtools_cpu'}"; 
	$opts{'samtools_mem'} ne '' and $smT_para .= " -m $opts{'samtools_mem'}"; 
	&if_redo("step2/both.s.bam") and &exeCmd("$opts{'path_samtools'} sort $smT_para -n step1/both.bam step2/both.s"); 
	# filterBam --uniq --paired --in both.s.bam --out both.sf.bam
	#  --paired should be omitted when using single-end reads. 
	my $tmp_pair = ($opts{'mode'} eq 'single') ? '' : '--paired' ; 
	&if_redo("step2/both.sf.bam") and &exeCmd("$opts{'dir_aug'}/auxprogs/filterBam/bin/filterBam --uniq --in step2/both.s.bam --out step2/both.sf.bam"); 
	# samtools sort both.sf.bam both.ssf
	# bam2hints --intronsonly --in=both.ssf.bam --out=hints.2.gff
	&if_redo("step2/both.ssf.bam") and &exeCmd("$opts{'path_samtools'} sort $smT_para step2/both.sf.bam step2/both.ssf"); 
	&if_redo("step2/hints.2.gff.ip") and &exeCmd("$opts{'dir_aug'}/auxprogs/bam2hints/bam2hints --intronsonly --in=step2/both.ssf.bam --out=step2/hints.2.gff.ip"); 
	$smT_para = ''; 
	&if_redo("step2/both.ssf.wig") and &exeCmd("$opts{'dir_aug'}/auxprogs/bam2wig/bam2wig step2/both.ssf.bam > step2/both.ssf.wig"); 
	&if_redo("step2/hints.ep.gff") and &exeCmd("cat step2/both.ssf.wig | perl $opts{'dir_aug'}/scripts/wig2hints.pl --width=10 --margin=10 --minthresh=2 --minscore=4 --prune=0.1 --src=W --type=ep --UCSC=unstranded.track --radius=4.5 --pri=4 --strand=\".\" > step2/hints.ep.gff"); 
	&if_redo("step2/hints.2.gff") and &exeCmd("cat step2/hints.ep.gff step2/hints.2.gff.ip > step2/hints.2.gff"); 

	### Maybe include "exonparthints" and "UTR" models to improve the prediciton. 
	### Later. 
	&if_redo("step2/extrinsic.cfg") and &exeCmd("cp $opts{'cfgEx_aug'} step2/extrinsic.cfg"); 
	
	# augustus --species=yourSpecies --extrinsicCfgFile=extrinsic.cfg 
	#  --alternatives-from-evidence=true 
	#  --hintsfile=hints.2.gff 
	#  --allow_hinted_splicesites=atac 
	#  genome.fa > aug2.out
	my $step3_oAug = &multiRunAug (
	  'cpuN'      => $opts{'cpuN'} , 
	  'dir_aug'   => $opts{'dir_aug'}, 
	  'path_aug'  => "$opts{'dir_aug'}/bin/augustus",
	  'species'   => $opts{'species'}, 
	  'genomFas'  => $opts{'genomFasRaw'}, 
	  'outAug'    => 'step2/aug2.out' , 
	  'para'      => " --extrinsicCfgFile=step2/extrinsic.cfg --alternatives-from-evidence=true "
	    . " --hintsfile=step2/hints.2.gff --allow_hinted_splicesites=atac " , 
	  'min_ctg_len' => $opts{'min_ctg_len'}, 
	  'file_chr_list' => "step2/chr.lst", 
	  'chunk_overlap' => $opts{'chunk_overlap'}, 
	  'chunk_size' => $opts{'chunk_size'}, 
	  'dir_SoutAug' => 'step2/SoutAugDir', 
	  'joblist' => 'step2/jobs.lst'
	); 
	&tsmsg("[Rec] Finish stepLis=3\n"); 
}

##########################################################################################
#  Sub-routines 
##########################################################################################


sub multiRunAug {
	my %parm = @_; 
	$parm{'cpuN'} = int($parm{'cpuN'}) // 1; 
	$parm{'path_aug'} = $parm{'path_aug'} // 'augustus'; 
	$parm{'species'} = $parm{'species'} // 'generic'; 
	( defined $parm{'genomFas'} and -f $parm{'genomFas'} ) or &stopErr("[Err] Genome fa file [$parm{'genomFas'}] not found.\n"); 
	$parm{'outAug'} = $parm{'outAug'} // 'aug.out'; 
	$parm{'para'} = $parm{'para'} // ' --extrinsicCfgFile=step1/extrinsic.cfg --alternatives-from-evidence=true --hintsfile=step1/hints.gff --allow_hinted_splicesites=atac '; 
	$parm{'min_ctg_len'} = $parm{'min_ctg_len'} // 200; 
	&if_redo($parm{'outAug'}) or return $parm{'outAug'}; 
	if ( $parm{'cpuN'} == 1 ) {
		&if_redo( $parm{'outAug'} ) and &exeCmd( "$parm{'path_aug'} --species=$parm{'species'} $parm{'para'} $parm{'genomFas'} > $parm{'outAug'}" ); 
		return $parm{'outAug'}; 
	}
	
	$parm{'hr_chr_list'}   = $parm{'hr_chr_list'} // &splitMfasta(%parm); 
	$parm{'file_chr_list'} = &write_chr_list(%parm); 
	
	#scripts/createAugustusJoblist.pl 
    #  --sequences=chr.lst 
	#  --wrap="#" ## This is specific for SGE (qsub) jobs. Should be "" on single-node server. 
	#  --overlap=20000 --chunksize=1252500 
	#  --outputdir=PG1All_v2_Scf_unmask_augDir/ 
	#  --joblist=jobs.lst 
	#    --jobprefix=jobCmds/myPrefix_ ## Should not assign jobprefix for single-node server. 
	#  --command "augustus  --species=spinachV2Scf_Annot01 
	#    --extrinsicCfgFile=/data/feizj/spinach/Annot/rna_augustus/01.byTopHat2_wiSample75/ParallelPred/extrinsic.cfg 
	#    --alternatives-from-evidence=true --allow_hinted_splicesites=atac 
	#    --introns=on --genemodel=complete"
	$parm{'chunk_overlap'} = $parm{'chunk_overlap'} // 50000; 
	$parm{'chunk_size'} = $parm{'chunk_size'} // 1000000; 
	$parm{'dir_SoutAug'} = $parm{'dir_SoutAug'} // 'SoutAugDir'; 
	if ( -d $parm{'dir_SoutAug'} ) {
		&if_redo( $parm{'dir_SoutAug'} ) and File::Path::remove_tree( $parm{'dir_SoutAug'}, {'safe'=>1, 'keep_root'=>1} ); 
	} else {
		mkdir( $parm{'dir_SoutAug'} ); 
	}
	
	$parm{'joblist'} = $parm{'joblist'} // 'jobs.lst'; 
	$parm{'joblist'} =~ s/(^\s+|\s+$)//g; 
	&if_redo( $parm{'joblist'} ) and &exeCmd( "$parm{'dir_aug'}/scripts/createAugustusJoblist.pl "
	  . " --sequences=$parm{'file_chr_list'}" 
	  . " --overlap=$parm{'chunk_overlap'} --chunksize=$parm{'chunk_size'} "
	  . " --outputdir=$parm{'dir_SoutAug'}/ "
	  . " --joblist=$parm{'joblist'} " 
	  . " --command \" "
	    . " $parm{'path_aug'} --species=$parm{'species'} $parm{'para'} "
	  . " \" "
	); 
	
	
	my $jobFh = &openFH("$parm{'joblist'}", "<"); 
	my @joblist = <$jobFh>; 
	close ($jobFh); 
	chomp(@joblist); 
	&if_redo( "$parm{'joblist'}.ok" ) or @joblist = (); 
	my $MAX_PROCESSES = $parm{'cpuN'} ; # Sometimes $parm{'cpuN'} - 1 may be better. 
	my $pm = new Parallel::ForkManager($MAX_PROCESSES); 
	foreach my $jobcmd (@joblist) {
		if ( -f "Nproc" ) {
			my $pFh = &openFH("Nproc", '<'); 
			my $new_maxP = <$pFh>; 
			close ($pFh); 
			chomp($new_maxP); 
			$new_maxP = int($new_maxP); 
			if ( $new_maxP > 0 and $new_maxP != $MAX_PROCESSES ) {
				$pm->set_max_procs($new_maxP); 
				&tsmsg("[Rec] Chaning MAX_PROCESSES from $MAX_PROCESSES to $new_maxP\n"); 
				$MAX_PROCESSES = $new_maxP; 
			}
		}
		
		my $pid = $pm->start and next; 
		&exeCmd($jobcmd); 
		$pm->finish; 
	}#End foreach ( parallel joblist )
	$pm->wait_all_children; 
	&if_redo( "$parm{'joblist'}.ok" ) and &exeCmd("echo \'\' > $parm{'joblist'}.ok"); 
	
	if ( &if_redo( "$parm{'outAug'}" ) ) { 
		opendir(DD, "$parm{'dir_SoutAug'}"); 
		my $rawAugOut = "$parm{'outAug'}.rawAugPred"; 
		unlink($rawAugOut); 
		for my $fn (readdir(DD)) {
			$fn =~ m/\.gff$/ or next; 
			&exeCmd("cat $parm{'dir_SoutAug'}/$fn >> $rawAugOut"); 
		}
		closedir(DD); 
		&exeCmd( "cat $rawAugOut | perl $parm{'dir_aug'}/scripts/join_aug_pred.pl > $parm{'outAug'}" ); 
		# &exeCmd("cat $parm{'dir_SoutAug'}/*.gff | perl $parm{'dir_aug'}/scripts/join_aug_pred.pl > $parm{'outAug'}"); 
	}

	return $parm{'outAug'}; 
}# End sub multiRunAug() 

sub write_chr_list {
	my %parm = @_; 
	$parm{'file_chr_list'} = $parm{'file_chr_list'} // 'chr.lst'; 
	&if_redo( $parm{'file_chr_list'} ) or return $parm{'file_chr_list'} ; 
	&tsmsg("[Rec] Writing sequence list [$parm{'file_chr_list'}]\n"); 
	
	my $outFh = &openFH( $parm{'file_chr_list'}, '>' ); 
	for my $ar ( @{$parm{'hr_chr_list'}{'chr_list'}} ) {
		print {$outFh} join("\t", $ar->[0], 1, $ar->[1])."\n"; 
	}
	close ($outFh); 
	return $parm{'file_chr_list'}; 
}

sub splitMfasta {
	my %parm = @_; 
	my %backH; 
	&tsmsg("[Rec] Splitting $parm{'genomFas'} into single-sequence fasta files.\n"); 
	
	$parm{'min_ctg_len'} = $parm{'min_ctg_len'} // 0; 
	
	$parm{'dir_Sfasta'} = $parm{'dir_Sfasta'} // 'Sfasta.1'; 
	if ( -d $parm{'dir_Sfasta'} ) {
		File::Path::remove_tree( $parm{'dir_Sfasta'}, {'safe'=>1, 'keep_root'=>1} ); 
	} else {
		mkdir( $parm{'dir_Sfasta'} ); 
	}
	$backH{'dir_Sfasta'} = $parm{'dir_Sfasta'}; 
	
	my $faFh = &openFH($parm{'genomFas'}, '<'); 
	for ( my ($relHR, $get) = &get_fasta_seq($faFh); defined $relHR; ($relHR, $get) = &get_fasta_seq($faFh) ) {
		$relHR->{'seq'} =~ s/\s//g; 
		my $ll = length($relHR->{'seq'}); 
		$ll >= $parm{'min_ctg_len'} or next; 
		
		$relHR->{'seq'} =~ s!(.{100})!$1\n!g; chomp($relHR->{'seq'}); 
		my $singleFasFile = "$parm{'dir_Sfasta'}/$relHR->{'key'}.fa"; 
		if ( &if_redo( $singleFasFile ) ) {
			my $outFh = &openFH( $singleFasFile, '>' ); 
			print {$outFh} ">$relHR->{'key'}\n$relHR->{'seq'}\n"; 
			close $outFh; 
		}
		push( @{$backH{'chr_list'}}, [$singleFasFile, $ll] ); 
	}
	
	return \%backH; 
}

# input a fasta file's handle and a signal show whether it should has a head line; 
# return (undef(),undef()) if input is not enough data. 
# return (\% , has_get_next) 
# 2007-12-14 13:02:52 about \%, {seq, key, head}
# 2013-11-01 13:02:52 about \%, {seq, key, head, definition}
sub get_fasta_seq {
	my $fh = shift;
	my $has_head = shift; 
	my $has_get = 0; # 检测是否占用了下一条序列的">"号; 
	my %backH; 
	( defined $has_head and $has_head =~ /^0+$/ ) or $has_head = 1; 
	ref($fh) eq 'GLOB' or ref($fh) eq '' or die "Wrong input!\n"; 
	my %back; 
	if ( $has_head == 1 ) 
	{
		defined ( $backH{head} = readline($fh) ) or return (undef(),undef()); 
		$backH{head} =~ s/^>//g; chomp $backH{head}; 
		$backH{key} = (split(/\s+/,$backH{head}))[0]; 
		( $backH{definition} = $backH{head} ) =~ s/^(\S+)//; 
	}
	my $r = $/; local $/ = "$r>"; 
	defined ( $backH{seq} = readline($fh) ) or do { warn "[Err]The last sequence [$backH{head}] is empty, and it is not calculated!\n"; return (undef(),undef()); } ; 
	chomp $backH{seq} > length($r) and $has_get = 1; 
	local $/ = $r; chomp $backH{seq}; 
	# check if this sequence is a NULL one. 2008-1-4 13:01:22 
	# print "head=+$backH{head}+\nseq=+$backH{seq}+\n"; 
	while ($backH{seq} =~ s/^>//gs) { 
		warn "[Err]Sequence [$backH{head}] is empty, and it is not calculated!\n"; 
		if ($backH{seq} =~ s/^([^$r]+)(?:$r|$)//s) {
			$backH{head} = $1; $backH{key} = (split(/\s+/, $backH{head}))[0]; 
			( $backH{definition} = $backH{head} ) =~ s/^(\S+)//; 
		}
	}
	# check if this sequence is a NULL one. 2008-1-4 13:01:25 
	return (\%backH, $has_get); 
}# end sub get_fasta_seq


sub if_redo {
	my @files = @_; 
	my $is_redo = 0; 
	if ( $opts{'continue'} ) {
		for my $fn ( @files ) {
			-e $fn or do { $is_redo = 1; last; }; 
		}
	} else {
		$is_redo = 1; 
	}
	return $is_redo; 
}# sub if_redo() 

