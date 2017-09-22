#!/usr/bin/perl -w
use strict; 
use warnings; 
use ConfigSunhh; 
my $cfgs_obj = ConfigSunhh->new(); 
use LogInforSunhh; 
use fileSunhh; 

use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"conf_file:s",     # Required. This file tells the path information of softwares. 
	"in_pref_list:s",  # Required. This file may be output of pipe_gatk.pl -do_step5b . 
	                   # Format : SAMPLE_NAME <tab> READ_GROUP_NAME <tab> LIBRARY_NAME <tab> dataPrefix   <tab> in_fq1 <tab> in_fq2 
	                   #          GS101       <tab> NA              <tab> NA           <tab> GS101_merged <tab> NA     <tab> NA
	                   #          GS102       <tab> GS102_Time2     <tab> GS102_Time2  <tab> GS102_Time2  <tab> Fq1    <tab> Fq2
); 

&chk_para(); 


my %cfg; 
$cfgs_obj->getConfig('cfg_file'=>$opts{'conf_file'}, 'replace'=>1, 'hash_r'=>\%cfg); 

my @fq_infor = @{ &load_prefList($opts{'in_pref_list'}) }; 
# Input files should be ${dataPrefix}_merged_dedup_pipe2.ba[i|m]; 
for my $h1 (@fq_infor) {
	my %th1 = %$h1; 
	print STDOUT "$cfg{'exe_samtools'} mpileup -f $cfg{'ref_fasta'} $cfg{'para_maoPU'} $th1{'pref'}_dedup_pipe2.bam > $th1{'pref'}_dedup_pipe2.bam.pu 2>sMao.err.$th1{'pref'}.s1\n"; 
	print STDOUT "$cfg{'exe_reSeqPrintSample'} $cfg{'ref_fasta'} $th1{'pref'}.1col $th1{'pref'}_dedup_pipe2.bam.pu $th1{'SM'} $cfg{'para_reSeqPrintSample'} 1>sMao.std.$th1{'pref'}.s2 2>sMao.err.$th1{'pref'}.s2\n"; 
}


################################################################################
#  Inner sub-routines. 
################################################################################
sub chk_para {
	for my $t1 (qw/conf_file in_pref_list/) {
		( defined $opts{$t1} and -e $opts{$t1} ) or &usage("[Err] bad para -$t1\n"); 
	}
	return; 
}# chk_para() 


# I would lines heading with '#'
# Return for paired: ([ {'SM'=>SM, 'RG'=>RG, 'LB'=>LB, 'pref'=>pref, 'fq1'=>fq1, 'fq2'=>fq2}, {}, ... ])
# Return for single: ([ {'SM'=>SM, 'RG'=>RG, 'LB'=>LB, 'pref'=>pref, 'fq1'=>fq1, 'fq2'=>''}, {}, ... ])
sub load_prefList {
	my $fn = shift;
	my $fh = &openFH($fn, '<');
	# Format : SAMPLE_NAME <tab> READ_GROUP_NAME <tab> LIBRARY_NAME <tab> dataPrefix <tab> in_fq1 <tab> in_fq2
	my @back;
	while (<$fh>) {
		m/^\s*($|#)/ and next;
		chomp;
		my @ta = split(/\t/, $_);
		my ($sm, $rg, $lb, $pref, $fq1, $fq2) = @ta;
		my %th;
		$opts{'singleFq'} and $fq2 = '';
		$fq2 //= '';
		$th{'SM'} = $sm;
		$th{'RG'} = $rg;
		$th{'LB'} = $lb;
		$th{'pref'} = $pref;
		$th{'fq1'} = $fq1;
		$th{'fq2'} = $fq2;
		push(@back, \%th);
	}
	close($fn);
	return(\@back);
}# load_prefList()

sub usage {
	&tsmsg(@_); 
my $help_txt = <<HH; 
################################################################################
# perl $0   -conf_file pipe_gatk_conf   -in_pref_list step5b_out_pref_list > cmd_list.sMao
#
#
################################################################################
HH
	&LogInforSunhh::usage($help_txt); 
}# usage() 


