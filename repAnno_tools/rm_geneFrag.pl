#!/usr/bin/perl 
use strict; 
use warnings; 
use Getopt::Long; 
use fileSunhh; 
use LogInforSunhh; 

my %opts;
GetOptions(\%opts,
	"rawLibLis:s", 
	"rawLib:s", 
	"dbProt:s", 
	"cpuN:i", 
	"evalue:f", 
	"pl_ProtExcluder:s", "pl_dealFa:s", 
"help!", 
); 

$opts{'help'} and &usage(); 
defined $opts{'rawLibLis'} or defined $opts{'rawLib'} or &usage(); 

$opts{'pl_dealFa'} = $opts{'pl_dealFa'} // `echo \$HOME/tools/github/NGS_data_processing/deal_fasta.pl` ;
chomp($opts{'pl_dealFa'}); 
$opts{'pl_ProtExcluder'} = $opts{'pl_ProtExcluder'} // `echo \$HOME/tools/github/NGS_data_processing/repAnno_tools/ProtExcluder1.1/ProtExcluder.pl`; 
chomp($opts{'pl_ProtExcluder'}); 

$opts{'evalue'} = $opts{'evalue'} // 1e-5; 
$opts{'cpuN'} = $opts{'cpuN'} // 1; 
$opts{'dbProt'} = $opts{'dbProt'} // 'uniprot_sprot_plants_rmTransProt.fa'; 

sub usage {
	print STDOUT <<HH;
################################################################################
# perl $0 -rawLibLis inLib.list -pl_ProtExcluder ProtExcluder.pl -dbProt uniprot_sprot_plants_rmTransProt.fa -pl_dealFa deal_fasta.pl 
# 
#   -rawLib
#   -cpuN
#   -evalue
#
#   -help 
################################################################################
HH
exit 1; 
}

my $tmp_file = 'tmp_inLib.fa'; 
my @inLibs; 
defined $opts{'rawLib'} and push(@inLibs, $opts{'rawLib'}); 
defined $opts{'rawLib'} and $opts{'rawLib'} eq $tmp_file and &stopErr("[Err] Do not use filename $tmp_file\n"); 

-d "Final" or mkdir("Final", 0755); 
if ( defined $opts{'rawLibLis'} ) {
	my $fh = &openFH($opts{'rawLibLis'}, "<"); 
	while (<$fh>) {
		chomp; 
		my @ta = split(/\s/, $_); 
		$ta[0] eq $tmp_file and &stopErr("[Err] Do not use filename $tmp_file\n"); 
		push(@inLibs, $ta[0]); 
	}
	close($fh); 
}
scalar(@inLibs) > 0 or &usage(); 


for my $fn (@inLibs) {
	my $resultFn = "${fn}.noProtFinal"; 
	unlink($resultFn); 
	unlink($tmp_file); 
	&exeCmd("perl $opts{'pl_dealFa'} -frag_head -frag_width 80 -frag 0-0 $fn | perl $opts{'pl_dealFa'} -chopKey \':\\d+\-\\d+\$' > $tmp_file"); 
	my $circN = 0; 
	&exeCmd("perl $opts{'pl_dealFa'} -attr key:len $tmp_file > prev_kl"); 
	my $no_diff = 0; 
	until ( $no_diff == 1) {
		$circN ++; 
		&tsmsg("[Rec] circle $circN for [$fn]\n"); 
		&exeCmd("blastx -evalue $opts{'evalue'} -db $opts{'dbProt'} -num_threads $opts{'cpuN'} -query $tmp_file -out ${tmp_file}.toPDB.bx0"); 
		&exeCmd("perl $opts{'pl_ProtExcluder'} ${tmp_file}.toPDB.bx0 $tmp_file"); 
		&exeCmd("perl $opts{'pl_dealFa'} -attr key:len ${tmp_file}noProtFinal > curr_kl"); 
		my $f1h = &openFH("prev_kl", '<'); 
		my $f2h = &openFH('curr_kl', '<'); 
		my (%prev_kl, %curr_kl); 
		while (<$f1h>) {
			chomp; 
			my @ta = split(/\t/, ); 
			$prev_kl{$ta[0]} = $ta[1]; 
		}
		while (<$f2h>) {
			chomp; 
			my @ta = split(/\t/, ); 
			$curr_kl{$ta[0]} = $ta[1]; 
		}
		close($f1h); 
		close($f2h); 
		open SS,'>',"saved.kl" or &stopErr("[Err] Failed to open file saved.kl. $!\n"); 
		my $saveN = 0; 
		for my $tk (keys %curr_kl) {
			defined $prev_kl{$tk} or &stopErr("[Err] Changed seq id [$tk].\n"); 
			if ( $prev_kl{$tk} eq $curr_kl{$tk} ) {
				delete $prev_kl{$tk}; 
				$saveN ++; 
				print SS "$tk\n"; 
			}
		} 
		close SS; 
		my $diffN = scalar(keys %prev_kl); 
		$diffN == 0 and $no_diff = 1; 

		if ( $saveN > 0 ) {
			&exeCmd("perl $opts{'pl_dealFa'} ${tmp_file}noProtFinal -drawByList -drawList saved.kl -drawWhole -drawLcol 0 >> $resultFn"); 
			&tsmsg("[Msg] Rest $diffN differences.\n"); 
			$diffN > 0 and &exeCmd("perl $opts{'pl_dealFa'} ${tmp_file}noProtFinal -drawByList -drawList saved.kl -drawWhole -drawLcol 0 -dropMatch > $tmp_file"); 
			$diffN > 0 and &exeCmd("perl $opts{'pl_dealFa'} ${tmp_file} -attr key:len > prev_kl"); 
		}
		&exeCmd("rm ${tmp_file}.toPDB.bx0* ${tmp_file}.ssi ${tmp_file}nPr temp"); 
	}
	&exeCmd("rm prev_kl curr_kl"); 
	# &exeCmd("mv ${tmp_file}noProtFinal Final/${fn}noProtFinal"); 
	&exeCmd("mv $resultFn Final/"); 
}



