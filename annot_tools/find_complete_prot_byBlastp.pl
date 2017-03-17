#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 

use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"prot_db:s", # Default Nr 
	"prot_qry:s", # required. 
	"opref:s", # out
	"blastp_para:s", # Do not include -outfmt
	"dist2left:i", # Max distance allowed to left. 
	"dist2right:i", # Max distance allowed to right. 
	"minIdentity:f", # Minimum identity% required to accept an alignment. 
	"maxSizeDiffR:f", # Maximum size difference between qry and sbj protein length based on the shorter one. 
	"topNhits:i", 
	"help!", 
); 

sub usage {
	print <<HH; 
##########################################################################################
# perl $0 -prot_qry incomplete_proteins.fasta -opref out_prefix
#
# -help 
#
# -opref          [out]
# -prot_db        [nr]
# -topNhits       [5] I don't want to use too high because the protein database may 
#                     not be very reliable. 
# -blastp_para    [' -evalue 1e-5 -seg yes -num_threads 10']
# -dist2left      [2] Max distance allowed to the left end of query / subject. 
# -dist2right     [2] Max distance allowed to the right end of query / subject. 
#
# -maxSizeDiffR   [-1] Maximum size difference between qry and sbj protein length based on the shorter one. 
# 
# -minIdentity    [0] 0-100. Min similarity% accepted for a valid alignment. 
##########################################################################################
HH
	exit 1; 
}

$opts{'prot_db'} //= 'nr'; 
defined $opts{'prot_qry'} or &usage(); 
$opts{'blastp_para'} //= " -evalue 1e-5 -seg yes -num_threads 10"; 
$opts{'dist2left'} //= 2; 
$opts{'dist2right'} //= 2; 
$opts{'maxSizeDiffR'} //= -1; 
$opts{'minIdentity'} //= 0; 
$opts{'topNhits'} //= 5; 
$opts{'opref'} //= 'out'; 
$opts{'help'} and &usage(); 

my $outfmt = "-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen frames stitle'"; 


&exeCmd_1cmd("blastp -query $opts{'prot_qry'} -out $opts{'opref'}.bp6 -db $opts{'prot_db'} $opts{'blastp_para'} $outfmt"); 

open F,'<',"$opts{'opref'}.bp6" or &stopErr("[Err] [$opts{'opref'}.bp6]: $!\n"); 
open O_BP,'>',"$opts{'opref'}.bp6.complete" or &stopErr("[Err] [$opts{'opref'}.bp6.complete]: $!\n"); 
my %complete_keys; 
my %cnt; 
while (<F>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	my ($qid, $sid) = @ta[0,1]; 
	defined $cnt{$qid}{'hit'}{$sid} and next; 
	$cnt{$qid}{'cnt'} ++; 
	$cnt{$qid}{'cnt'} > $opts{'topNhits'} and next; 
	$cnt{$qid}{'hit'}{$sid} = 1; 
	$ta[2] >= $opts{'minIdentity'} or next; 
	unless ( $opts{'maxSizeDiffR'} < 0 ) {
		my $dlen = $ta[12]; 
		$ta[13] < $dlen and $dlen = $ta[13]; 
		my $ddiff = abs($ta[12]-$ta[13]); 
		$ddiff <= $dlen * $opts{'maxSizeDiffR'} or next; 
	}
	@ta[6,7] = sort { $a<=>$b } @ta[6,7]; 
	@ta[8,9] = sort { $a<=>$b } @ta[8,9]; 
	$ta[6]-1 <= $opts{'dist2left'} or next; 
	$ta[12]-$ta[7] <= $opts{'dist2right'} or next; 
	$ta[8]-1 <= $opts{'dist2left'} or next; 
	$ta[13]-$ta[9] <= $opts{'dist2right'} or next; 
	print O_BP "$_\n"; 
	$complete_keys{$qid} = 1; 
}
close O_BP; 
close F; 


