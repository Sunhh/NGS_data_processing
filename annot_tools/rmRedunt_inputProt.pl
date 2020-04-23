#!/usr/bin/perl
# This is used to remove redundant protein sequences from given file. 
# Reference : 
#   http://www.molecularevolution.org/molevolfiles/exercises/augustus/training.html
# Maybe next : 
#   https://www.biostars.org/p/14100/
use strict; 
use warnings; 
use LogInforSunhh; 

use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"prot_qry:s", # required. 
	"prot_db:s",  # This is prot_qry if not given; 
	"opref:s", # self
	"blastp_para:s", # Do not include -outfmt
	"minIdentity:f", # Minimum identity% required to accept an alignment. 
	"help!", 
); 

sub usage {
	print <<HH; 
##########################################################################################
# perl $0 -prot_qry incomplete_proteins.fasta -opref out_prefix
#
# -help 
#
# -opref          [self]
# -blastp_para    [' -evalue 1e-5 -seg yes -num_threads 10']
# 
# -minIdentity    [80] 0-100. Min similarity% accepted for a valid alignment. 
#
# -prot_db        [-prot_qry if not given] 
##########################################################################################
HH
	exit 1; 
}

defined $opts{'prot_qry'} or &usage(); 
$opts{'blastp_para'} //= " -evalue 1e-5 -seg yes -num_threads 10"; 
$opts{'minIdentity'} //= 80; 
$opts{'opref'} //= 'self'; 
$opts{'prot_db'} //= $opts{'prot_qry'}; 
$opts{'help'} and &usage(); 

my $outfmt = "-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen frames stitle'"; 

&exeCmd_1cmd("makeblastdb -dbtype prot -in $opts{'prot_db'}"); 
&exeCmd_1cmd("blastp -query $opts{'prot_qry'} -db $opts{'prot_db'} -out $opts{'opref'}.bp6 $opts{'blastp_para'} $outfmt"); 
&exeCmd_1cmd("rm $opts{'prot_db'}.phr $opts{'prot_db'}.pin $opts{'prot_db'}.psq"); 

# Find good, redundant alignments. 
open F,'<',"$opts{'opref'}.bp6" or &stopErr("[Err] [$opts{'opref'}.bp6]: $!\n"); 
open O_BP,'>',"$opts{'opref'}.bp6.redundAln" or &stopErr("[Err] [$opts{'opref'}.bp6.redundAln]: $!\n"); 
my (%cnt, %skip, %rm); 
while (<F>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	my ($qid, $sid) = @ta[0,1]; 
	$qid eq $sid and next; 
	# This is recommended by Augustus somewhere. 
	# http://www.molecularevolution.org/molevolfiles/exercises/augustus/training.html
	$ta[2] >= $opts{'minIdentity'} or next; 
	$ta[3] * $ta[2] >= $ta[12] * $opts{'minIdentity'} or next; 
	# @ta[6,7] = sort { $a<=>$b } @ta[6,7]; 
	# @ta[8,9] = sort { $a<=>$b } @ta[8,9]; 
	print O_BP "$_\n"; 
	my $tk = "$ta[0]\t$ta[1]"; 
	$ta[11] =~ s!^\s+|\s+$!!g; 
	$cnt{$tk} = [ @ta[12,13,10,11] ]; # [ qlen, slen, evalue, bitscore ]
}
close O_BP; 
close F; 

# Output keys to be skipped. 
for my $tk ( sort { $cnt{$a}[2] <=> $cnt{$b}[2] || $cnt{$b}[3] <=> $cnt{$a}[3] || $cnt{$b}[0] <=> $cnt{$a}[0] } keys %cnt ) {
	my ($id1, $id2) = split(/\t/, $tk); 
	my ($len1, $len2, $evalue, $score) = @{$cnt{$tk}}; 
	defined $skip{$tk} and next; 
	defined $skip{$id1} and next; 
	if ( $opts{'prot_qry'} eq $opts{'prot_db'} ) {
		if ( $len2 > $len1 ) {
			if ( !defined $rm{$id2} ) {
				$rm{$id1} = $cnt{$tk}[0]; 
			}
		} else {
			if ( !defined $rm{$id1} ) {
				$rm{$id2} = $cnt{$tk}[1]; 
			}
		}
	} else {
		$skip{$id1} = 1; 
		$rm{$id1} = $cnt{$tk}[0]; 
	}
	$skip{$tk} = 1; 
	$skip{"$id2\t$id1"} = 1; 
}
open O_RM,'>',"$opts{'opref'}.bp6.redund_list" or &stopErr("[Err] [$opts{'opref'}.bp6.redundID]: $!\n"); 
for my $id ( sort keys %rm ) {
	print O_RM "$id\n"; 
}
close O_RM; 



