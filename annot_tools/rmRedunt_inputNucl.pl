#!/usr/bin/perl
# 2015-09-26 
# This is used to remove redundant Nucleotide sequences from given file. 
# This is edited from rmRedunt_inputProt.pl 
# I probably fail when there are N gaps, but I don't want to use index method. 
# Reference : 
#   http://www.molecularevolution.org/molevolfiles/exercises/augustus/training.html
# Maybe next : 
#   https://www.biostars.org/p/14100/
use strict; 
use warnings; 
use LogInforSunhh; 
use fastaSunhh; 
my $fasS_obj = fastaSunhh->new(); 

use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"nucl_qry:s", # required. 
	"opref:s", # self
	"blastn_para:s", # Do not include -outfmt
	"minIdentity:f", # Minimum identity% required to accept an alignment. 
	"help!", 
); 

sub usage {
	print <<HH; 
##########################################################################################
# perl $0 -nucl_qry redundant_nucleotides.fasta -opref out_prefix
#
# -help 
#
# -opref          [self]
# -blastn_para    [' -task blastn -evalue 1e-5 -dust yes -max_target_seqs 5 -num_threads 10']
# 
# -minIdentity    [100] 0-100. Min similarity% accepted for a valid alignment. 
##########################################################################################
HH
	exit 1; 
}

defined $opts{'nucl_qry'} or &usage(); 
$opts{'blastn_para'} //= " -task blastn -evalue 1e-5 -dust yes -max_target_seqs 5 -num_threads 10"; 
$opts{'minIdentity'} //= 100; 
$opts{'opref'} //= 'self'; 
$opts{'help'} and &usage(); 

my $outfmt = "-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen frames stitle'"; 
$outfmt = "-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand'"; 

&exeCmd_1cmd("makeblastdb -dbtype nucl -in $opts{'nucl_qry'}"); 
&exeCmd_1cmd("blastn -query $opts{'nucl_qry'} -db $opts{'nucl_qry'} -out $opts{'opref'}.bn6 $opts{'blastn_para'} $outfmt"); 
&exeCmd_1cmd("rm $opts{'nucl_qry'}.nsq $opts{'nucl_qry'}.nin $opts{'nucl_qry'}.nhr"); # Sometimes the input is so big that it produces input.??.n?? files, but I'd like manually remove them in that case. 

# Find good, redundant alignments. 
open F,'<',"$opts{'opref'}.bn6" or &stopErr("[Err] [$opts{'opref'}.bn6]: $!\n"); 
open O_BP,'>',"$opts{'opref'}.bn6.redundAln" or &stopErr("[Err] [$opts{'opref'}.bn6.redundAln]: $!\n"); 
# Edit here. 
my (%cnt, %skip, %rm); 
while (<F>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	my ($qid, $sid) = @ta[0,1]; 
	$qid eq $sid and next; 
	$ta[2] >= $opts{'minIdentity'} or next; 
	$ta[3] * $ta[2] >= $ta[12] * $opts{'minIdentity'} or next; 
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
	if ( $len2 > $len1 ) {
		if ( !defined $rm{$id2} ) {
			$rm{$id1} = $cnt{$tk}[0]; 
		}
	} else {
		if ( !defined $rm{$id1} ) {
			$rm{$id2} = $cnt{$tk}[1]; 
		}
	}
	$skip{$tk} = 1; 
	$skip{"$id2\t$id1"} = 1; 
}
open O_RM,'>',"$opts{'opref'}.bn6.redund_list" or &stopErr("[Err] [$opts{'opref'}.bn6.redundID]: $!\n"); 
for my $id ( sort keys %rm ) {
	print O_RM "$id\n"; 
}
close O_RM; 

my %seq = %{ $fasS_obj->save_seq_to_hash('faFile' => $opts{'nucl_qry'}) }; 
open O_FA,'>',"$opts{'opref'}.NR.fa" or die "$!\n"; 
for my $id (sort { $seq{$a}{'Order'} <=> $seq{$b}{'Order'} } keys %seq) {
	defined $rm{$id} and next; 
	chomp($seq{$id}{'seq'}); 
	print O_FA ">$seq{$id}{'head'}\n$seq{$id}{'seq'}\n"; 
}
close O_FA; 

