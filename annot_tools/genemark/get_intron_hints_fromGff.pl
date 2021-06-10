#!/usr/bin/perl -w
### Designed referring filterIntronsFindStrand.pl in BRAKER; 
use strict; 
use warnings; 
use fileSunhh; 
use LogInforSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts,
	"help!",
	"inGff:s",       
	"coverStrand!", 
	"check_ya!",     # Boolean, check intron-exon boundary or not. 
	"good_ya:s@", # GTAG,GCAG,ATAC
	"genome_fa:s",   # Required for -coverStrand and -good_ya; 
	"pl_dealGff:s",  # /home/Sunhh/tools/github/NGS_data_processing/temp/deal_gff3.pl
); 

$opts{'pl_dealGff'} = '/home/Sunhh/tools/github/NGS_data_processing/temp/deal_gff3.pl'; 

my $help_txt = <<HH; 
######################################################################
# perl $0 -inGff WM97_v1.scf.annot_prot.gff3 
#
# -pl_dealGff     [$opts{'pl_dealGff'}]
# -check_ya       [Boolean] Check if the intron-exon boundary is good; 
#   -good_ya      [String] Good intron-exon boundary combinations accepted; GTAG,GCAG,ATAC; 
#   -coverStrand  [Boolean] Overwrite strand information according to the -good_ya; 
#   -genome_fa    [String] genome.fasta file containing the scaffold sequences;
#
HH

defined $opts{'inGff'} or &LogInforSunhh::usage($help_txt); 
defined $opts{help} and &LogInforSunhh::usage($help_txt); 



my %genome; 
my @good_comb; 
defined $opts{'coverStrand'} and $opts{'check_ya'} = 1; 
if (defined $opts{'check_ya'}) {
	$opts{'good_ya'} //= ['GTAG,GCAG,ATAC']; 
	my %t1h; 
	for my $t1 (split(/[\s,]/, join(',', @{$opts{'good_ya'}}))) {
		defined $t1h{uc($t1)} and next; 
		$t1h{uc($t1)} = 1; 
		push(@good_comb, uc($t1)); # Not using %hash because the order matters. 
	}
	defined $opts{'genome_fa'} or die "-genome_fa required\n"; 
	use fastaSunhh; 
	my $fs_obj = fastaSunhh->new(); 
	%genome = %{$fs_obj->save_seq_to_hash( 'faFile'=>$opts{'genome_fa'})}; 
	for (keys %genome) {
		$genome{$_}{'seq'} =~ s!\s!!g; 
		$genome{$_}{'seq'} = uc($genome{$_}{'seq'}); 
	}
}

my $tmpDir = &fileSunhh::new_tmp_dir('create'=>1); 
&exeCmd_1cmd("perl $opts{'pl_dealGff'} -inGff $opts{'inGff'} -list_intron -intron_byFeat CDS > $tmpDir/intron.tbl"); 
open F,'<',"$tmpDir/intron.tbl" or die "$tmpDir/intron.tbl, $!\n"; 
my @lines; 
my %multV; 
while (<F>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	$ta[0] eq 'ParentID' and $. == 1 and next; 
	push(@lines, [@ta[2,3,4,5]]); 
	my $k = "$ta[2]\t$ta[3]\t$ta[4]\t$ta[5]"; 
	$multV{$k} ++; 
	$multV{$k} > 1 and next; 
	# print join("\t", $ta[2], "fromGff", 'intron', $ta[3], $ta[4], 1, $ta[5], "mult=1;pri=4;src=E")."\n"; 
}
close F; 

unless (defined $opts{'check_ya'}) {
	for my $t1 (sort { $a->[0] cmp $b->[0] || $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2] } @lines) {
		my $k = join("\t", @{$t1}[0,1,2,3]); 
		print STDOUT join("\t", $t1->[0], qw/fromGff intron/, $t1->[1], $t1->[2], $multV{$k}, $t1->[3], '.', "mult=$multV{$k};pri=4;src=E")."\n"; 
	}
	exit(0); 
}

if (defined $opts{'check_ya'}) {
	my %multV2; 
	my %hasLine2; 
	my @lines2; 
	for my $t1 (sort { $a->[0] cmp $b->[0] || $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2] } @lines) {
		my $k = join("\t", @{$t1}[0,1,2,3]); # scfID, intron_S, intron_E, intron_Strand
		defined $genome{$t1->[0]} or do { &tsmsg("[Wrn] Failed to find seq for ScfID [$t1->[0]], which is skipped\n"); next; }; 
		my $chk = &chk_boundary($t1->[0], $t1->[1], $t1->[2]); 
		$chk eq '' and next; 
		if (defined $opts{'coverStrand'}) {
			$t1->[3] = $chk; 
		} else {
			$t1->[3] =~ m!^(\+|\-)$! or $t1->[3] = $chk; 
		}
		my $k2 = join("\t", @{$t1}[0,1,2,3]); 
		$multV2{$k2} += $multV{$k}; 
		defined $hasLine2{$k2} and next; 
		$hasLine2{$k2} = 1; 
		push(@lines2, [@$t1]); 
	}
	for my $t1 (@lines2) {
		my $k2 = join("\t", @{$t1}[0,1,2,3]); 
		print STDOUT join("\t", $t1->[0], qw/fromGff intron/, $t1->[1], $t1->[2], $multV2{$k2}, $t1->[3], '.', "mult=$multV2{$k2};pri=4;src=E")."\n"; 
	}
}

sub chk_boundary {
	my ($scfID, $intronS, $intronE) = @_; 
	my $cc = substr($genome{$scfID}{'seq'}, $intronS-1, 2) . substr($genome{$scfID}{'seq'}, $intronE-2, 2); # All upper-case; 
	my $rr = reverse($cc); $rr =~ tr/acgturykmbvdhACGTURYKMBVDHwWsSnN/tgcaayrmkvbhdTGCAAYRMKVBHDwWsSnN/; 
	for my $c1 (@good_comb) {
		if ( $c1 eq $cc ) {
			return('+'); 
		} elsif ( $c1 eq $rr ) {
			return('-'); 
		}
	}
	return(''); 
}# chk_boundary() 

&fileSunhh::_rmtree($tmpDir); 

