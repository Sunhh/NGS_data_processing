#!/usr/bin/perl
use strict; 
use warnings; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"addTag:s", # 3:, added to the head of each gene. 
	"notPasa!", # If given, do not add gene line for pasa gff. 
); 

$opts{'addTag'} //= "3:"; 

my $help_txt = <<HH; 
perl $0 -addTag '$opts{'addTag'}' P1denovoAndGG_pasa.pasa_assemblies.named.gff3 > P1denovoAndGG_pasa.pasa_assemblies.named.fmt.gff3

-notPasa        [Bool] If input is not pasa.gff3

HH

!@ARGV and die "$help_txt"; 

my $inF=shift; 
my %se; 
open F,'<',"$inF" or die; 
while (<F>) {
	chomp; 
	if (m!^\s*(#|$)!) {
		print STDOUT "$_\n"; 
		next; 
	}
	my @ta = split(/\t/, $_); 
	# if ( $ta[2] =~ m!^(gene|mrna|exon|cds|five_prime_UTR|three_prime_UTR|match|match_part)$!i ) {
	# }
	$ta[0] =~ s!\s.*$!!; 
	my ($id, $tgt); 
	if ( $opts{'notPasa'} ) {
		if ( $ta[8] =~ s!(^|\s|;)ID=([^\s;]+)!$1ID=$opts{'addTag'}$2! ) {
			$id = "$opts{'addTag'}$2"; 
		} 
		if ( $ta[8] =~ s!(^|\s|;)Parent=([^\s;]+)!$1Parent=$opts{'addTag'}$2! ) {
			$id = "$opts{'addTag'}$2"; 
		}
		defined $id or die "Need ID/Parent: $_\n"; 
		$tgt = 'NA'; 
	} else {
		$ta[8] =~ s!^ID=([^\s;]+);Target=(\S+)!Parent=$opts{'addTag'}$1;Target=$2! or die "$_\n"; 
		$id="$opts{'addTag'}$1"; 
		$tgt = $2; 
		$ta[2] = "match_part"; 
	}
	$se{$id}{tgt} //= $tgt; 
	$se{$id}{tgt} eq $tgt or die "$_\n"; 
	$se{$id}{s} //= $ta[3]; 
	$se{$id}{e} //= $ta[4]; 
	$se{$id}{str} //= $ta[6]; 
	$se{$id}{ln} //= $.; 
	$se{$id}{s} > $ta[3] and $se{$id}{s} = $ta[3]; 
	$se{$id}{e} < $ta[4] and $se{$id}{e} = $ta[4]; 
	# $ta[1] = "P1denovoAndGG_pasa"; 
	$ta[8] =~ s!;$!!; $ta[8] = "$ta[8];"; 
	push(@{$se{$id}{lines}}, [@ta]); 
}
close F; 
for my $id (sort {$se{$a}{ln} <=> $se{$b}{ln}} keys %se) {
	my %th = %{$se{$id}}; 
	unless ( $opts{'notPasa'} ) {
		print STDOUT join("\t", $th{lines}[0][0], $th{lines}[0][1], "match", $th{s}, $th{e}, '.', $th{str}, '.', "ID=$id;Name=$th{tgt};")."\n"; 
		if ( $th{str} eq '-' ) {
			@{$th{lines}} = sort { $b->[3] <=> $a->[3] } @{$th{lines}}; 
		} elsif ( $th{str} eq '+' ) {
			@{$th{lines}} = sort { $a->[3] <=> $b->[3] } @{$th{lines}}; 
		} else {
			die "str=$th{str}\n"; 
		}
	}
	for my $tr ( @{$th{lines}} ) {
		print STDOUT join("\t", @$tr)."\n"; 
	}
}
