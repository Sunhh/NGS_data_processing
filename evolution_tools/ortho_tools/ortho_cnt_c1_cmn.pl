#!/usr/bin/perl
# Get ortho-groups satisfying the rules C1 : At least \$c1_required_taxaNum in \%c1_required_taxaLis; 
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 

use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"in_file:s", # all_orthomcl.out , OrthologousGroups.csv , Orthogroups.csv
	"in_fmt:s",  # Default orthomcl , could be 'orthofinder'
	"required_taxaLis:s@", 
	"required_taxaNum:i@", # Minimum numebr required. 
	"required_taxaNum_max:i@", # Maximum number required. 
	"avoid_taxaLis:s@", 
	"avoid_taxaNum:i@",    # Maximum number allowed. 
	"class_ID:s", # c1 
	"rest_file:s", 
); 

$opts{'class_ID'} //= 'c1'; 
$opts{'in_fmt'}   //= 'orthofinder'; 

my $help_txt = <<HH; 
perl $0 -in_file Orthogroups.csv   -required_taxaLis taxa_list_1col -required_taxaNum taxa_min_num 

-in_fmt     [$opts{'in_fmt'}]

-class_ID   [$opts{'class_ID'}]
-rest_file  [in_file.rest]

-required_taxaLis     [filename] Format : taxaID_1 \\n taxaID_2 \\n taxaID_3 \\n ... 
-required_taxaNum     [number] Minimum number of taxa required in OG. 
-required_taxaNum_max [number] Maximum number of taxa allowed in OG. -1 (default) means no limit. 

-avoid_taxaLis        [filename] 
-avoid_taxaNum        [number] Use 0 for exclude any taxa. 

HH

$opts{'help'} and &LogInforSunhh::usage($help_txt); 

defined $opts{'in_file'} or &LogInforSunhh::usage($help_txt); 
$opts{'rest_file'} //= "$opts{'in_file'}.rest"; 
my $ofh_rest = &openFH($opts{'rest_file'}, '>'); 
my $has_header_rest = 0; 

my %glob; 
&prepare_glob(); 
my ( @ortho_grps , %taxGen2grpID);
{
	my ($r1, $r2) = &load_grps( $opts{'in_file'} , $glob{'in_fmt'} );
	@ortho_grps = @$r1;
	%taxGen2grpID = %$r2;
}


$glob{'required_taxaLis'}     //= [{}]; 
$glob{'required_taxaNum'}     //= [0]; 
$glob{'required_taxaNum_max'} //= [-1]; 

if (defined $opts{'required_taxaLis'}) {
	for (my $i=0; $i<@{$opts{'required_taxaLis'}}; $i++) {
		( $glob{'required_taxaLis'}[$i] ) =   &load_taxaLis( $opts{'required_taxaLis'}[$i] );
		$glob{'required_taxaAll'}[$i]     =   scalar(keys %{$glob{'required_taxaLis'}[$i]});
		$opts{'required_taxaNum'}[$i]     //= $glob{'required_taxaAll'}[$i];
		$glob{'required_taxaNum'}[$i]     =   $opts{'required_taxaNum'}[$i];
		$opts{'required_taxaNum_max'}[$i] //= -1; 
		$glob{'required_taxaNum_max'}[$i] =   $opts{'required_taxaNum_max'}[$i]; 
	}
}
if (defined $opts{'avoid_taxaLis'}) {
	for (my $i=0; $i<@{$opts{'avoid_taxaLis'}}; $i++) {
		( $glob{'avoid_taxaLis'}[$i] ) =   &load_taxaLis( $opts{'avoid_taxaLis'}[$i] ); 
		$glob{'avoid_taxaAll'}[$i]     =   scalar(keys %{$glob{'avoid_taxaLis'}[$i]}); 
		$opts{'avoid_taxaNum'}[$i]     //= 0; 
		$glob{'avoid_taxaNum'}[$i]     =   $opts{'avoid_taxaNum'}[$i]; 
	}
}

my %cnt; 
print STDOUT join("\t", qw/Class_ID OG_ID taxID geneID/)."\n"; 
for my $r1 (@ortho_grps) {
	if ($r1->{'grpID'} eq '') {
		# print STDOUT    $r1->{'line'}."\n"; 
		$glob{'rest_header'} = "$r1->{'line'}"; 
		next; 
	}
	my @taxa_ids = sort keys %{$r1->{'tax2gene'}};
	my $is_good = 1; 
	if ( defined $opts{'required_taxaLis'} ) {
		for ( my $i=0; $i< @{$glob{'required_taxaLis'}}; $i++ ) {
			$cnt{'curr_requiredN'} = 0;
			for my $t1 (@taxa_ids) {
				defined $glob{'required_taxaLis'}[$i]{$t1} and $cnt{'curr_requiredN'}++;
			}
			$cnt{'curr_requiredN'} >= $glob{'required_taxaNum'}[$i] or do { $is_good = 0; last; };
			if ( $glob{'required_taxaNum_max'}[$i] >= 0 ) {
				$cnt{'curr_requiredN'} <= $glob{'required_taxaNum_max'}[$i] or do { $is_good = 0; last; }; 
			}
		}
	}
	if ( defined $opts{'avoid_taxaLis'} ) {
		for ( my $i=0; $i< @{$glob{'avoid_taxaLis'}}; $i++ ) {
			$cnt{'curr_avoidN'} = 0; 
			for my $t1 (@taxa_ids) {
				defined $glob{'avoid_taxaLis'}[$i]{$t1} and $cnt{'curr_avoidN'}++; 
			}
			$cnt{'curr_avoidN'}    <= $glob{'avoid_taxaNum'}[$i]    or do { $is_good = 0; last; }; 
		}
	}
	if ( $is_good == 1 ) {
		for my $taxID (@taxa_ids) {
			for my $gene_ID (sort keys %{$r1->{'tax2gene'}{$taxID}}) {
				print STDOUT join("\t", $opts{'class_ID'}, $r1->{'grpID'}, $taxID, $gene_ID)."\n"; 
			}
		}
	} else {
		if ($has_header_rest == 0) {
			$has_header_rest = 1; 
			print $ofh_rest "$glob{'rest_header'}\n"; 
		}
		print $ofh_rest $r1->{'line'}."\n"; 
	}
}



################################################################################
# Inner sub-routines. 
################################################################################

=head1 load_grps( $opts{'in_file'} )

Return : ( [ \%grp_1, \%grp_2, ... ], %gene2grpID )
  keys %grp = qw/grpID geneN taxaN tax2gene/
    $grp{'grpID'} = ORTHOMCL\d+
    $grp{'geneN'} = \d+
    $grp{'taxaN'} = \d+
    $grp{'tax2gene'} = { {taxa_1}{gene_1} => $order, {taxa_1}{gene_2} => $order, {taxa_2}{gene_3} => $order, ... }
  keys %gene2grpID = ( gene_1, gene_2, gene_3, ... )
    $gene2grpID{$taxa_1}{$gene_1} = $grpID

=cut
sub load_grps {
	my ($fn, $fmt) = @_;
	$fmt //= 'orthocml';
	if ($fmt =~ m!^\s*orthomcl\s*$!i) {
		return( &load_grps_fmt_orthocml($fn) );
	} elsif ($fmt =~ m!^\s*orthofinder\s*$!i) {
		return( &load_grps_fmt_orthofinder($fn) );
	} else {
		&stopErr("[Err] Unknown format of input Orthologous Group file.\n");
	}
	return;
}# load_grps()

sub load_grps_fmt_orthofinder {
	my $fn = shift;
	my $fh = &openFH($fn, '<');
	my (@back, %gene2grpID);
	# For OrthologousGroups.csv
	my @hh;
	while (&wantLineC($fh)) {
		chomp;
		my @ta = &splitL("\t", $_);
		if ($. == 1) {
			$ta[0] eq '' or $ta[0] eq 'ID' or &stopErr("[Err] Bad 1st line: $_\n");
			@hh = @ta;
			push(@back, {}); 
			$back[-1]{'grpID'} = ''; 
			$back[-1]{'line'}  = $_; 
			$back[-1]{'geneN'} = 0; 
			$back[-1]{'taxaN'} = $#hh; 
			next;
		}
		my ($grpID, $geneN, $taxaN) = ($ta[0], 0, 0);
		push(@back, {});
		$back[-1]{'grpID'} = $grpID;
		$back[-1]{'line'}  = $_; 
		for (my $i=1; $i<@ta; $i++) {
			$ta[$i] eq '' and next;
			$ta[$i] =~ m!^(\s*,*\s*)*$! and next;
			my $tax_name = $hh[$i];
			$taxaN ++;
			for my $tax_gene (split(/, /, $ta[$i])) {
				$geneN ++;
				$back[-1]{'tax2gene'}{$tax_name}{$tax_gene} = $geneN;
				defined $gene2grpID{$tax_name}{$tax_gene} and &stopErr("[Err] repeat ID {$tax_name}{$tax_gene}\n");
				$gene2grpID{$tax_name}{$tax_gene} = $back[-1]{'grpID'};
			}
		}
		$back[-1]{'geneN'} = $geneN;
		$back[-1]{'taxaN'} = $taxaN;
	}
	close($fh);
	return(\@back, \%gene2grpID);
}# load_grps_fmt_orthofinder()

sub load_grps_fmt_orthocml {
	my $fn = shift;
	my $fh = &openFH($fn, '<');
	my @back;
	my %gene2grpID;
	# For all_orthomcl.out
	while (<$fh>) {
		chomp;
		my @ta=split(/\t/, $_);
		$ta[0] =~ m!^(ORTHOMCL\d+)\s*\(\s*(\d+)\s*genes\,\s*(\d+)\s*taxa\)\s*:\s*$! or die "$ta[0]\n";
		my $grpID = $1;
		my $geneN = $2;
		my $taxaN = $3;
		push(@back, {});
		$back[-1]{'grpID'} = $grpID;
		$back[-1]{'line'}  = $_; 
		$back[-1]{'geneN'} = $geneN;
		$back[-1]{'taxaN'} = $taxaN;
		$ta[1] =~ s!^\s+|\s+$!!g;
		my $i=0;
		for my $tb (split(/\s+/, $ta[1])) {
			$i ++;
			$tb =~ m!^(\S+)\(([^\s()]+)\)$! or die "|$tb|\n";
			my $tax_gene = $1;
			my $tax_name = $2;
			$back[-1]{'tax2gene'}{$tax_name}{$tax_gene} = $i;
			defined $gene2grpID{$tax_name}{$tax_gene} and &stopErr("[Err] repeat ID {$tax_name}{$tax_gene}\n");
			$gene2grpID{$tax_name}{$tax_gene} = $back[-1]{'grpID'};
		}
	}
	close($fh);
	return(\@back, \%gene2grpID);
}# load_grps_fmt_orthocml ()

sub prepare_glob {
        $glob{'in_fmt'} = 'orthofinder';
        defined $opts{'in_fmt'} and $glob{'in_fmt'} = $opts{'in_fmt'};
}# prepare_glob

=head1 load_taxaLis ($filename)

Return   : ( \%taxaLis )

  $taxaLis{ $taxaID } = $lineN

=cut
sub load_taxaLis {
	my $fn = shift;
	my %back;
	my $fh = &openFH($fn, '<');
	while (&wantLineC($fh)) {
		my @ta = &splitL("\t", $_);
		$back{$ta[0]} //= $.;
	}
	close($fh);
	return(\%back);
}# load_taxaLis ()

