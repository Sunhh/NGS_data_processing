#!/usr/bin/perl -w 
# 2013-08-20 I think it is good enough. 
# 2014-03-12 Add -NMperc and -maxRdLen to filter "NM" tag. Cause "bwa aln -n float_number" doesn't work properly. 
# 2014-09-29 Add -showNumber to get include/exclude flag numbers from the selection requirement. 
use strict; 
use Getopt::Long; 

use SeqAlnSunhh; # canonical SAM-FLAG keep/drop engine (mk_flag)
my %opts; 

GetOptions(\%opts, 
	"drop:s", 
	"keep:s", 
	"NMperc:f", "maxRdLen:i", 
	"XT_tag:s", # U/R/N/M
	"anyEnd_pair!", 
	"anyEnd_info!", 
	"bothEnd_pair!", 
	"h2diff_pair!", 
	"h2diff_F!", 
	"onlyMapRd!", 
	"plusRd!", 
	"minusRd!", 
	"sameRef!", 
	"diffRef!", 
	"showNumber!", 
	"help!", 
); 

my $help_doc = <<HELP; 
###################################################################
perl $0 in_PE_aln.sam 
# I think it is enough for use. 

At least set one of -drop and -keep parameters. The alignment must support 
both paramters if they are assigned at the same time. 

-drop        Required if -keep not assigned. Example: 
               -drop 2=1,3=1;2=0,3=0,4=0,5=0;2=0,3=0,4=1,5=1
               This rule will remove read pairs without any end mapped to the reference, 
               as well as paired mapping with the same direction for both sides. 
               \";\" is calculated as logical OR  , while 
               \",\" is calculated as logical AND . 
-keep        Required if -drop not assigned. The format is same to -drop. 

-sameRef     [Boolean] Need the read pair aligned to a same reference sequence. 
-diffRef     [Boolean] Need the read pair aligned to different reference sequences. 

-XT_tag      ['U/R/M/N'], separated by \",\". I don't like this parameter. But sometimes useful. 

-NMperc      [float number]. Maximum mismatch\% allowed in the main alignment. 
-maxRdLen    [INT]. Maximum read length in query. 

------------- Predefined sets. They should not be assigned at the same time. 
------------- This will overlap the parameters of -keep and -drop. 
-anyEnd_pair   -keep 0=1 -drop 2=1,3=1
-anyEnd_info   -drop 2=1,3=1
                 # this will keep single-end, aligned reads too. 
-bothEnd_pair  -keep 0=1,2=0,3=0 
-h2diff_pair   -keep 0=1,2=0,3=0,4=0,5=1;0=1,2=0,3=0,4=1,5=0
                 # Head to head mapping discarding the position order. 
-h2diff_F      -keep 0=1,2=0,3=0,4=0,5=1
                 # Used to get alignments with insert size. 
-onlyMapRd     -drop 2=1
-plusRd        -keep 2=0,4=0
-minusRd       -keep 2=0,4=1

-showNumber    Only show the flag numbers with Include/Exclude information. 

###################################################################
# BitPos	Flag	Chr	Description
# 0	0x0001	p	the read is paired in sequencing
# 1	0x0002	P	the read is mapped in a proper pair
# 2	0x0004	u	the query sequence itself is unmapped
# 3	0x0008	U	the mate is unmapped
# 4	0x0010	r	strand of the query (1 for reverse)
# 5	0x0020	R	strand of the mate
# 6	0x0040	1	the read is the first read in a pair
# 7	0x0080	2	the read is the second read in a pair
# 8	0x0100	s	the alignment is not primary
# 9	0x0200	f	the read fails platform/vendor quality checks
# 10	0x0400	d	the read is either a PCR or an optical duplicate
###### FLAG list prepared. 
HELP

if        ($opts{anyEnd_pair}) {
	$opts{keep} = '0=1'; 
	$opts{drop} = '2=1,3=1'; 
} elsif   ($opts{anyEnd_info}) {
	$opts{keep} = ''; 
	$opts{drop} = '2=1,3=1'; 
} elsif   ($opts{bothEnd_pair}) {
	$opts{keep} = '0=1,2=0,3=0'; 
	$opts{drop} = ''; 
} elsif   ($opts{h2diff_pair}) {
	$opts{keep} = '0=1,2=0,3=0,4=0,5=1;0=1,2=0,3=0,4=1,5=0'; 
	$opts{drop} = ''; 
} elsif   ($opts{h2diff_F}) { 
	$opts{keep} = '0=1,2=0,3=0,4=0,5=1'; 
	$opts{drop} = ''; 
} elsif   ($opts{onlyMapRd}) {
	$opts{keep} = ''; 
	$opts{drop} = '2=1'; 
} elsif   ($opts{plusRd}) {
	$opts{keep} = '2=0,4=0'; 
	$opts{drop} = ''; 
} elsif   ($opts{minusRd}) {
	$opts{keep} = '2=0,4=1'; 
	$opts{drop} = ''; 
#} elsif   ($opts{}) {
#	$opts{keep} = ''; 
#	$opts{drop} = ''; 
}

sub usage {
	print "$help_doc\n"; 
	exit 1; 
}
if (-t and !@ARGV) {
	&usage(); 
} elsif ( $opts{help} ) {
	&usage(); 
}


###### Make TAG list 
my $need_xt_tag = 0; 
my %good_xt_tag; 
if (defined $opts{XT_tag}) {
	$need_xt_tag = 1; 
	for (split(/,/, $opts{XT_tag})) {
		$_ = uc($_); 
		$good_xt_tag{"XT:A:$_"} = 1; 
	}
}
###### End make TAG list 

###### Set maximum NM% allowed. 
my $need_nm_tag = 0; 
my %max_nmN; 
my $maxRdLen = (defined $opts{maxRdLen}) ? $opts{maxRdLen} : 1e3 ; 
# $maxRdLen < 1 and $maxRdLen = 1; 
$maxRdLen < 1 and die "-maxRdLen=$maxRdLen?\n"; 
if (defined $opts{NMperc}) {
	$opts{NMperc} >= 0 or die "-NMperc=$opts{NMperc} < 0 ?\n"; 
	$need_nm_tag = 1; 
	for my $rdLen ( 1..$maxRdLen ) {
		$max_nmN{$rdLen} = int( $rdLen * $opts{NMperc} ); 
	}
}
###### End set maximum NM% allowed. 


###### Build the include/exclude flag set via MyPM/SeqAlnSunhh::mk_flag (canonical engine).
my ($have_drop, $have_keep) = (0, 0);
defined $opts{drop} and $opts{drop} ne '' and $have_drop = 1;
defined $opts{keep} and $opts{keep} ne '' and $have_keep = 1;
($have_drop or $have_keep) or die "\n\nNo filter parameter accepted! Please check your input for -drop/-keep!\n\n$help_doc\n";
my $good_flag = &SeqAlnSunhh::mk_flag( 'keep'=>$opts{keep}, 'drop'=>$opts{drop} );
# Preserve the original %is_output representation: 1 for kept, 0 for excluded, defined for all 0..4095.
my %is_output;
for my $f ( 0 .. 4095 ) { $is_output{$f} = $good_flag->{$f} ? 1 : 0; }
###### Set FLAGs OK.

###### Show flag number usage: If we only need to show flag number usage. 
if ( $opts{showNumber} ) {
	for ( 1 .. 4095 ) {
		if ( defined $is_output{$_} ) {
			print STDOUT join("\t", $_, 'Include')."\n"; 
		} else {
			print STDOUT join("\t", $_, 'Exclude')."\n"; 
		}
	}
	exit(0); 
}
###### End showing flag number usage. 

# exit; 

while (<>) {
	if (m/^@/) {
		print; 
	}else{
		chomp; 
		my $is_op = 1; 
		m/^[^\t]+\t(\d+)/ or die "Failed to parse line:\n$_\n"; 
		my $tmp_flag = $1; 

		$is_output{$tmp_flag} or $is_op = 0; 

		if ($opts{sameRef}) {
			/^(?:[^\t]+\t){6}\=\t/ or $is_op = 0; 
		} elsif ($opts{diffRef}) {
			/^(?:[^\t]+\t){6}\=\t/ and $is_op = 0; 
		}

		# Checking 
		my @ta; 
		if ($is_op == 1 and $need_xt_tag == 1) {
			@ta = split(/\t/, $_); 
		} elsif ($is_op == 1 and $need_nm_tag == 1) {
			@ta = split(/\t/, $_); 
		}

		# Check XT: tag
		if ($is_op == 1 and $need_xt_tag == 1) {
			if ($need_xt_tag == 1) {
				my $is_xt_ok = 0; 
				CHK_XT: 
				for my $tb (@ta[11 .. $#ta]) {
					defined $good_xt_tag{$tb} and do { $is_xt_ok = 1; last CHK_XT; }; 
				}
				$is_xt_ok == 1 or $is_op = 0; 
			}#End if CHK_XT; 
		}

		# Check NM:i: tag 
		if ($is_op == 1 and $need_nm_tag == 1) {
			my $rdLen = length($ta[9]); 
			my $cur_nmN = 0; 
			CHK_NM: 
			for my $tb (@ta[11 .. $#ta]) {
				if ( $tb =~ m/^NM:i:(\d+)$/ ) {
					$1 <= $max_nmN{$rdLen} or $is_op = 0; 
					last CHK_NM; 
				}
			}
		}

		$is_op == 1 and print STDOUT "$_\n"; 
		last; 
	}
}
while (<>) {
	chomp; 
	my $is_op = 1; 
	m/^[^\t]+\t(\d+)/ or die "Failed to parse line:\n$_\n"; 
	my $tmp_flag = $1; 

	$is_output{$tmp_flag} or $is_op = 0; 

	if ($opts{sameRef}) {
		/^(?:[^\t]+\t){6}\=\t/ or $is_op = 0; 
	} elsif ($opts{diffRef}) {
		/^(?:[^\t]+\t){6}\=\t/ and $is_op = 0; 
	}

	# Checking 
	my @ta; 
	if ($is_op == 1 and $need_xt_tag == 1) {
		@ta = split(/\t/, $_); 
	} elsif ($is_op == 1 and $need_nm_tag == 1) {
		@ta = split(/\t/, $_); 
	}

	# Check XT: tag
	if ($is_op == 1 and $need_xt_tag == 1) {
		if ($need_xt_tag == 1) {
			my $is_xt_ok = 0; 
			CHK_XT: 
			for my $tb (@ta[11 .. $#ta]) {
				defined $good_xt_tag{$tb} and do { $is_xt_ok = 1; last CHK_XT; }; 
			}
			$is_xt_ok == 1 or $is_op = 0; 
		}#End if CHK_XT; 
	}

	# Check NM:i: tag 
	if ($is_op == 1 and $need_nm_tag == 1) {
		my $rdLen = length($ta[9]); 
		my $cur_nmN = 0; 
		CHK_NM: 
		for my $tb (@ta[11 .. $#ta]) {
			if ( $tb =~ m/^NM:i:(\d+)$/ ) {
				$1 <= $max_nmN{$rdLen} or $is_op = 0; 
				last CHK_NM; 
			}
		}
	}
	
	$is_op == 1 and print STDOUT "$_\n"; 

}

