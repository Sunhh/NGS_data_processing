#!/usr/bin/perl -w 
# 2013-08-20 I think it is good enough. 
#
use strict; 
use Getopt::Long; 

my %opts; 

GetOptions(\%opts, 
	"drop:s", 
	"keep:s", 
	"NM_tag:s", # Not used now. 
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

-XT_tag      'U/R/M/N'. I don't like this parameter. But sometimes useful. 
-sameRef
-diffRef

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



###### Make FLAG list 
my %flag; 
for (0 .. 2047) {
	my $binum = unpack("B32", pack("N", $_)); 
	$flag{$_} = [ reverse( split(//, sprintf("%011d", $binum) ) ) ]; 
}


######### Set FLAGs in or out. 
# Writing format for -drop/-keep parameter : 
# -drop 2=1;4=0,5=0
# -keep 4=1,5=1;3=1
# drop and keep are both considered at the same time if assigned. 
#   rules separated by ";" are calculated as logical "OR" , 
#   rules separated by "," are calculated as logical "AND" . 
my ($have_drop, $have_keep) = (0, 0); 
defined $opts{drop} and $opts{drop} ne '' and $have_drop = 1; 
defined $opts{keep} and $opts{keep} ne '' and $have_keep = 1; 

my %is_output; 
FLAG_NUM: 
for (sort { $a<=>$b } keys %flag) {
	my @ta = @{$flag{$_}}; 
	my (%drop_flag, %keep_flag); # Actually these variables are not very useful, but I add them here for future extending functions. 

	if ($have_drop) {
		# logical OR, so any of TRUE will make the dropping [CHK_TTL_DROP] TRUE. 
		CHK_TTL_DROP: 
		for my $exp1 ( split(/;/, $opts{drop}) ) {

			my $should_drop = 1; 
			# logical AND, so any of FALSE will make the $should_drop [CHK_EACH_DROP] FALSE. 
			CHK_EACH_DROP: 
			for my $exp2 ( split(/,/, $exp1) ) {
				$exp2 =~ m/^(\d+)=([01])$/ or die "Failed to parse |$exp2|\n"; 
				my ($coln, $colv) = ($1, $2); 
				$ta[$coln] != $colv and do { $should_drop = 0; last CHK_EACH_DROP; }; 
			}#End for my $exp2 [CHK_EACH_DROP]

			# logical OR, so any of TRUE will make the dropping [CHK_TTL_DROP] TRUE. 
			$should_drop == 1 and do { $drop_flag{$_} = 1; last CHK_TTL_DROP; }; 

		}# End for my $exp1 
		$is_output{$_} = ( defined $drop_flag{$_} and $drop_flag{$_} == 1 ) ? 0 : 1 ; 
	}#End if have_drop

	if ($have_keep) {
		# logical OR, so any of TRUE will make the keeping [CHK_TTL_KEEP] TRUE. 
		CHK_TTL_KEEP: 
		for my $exp1 ( split(/;/, $opts{keep}) ) {

			my $should_keep = 1; 
			# logical AND, so any of FALSE will make the $should_keep [CHK_EACH_KEEP] FALSE. 
			CHK_EACH_KEEP: 
			for my $exp2 ( split(/,/, $exp1) ) {
				$exp2 =~ m/^(\d+)=([01])$/ or die "Failed to parse |$exp2|\n"; 
				my ($coln, $colv) = ($1, $2); 
				$ta[$coln] != $colv and do { $should_keep = 0; last CHK_EACH_KEEP; }; 
			}#End for my $exp2 [CHK_EACH_KEEP]

			# logical OR, so any of TRUE will make the keeping [CHK_TTL_KEEP] TRUE. 
			$should_keep == 1 and do { $keep_flag{$_} = 1; last CHK_TTL_KEEP; }; 

		}# End for my $exp1 [CHK_TTL_KEEP]
		if (defined $keep_flag{$_} and $keep_flag{$_} == 1) {
			defined $is_output{$_} or $is_output{$_} = 1; 
		} else {
			$is_output{$_} = 0; 
		}
	}#End if have_keep 

	if (!defined $is_output{$_}) {
		die "\n\nNo filter parameter accepted! Please check your input for -drop/-keep!\n\n$help_doc\n"; 
#		$is_output{$_} = 1; 
	}
}# End for ( sort { $a<=>$b } keys %flag ) [FLAG_NUM]
###### Set FLAGs OK. 

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

		if ($is_op == 1 and $need_xt_tag == 1) {
			my @ta = split(/\t/, $_); 
			if ($need_xt_tag == 1) {
				my $is_xt_ok = 0; 
				CHK_XT: 
				for my $tb (@ta[11 .. $#ta]) {
					defined $good_xt_tag{$tb} and do { $is_xt_ok = 1; last CHK_XT; }; 
				}
				$is_xt_ok == 1 or $is_op = 0; 
			}#End if CHK_XT; 
		}

		if ($opts{sameRef}) {
			/^(?:[^\t]+\t){6}\=\t/ or $is_op = 0; 
		} elsif ($opts{diffRef}) {
			/^(?:[^\t]+\t){6}\=\t/ and $is_op = 0; 
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

	if ($is_op == 1 and $need_xt_tag == 1) {
		my @ta = split(/\t/, $_); 
		if ($need_xt_tag == 1) {
			my $is_xt_ok = 0; 
			CHK_XT: 
			for my $tb (@ta[11 .. $#ta]) {
				defined $good_xt_tag{$tb} and do { $is_xt_ok = 1; last CHK_XT; }; 
			}
			$is_xt_ok == 1 or $is_op = 0; 
		}#End if CHK_XT; 
	}

	$is_op == 1 and print STDOUT "$_\n"; 

}

