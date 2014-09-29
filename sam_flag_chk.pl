#!/usr/bin/perl
# 2013-08-20 I think it is good enough. 
# 2014-03-12 Add -NMperc and -maxRdLen to filter "NM" tag. Cause "bwa aln -n float_number" doesn't work properly. 
use strict; 
use warnings; 
use LogInforSunhh; 
use Getopt::Long; 

my %opts; 

GetOptions(\%opts, 
	"help!", 
); 

my $help_doc = <<HELP; 
###################################################################
perl $0 sam_flag_number
# I think it is enough for use. 
sam_flag_number <= 2047

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

sub usage {
	print "$help_doc\n"; 
	exit 1; 
}
if (-t and !@ARGV) {
	&usage(); 
} elsif ( $opts{help} ) {
	&usage(); 
}

###### Make FLAG list
my %infor_flag; 
$infor_flag{ 0 } = [qw/0     0x0001  p/, "the read is paired in sequencing"]; 
$infor_flag{ 1 } = [qw/1     0x0002  P/, "the read is mapped in a proper pair"]; 
$infor_flag{ 2 } = [qw/2     0x0004  u/, "the query sequence itself is unmapped"]; 
$infor_flag{ 3 } = [qw/3     0x0008  U/, "the mate is unmapped"]; 
$infor_flag{ 4 } = [qw/4     0x0010  r/, "strand of the query (1 for reverse)"]; 
$infor_flag{ 5 } = [qw/5     0x0020  R/, "strand of the mate"]; 
$infor_flag{ 6 } = [qw/6     0x0040  1/, "the read is the first read in a pair"]; 
$infor_flag{ 7 } = [qw/7     0x0080  2/, "the read is the second read in a pair"]; 
$infor_flag{ 8 } = [qw/8     0x0100  s/, "the alignment is not primary"]; 
$infor_flag{ 9 } = [qw/9     0x0200  f/, "the read fails platform/vendor quality checks"]; 
$infor_flag{ 10} = [qw/10    0x0400  d/, "the read is either a PCR or an optical duplicate"]; 
my $flag_num = shift; 
$flag_num =~ m/^\d+$/ or &usage(); 
$flag_num <= 2047 or &stopErr("sam_flag_number should not be larger than 2047.\n"); 
my $binum = unpack("B32", pack("N", $flag_num)); 
my @flag_arr = reverse( split(//, sprintf("%011d", $binum) ) ); 
for (my $i=0; $i<@flag_arr; $i++) {
	$flag_arr[$i] == 0 and next; 
	print STDOUT join("\t", @{$infor_flag{$i}})."\n"; 
}


