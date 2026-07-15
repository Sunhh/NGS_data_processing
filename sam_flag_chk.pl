#!/usr/bin/perl
# 2013-08-20 I think it is good enough.
# 2014-03-12 Add -NMperc and -maxRdLen to filter "NM" tag. Cause "bwa aln -n float_number" doesn't work properly.
# 2026-07-11 Delegate the FLAG bit-table + decode to MyPM/SeqAlnSunhh (single source
#            of truth: sam_flag_table()/sam_flag_infor()). Behaviour unchanged for
#            flags 0..2047; now also decodes bit 11 (supplementary) up to 4095 and no
#            longer dies on flags >= 2048.
use strict;
use warnings;
use LogInforSunhh;
use SeqAlnSunhh;
use Getopt::Long;

my %opts;
GetOptions(\%opts,
  "help!",
);

my $help_doc = <<HELP;
###################################################################
perl $0 sam_flag_number
# I think it is enough for use.
sam_flag_number <= 4095

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
# 11	0x0800	d	the alignment is supplementary alignment
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

my $flag_num = shift;
(defined $flag_num and $flag_num =~ m/^\d+$/) or &usage();
$flag_num <= 4095 or &stopErr("sam_flag_number should not be larger than 4095.\n");

# Canonical table + per-bit decode both come from MyPM/SeqAlnSunhh.
my $tbl   = &SeqAlnSunhh::sam_flag_table();          # { bit => [bitpos, hex, char, desc] }
my $infor = &SeqAlnSunhh::sam_flag_infor($flag_num); # [ [ bit_value, desc ], ... ]
for (my $i = 0; $i < @$infor; $i++) {
  $infor->[$i][0] == 0 and next;
  print STDOUT join("\t", @{$tbl->{$i}})."\n";
}
