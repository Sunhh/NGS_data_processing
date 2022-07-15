#!/usr/bin/perl
use strict;
use warnings;
use fileSunhh;
use LogInforSunhh;
use Getopt::Long;
my %opts;
GetOptions(\%opts,
  "qry_tbl:s", "ref_tbl:s",
  "out_prefix:s",
  "help!"
);

my $htxt = <<HH;
##########################################################################################
# perl $0 -qry_tbl nucdiff_query_struct.gff -ref_tbl nucdiff_ref_struct.gff -out_prefix shared
##########################################################################################
HH

$opts{'out_prefix'} //= "shared";

for (qw/qry_tbl ref_tbl/) {
  defined $opts{$_} or &LogInforSunhh::usage($htxt);
}

my $q_tbl = &load_tbl($opts{'qry_tbl'});
my $r_tbl = &load_tbl($opts{'ref_tbl'});
my (%qH, %rH);
for (@$q_tbl) { $qH{$_->[0]} = 1; }
for (@$r_tbl) { defined $qH{$_->[0]} and $rH{$_->[0]} = 1; }

&out_tbl("$opts{'out_prefix'}.qry", $q_tbl, \%rH);
&out_tbl("$opts{'out_prefix'}.ref", $r_tbl, \%rH);

sub out_tbl {
  my ($ofn, $tbl_arr, $id_hash) = @_;
  my $ofh = &openFH($ofn, '>');
  for (@$tbl_arr) {
    defined $id_hash->{$_->[0]} or next;
    print {$ofh} "$_->[2]\n";
  }
  close($ofh);
  return();
}# out_tbl()

sub load_tbl {
  my ($f) = @_;
  my $ifh = &openFH($f, '<');
  my @back;
  while (<$ifh>) {
    m!^\s*(#|$)! and next;
    chomp;
    # my @ta=split(/\t/, $_);
    m!(?:^|\s|;)ID=([^\s.;]+)(\.\d+)?;! or &stopErr("[Err] Bad 1 - ID not found: [$_]\n");
    push(@back, [$1, $2, $_]);
  }
  close($ifh);
  return(\@back);
}# load_tbl()

# C31_Chr01       NucDiff_v2.0    SO:0001873      9464    14349   .       .       .       ID=SV_1;Name=translocation-overlap;overlap_len=4886;ref_sequence_1=C39_Chr16;blk_1_query=1-14349;blk_1_ref=2139435-2153819;blk_1_query_len=14349;blk_1_ref_len=14385;blk_1_st_query=1;blk_1_st_ref=2139435;blk_1_end_query=14349;blk_1_end_ref=2153819;ref_sequence_2=C39_Chr01;blk_2_query=9464-11510801;blk_2_ref=1-12001526;blk_2_query_len=11501338;blk_2_ref_len=12001526;blk_2_st_query=9464;blk_2_st_ref=1;blk_2_end_query=11510801;blk_2_end_ref=12001526;color=#A0A0A0
# C31_Chr01       NucDiff_v2.0    SO:1000002      184573  185273  .       .       .       ID=SV_2;Name=substitution;subst_len=701;query_dir=1;ref_sequence=C39_Chr01;ref_coord=174594-175294;color=#42C042
# C31_Chr01       NucDiff_v2.0    SO:0000159      185273  185273  .       .       .       ID=SV_3;Name=deletion;del_len=540;query_dir=1;ref_sequence=C39_Chr01;ref_coord=175295-175834;color=#0000EE
# C31_Chr01       NucDiff_v2.0    SO:0000667      216258  216463  .       .       .       ID=SV_4;Name=insertion;ins_len=206;query_dir=1;ref_sequence=C39_Chr01;ref_coord=207278;color=#EE0000
# C31_Chr01       NucDiff_v2.0    SO:1000035      216464  216465  .       .       .       ID=SV_5;Name=duplication;ins_len=2;query_dir=1;ref_sequ\ence=C39_Chr01;ref_coord=207278;query_repeated_region=216256-216257;color=#EE0000
# C31_Chr01       NucDiff_v2.0    SO:0000159      221576  221576  .       .       .       ID=SV_6;Name=deletion;del_len=212;query_dir=1;ref_sequence=C39_Chr01;ref_coord=212479-212690;color=#0000EE
