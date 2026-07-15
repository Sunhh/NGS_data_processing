use strict; use warnings;
use FindBin; use lib "$FindBin::Bin/..";
use Test::More;
use File::Temp qw(tempfile);
use ReadInSeqSunhh;

# Legacy positional reader: get_fasta_seq($fh, $has_head) -> (\%rec, $has_next)
my ($fh,$fn) = tempfile(UNLINK=>1);
print {$fh} ">s1 desc one\nACGTACGT\n>s2\nTTTT\n"; close $fh;
open my $r,'<',$fn or die $!;
my ($rec,$got) = ReadInSeqSunhh::get_fasta_seq($r, 1);
close $r;
is( $rec->{key},        's1',        'ReadInSeq key' );
is( $rec->{seq},        'ACGTACGT',  'ReadInSeq seq' );
is( $rec->{definition}, ' desc one', 'ReadInSeq definition (space-prefixed)' );
is( $got, 1, 'ReadInSeq has_next flag' );

done_testing();
