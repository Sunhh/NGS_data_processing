use strict; use warnings;
use FindBin; use lib "$FindBin::Bin/..";
use Test::More;
use File::Temp qw(tempfile);
use fastaSunhh;
my $fa = fastaSunhh->new();

# rcSeq: in-place reverse-complement via scalar ref (tag 'rc')
my $s = 'ATGCaan'; fastaSunhh::rcSeq(\$s, 'rc');
is( $s, 'nttGCAT', 'rcSeq reverse-complements in place' );

# siteList: function-style (\$pattern, \$seq, mode); 'max' = overlapping matches
my ($pat,$seq) = ('[Nn]+', 'ACNNNGTnnA');
is_deeply( [ fastaSunhh::siteList(\$pat, \$seq, 'max') ],
           [[3,5,'NNN'],[4,5,'NN'],[5,5,'N'],[8,9,'nn'],[9,9,'n']],
           'siteList max (overlapping)' );

# bbb2aa: FUNCTION bbb2aa($codon,$tbl_num) -> ($aa,$frame); NOT a method
is_deeply( [ fastaSunhh::bbb2aa('ATG',1) ], ['M',1], 'bbb2aa ATG -> M' );
is( (fastaSunhh::bbb2aa('TAA',1))[0], '*', 'bbb2aa TAA -> * (stop)' );

# get_fasta_seq: OO method, named args; returns (\%rec, $has_next)
my ($fh,$fn) = tempfile(UNLINK=>1);
print {$fh} ">s1 desc one\nACGTACGT\n>s2\nTTTT\n"; close $fh;
my ($rec) = $fa->get_fasta_seq('faFile'=>$fn);
is( $rec->{key}, 's1',       'get_fasta_seq key'  );
is( $rec->{seq}, 'ACGTACGT', 'get_fasta_seq seq'  );

done_testing();
