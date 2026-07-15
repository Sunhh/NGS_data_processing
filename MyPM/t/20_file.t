use strict; use warnings;
use FindBin; use lib "$FindBin::Bin/..";
use Test::More;
use File::Temp qw(tempfile);
use fileSunhh;

# NOTE: splitL() CHOMPS its 2nd arg in place -> must pass a modifiable variable.
my $line = "x\ty\tz";
is_deeply( [ fileSunhh::splitL("\t", $line) ], [qw/x y z/], 'splitL on tab' );

my ($c,$b,$d) = ("# comment\n", "   \n", "data\n");
is( fileSunhh::isSkipLine($c), 1, 'isSkipLine: comment -> skip' );
is( fileSunhh::isSkipLine($b), 1, 'isSkipLine: blank -> skip' );
is( fileSunhh::isSkipLine($d), 0, 'isSkipLine: data -> keep' );

# openFH write+read roundtrip
my ($fh,$fn) = tempfile(UNLINK=>1);
my $ofh = fileSunhh::openFH($fn, '>');
print {$ofh} "l1\nl2\n"; close $ofh;
my $ifh = fileSunhh::openFH($fn, '<');
my @got = <$ifh>; close $ifh;
is_deeply( \@got, ["l1\n","l2\n"], 'openFH write then read roundtrip' );

# load_tabFile -> list of arrayrefs (skips comment/blank by default)
my ($fh2,$fn2) = tempfile(UNLINK=>1);
print {$fh2} "# header\na\t1\nb\t2\n"; close $fh2;
is_deeply( [ fileSunhh::load_tabFile($fn2) ], [['a','1'],['b','2']], 'load_tabFile' );

done_testing();
