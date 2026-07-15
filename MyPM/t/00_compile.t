use strict; use warnings;
use FindBin; use lib "$FindBin::Bin/..";
use Test::More;
# Layer 1: every module must load (require = compile + run BEGIN) without error.
my @mods = qw(
  LogInforSunhh mathSunhh fileSunhh fastaSunhh gffSunhh SeqAlnSunhh
  SNP_tbl mcsSunhh plotSunhh ReadInSeqSunhh ReadInAlnSunhh ConfigSunhh
  fromBraker wm97Sunhh
);
require_ok($_) for @mods;
done_testing();
