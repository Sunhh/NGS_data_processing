use strict; use warnings;
use FindBin; use lib "$FindBin::Bin/..";
use Test::More;
use mathSunhh;

# --- pure stats (function-style) ---
is( mathSunhh::min(3,1,2), 1, 'min' );
is( mathSunhh::max(3,1,2), 3, 'max' );
is_deeply( [ mathSunhh::minmax([3,1,2,9,4]) ], [1,9], 'minmax(\@) -> (min,max)' );
is( sprintf("%.4f", mathSunhh::log10(1000)), '3.0000', 'log10(1000)=3' );

# --- combinatorics (function-style, arrayref-first; canonical home for SNP_tbl too) ---
is_deeply( [ mathSunhh::permutations([1,2,3],2) ],
           [[1,2],[1,3],[2,1],[2,3],[3,1],[3,2]], 'permutations(3,2)' );
is_deeply( [ mathSunhh::combinations([1,2,3],2) ],
           [[1,2],[1,3],[2,3]], 'combinations(3,2)' );

# --- param helper (function form: non-object first arg is put back) ---
is_deeply( { mathSunhh::_setHashFromArr('a',1,'b',2) }, {a=>1,b=>2}, '_setHashFromArr' );

# --- interval / location (method-style: shift $self) ---
is_deeply( [ mathSunhh::ovl_region(1,10,5,20) ], [6,[5,10]], 'ovl_region overlap' );
is_deeply( [ mathSunhh::ovl_region(1,4,10,20) ], [0,[]],     'ovl_region disjoint' );
is_deeply( mathSunhh::mergeLocBlk([[1,5],[3,8],[20,25]]),
           [[1,8],[20,25]], 'mergeLocBlk merges touching blocks' );
is_deeply( mathSunhh::repArr(['a','b'], 'times'=>2), ['a','b','a','b'], 'repArr times=2' );

# --- Benjamini-Hochberg FDR (pure-Perl p.adjust(...,"BH")) ---
is_deeply(
  [ map { sprintf "%.6f", $_ } @{ mathSunhh::p_adjust_BH([0.01,0.5,0.9,0.001,0.2]) } ],
  [qw/0.025000 0.625000 0.900000 0.005000 0.333333/],
  "p_adjust_BH matches R p.adjust BH (monotone/cummin)" );
is_deeply(
  [ map { sprintf "%.6f", $_ } @{ mathSunhh::p_adjust_BH([0.5,0.6,0.7]) } ],
  [qw/0.700000 0.700000 0.700000/],
  "p_adjust_BH caps at 1 / enforces monotonicity" );
is_deeply( mathSunhh::p_adjust_BH([0.03]), [0.03], "p_adjust_BH single value" );

# --- chi-square upper-tail p-value (pure-Perl Statistics::Distributions::chisqrprob) ---
is( sprintf("%.4f", mathSunhh::chisqrprob(1, 3.841459)), "0.0500", "chisqrprob df=1 @ 3.8415 = 0.05" );
is( sprintf("%.4f", mathSunhh::chisqrprob(2, 9.210340)), "0.0100", "chisqrprob df=2 @ 9.2103 = 0.01" );
is( sprintf("%.7f", mathSunhh::chisqrprob(3, 5)),        "0.1717971", "chisqrprob df=3 @ 5 == R pchisq" );
is( mathSunhh::chisqrprob(2, 0), 1, "chisqrprob x<=0 -> 1" );
is( mathSunhh::chisqrprob(0, 5), 1, "chisqrprob df<=0 -> 1" );

done_testing();
