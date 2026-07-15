use strict; use warnings;
use FindBin; use lib "$FindBin::Bin/..";
use Test::More;
use SNP_tbl;

# IUPAC degenerate code <-> bases
is( join('', SNP_tbl::dna_d2b('R')), 'AG',   'dna_d2b R -> A,G' );
is( join('', SNP_tbl::dna_d2b('N')), 'ACGT', 'dna_d2b N -> A,C,G,T' );
is( SNP_tbl::dna_b2d('AG'), 'R', 'dna_b2d AG -> R' );
is( SNP_tbl::dna_b2d('GA'), 'R', 'dna_b2d GA -> R (order-insensitive)' );

# permutations/combinations here are thin wrappers delegating to mathSunhh;
# confirm they still work through SNP_tbl's own name.
is_deeply( [ SNP_tbl::combinations([1,2,3],2) ],
           [[1,2],[1,3],[2,3]], 'SNP_tbl::combinations delegates OK' );

done_testing();
