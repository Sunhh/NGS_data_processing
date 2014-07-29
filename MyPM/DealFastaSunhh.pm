package DealFastaSunhh; 

#Author: Honghe Sun, hs738@cornell.edu, 2014-07-15
# Date   : 2014-07-29

use strict;
use warnings;
use LogInforSunhh;
use Exporter qw(import);

our @EXPORT = qw(siteList);
our @EXPORT_OK;

#  Input : (\$expr_pattern, \$sequence_to_search, $check_mode)
#            $check_mode : 'min' - give the least number of matched patterns without overlap. 
#                          'max' - give the most number of matched patterns including overlap. 
# Return : ( [match_Start_1, match_End_1, match_Seq_1], [match_Start_2, match_End_2, match_Seq_2], ... )
sub siteList ($$$) {
	my ($siteR, $refR, $modeChk) = @_; 
	( defined $siteR and ref($siteR) eq 'SCALAR' ) or &stopErr("[Err] Input src_pattern wrong!\n"); 
	$$siteR eq '' and &stopErr("[Err] Input srch_pattern is empty. Exit!\n"); 
	my $qrSite = qr/$$siteR/s; 
	( defined $refR and ref($refR) eq 'SCALAR' ) or &stopErr("[Err] Input src_pattern wrong!\n"); 
	my $Is_min = 1; 
	if ( defined $modeChk ) {
		$modeChk = lc($modeChk); 
		if ( $modeChk eq 'max' ) {
			$Is_min = 0; 
		} elsif ( $modeChk eq 'min' ) {
			$Is_min = 1; 
		} else {
			&tsmsg("[Err] Third parameter of sub siteList should be 'Min' or 'Max', instead of [$modeChk].\n"); 
			&tsmsg("[Err] Using 'Min'.\n"); 
		}
	}

	my @posList = ();
	pos($$refR) = 0;
	while ($$refR =~ m/\G(?:.*?)($qrSite)/gs) {
		push ( @posList, [ $-[1]+1 , $+[1], $1 ] );
		$Is_min?(pos($$refR) = $+[1]):(pos($$refR) = $-[1]+1);
	}#End while 
	return @posList;
}# end siteList subroutine. 2013-10-30

1; 
