package ReadInSeqSunhh; 
# Author : Honghe Sun, hs738@cornell.edu
# Date   : 2014-07-29

use strict; 
use warnings; 
use LogInforSunhh;
use Exporter qw(import); 

our @EXPORT = qw(get_fasta_seq); 
our @EXPORT_OK; 

# input a fasta file's handle and a signal show whether it should has a head line;
# return (undef(),undef()) if input is not enough data.
# return (\% , has_get_next)
# 2013-11-01 13:02:52 about \%, {seq, key, head, definition}
# Here we require the fasta header must in a line start with '>', and '>' must be followed by non-blank key information. 
# So there should not be quality lines using '>' as a quality value. 
sub get_fasta_seq {
	my $fh = shift;
	my $has_head = shift; # Input parameter. Value=0 means there is no header information in the input file ($fh). 
	( defined $has_head and $has_head =~ /^0+$/ ) or $has_head = 1; # Value=1 means the first line of input file is header information line. 
	my $has_get = 0;      # Value=1 if I have found a fasta sequence. ; 
	my %backH;
	ref($fh) eq 'GLOB' or ref($fh) eq '' or die "Wrong input!\n"; # 'GLOB' for file handle. 
	if ( $has_head == 1 ) {
		defined ( $backH{head} = readline($fh) ) or return (undef(),undef());
		$backH{head} =~ s/^>//g; chomp $backH{head};
		( $backH{definition} = $backH{head} ) =~ s/^(\S+)// or &stopErr("[Err] Lack of Key in the header: $backH{head}\n"); 
		$backH{key} = $1; 
	}
	my $r = $/; local $/ = "$r>";
	defined ( $backH{seq} = readline($fh) ) or do { &tsmsg("[Err]The last sequence [$backH{head}] is empty, and it is not calculated!\n"); return (undef(),undef()); } ;
	chomp $backH{seq} > length($r) and $has_get = 1; # Here white space will also be treated as part of fasta sequence. 
	local $/ = $r; chomp $backH{seq};                # Remove the last "return" symbol. 
	# check if this sequence is a NULL one. 2008-1-4 13:01:22
	# print "head=+$backH{head}+\nseq=+$backH{seq}+\n";
	while ($backH{seq} =~ s/^>//s) {
		# This should not happen. 
		# The previous sequence record is empty, so it should be skipped. 
		# Now $bachH{seq} is storing the header and sequence of the next sequence. 
		&tsmsg("[Err]Sequence [$backH{head}] is empty, and it is not calculated!\n");
		if ($backH{seq} =~ s/^([^$r]+)(?:$r|$)//s) {
			$backH{head} = $1; 
			( $backH{definition} = $backH{head} ) =~ s/^(\S+)// or &stopErr("[Err] Lack of Key in the header: $backH{head}\n");
			$backH{key} = $1; 
		}
	}#End while 
	# check if this sequence is a NULL one. 2008-1-4 13:01:25
	return (\%backH, $has_get);
}# end sub get_fasta_seq

1; 


