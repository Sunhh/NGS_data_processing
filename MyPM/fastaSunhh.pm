package fastaSunhh; 

#Author: Honghe Sun, hs738@cornell.edu, 2014-07-15
# Date   : 2014-07-29

use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use mathSunhh; 
use Exporter qw(import);
our @EXPORT = qw(siteList);
our @EXPORT_OK;

##########################################################
#########  Basic settings.
##########################################################
my $mathObj = mathSunhh->new(); 
my %str2num = qw(
  +      1
  -     -1
  1      1
  -1    -1
  plus   1
  minus -1
  0     -1
); 
my %IUPAC_b2d; # bases to degenerate;
my %IUPAC_d2b; # degenerate to bases;
{
my @ta = qw(
AA   CC   GG   TT   UU
WAT  SCG  MAC  KGT  RAG   YCT
BCGT DAGT HACT VACG -ACGT NACGT
);
	for my $tb (@ta) {
		my @bb = split(//, $tb);
		my $dbase = shift(@bb);
		$IUPAC_b2d{$dbase} = $dbase;
		for (&mathSunhh::permutations(\@bb, scalar(@bb))) {
			$IUPAC_b2d{join('', @$_)} = $dbase;
		}
		for (&mathSunhh::combinations(\@bb, scalar(@bb))) {
			$IUPAC_d2b{$dbase} = [sort @$_];
		}
	}
}# End of IUPAC_xxx;

my %codon_tbl; 


##########################################################
#########  Methods.
##########################################################

=head1 fastaSunhh->new()

Function     : Open an object of fastaSunhh

=cut
sub new {
	my $class = shift; 
	my $self = {}; 
	bless $self, $class; 
	
	$self->_initialize(@_); 
	
	return $self; 
}#sub new

sub _initialize {
	my $self = shift; 
	
	my %parm = $mathObj->_setHashFromArr(@_); 
	for my $k ( keys %parm ) {
		$self->{$k} = $parm{$k}; 
	}
	
	return; 
}#sub _initialize

=head2 save_seq_to_hash( 'faFh'=>$faFileHandle, 'faFile'=>$faFileName, 'has_head'=>[1/0] )

Required     : 'faFh' , or 'faFile' when 'faFh' missed. 

Function     : Save all sequences' information in hash

'has_head' = 1 by default. 

Return       : (\%key_to_infor)

 \%key_to_infor : {key} => value

   {$seq_ID}      => \%seq_record_hash; 

     \%seq_record_hash : {key} => value 

       {'head'}       => the line heading with '>', but with heading '>' removed; 
       {'key'}        => sequence ID ; 
       {'definition'} => Same to 'head' except that 'key' is removed. 
       {'seq'}        => sequence information in the same format of input. 
       {'has_get'}    => Tell if there is another sequence in the next. 
       {'FH'}         => File handle being used. 
       {'Order'}      => Tell the order number in the fasta file. This is useful when sorting sequences by input order. 

=cut
sub save_seq_to_hash {
	my $self = shift; 
	my %parm = $mathObj->_setHashFromArr(@_); 
	my $fh = $parm{'faFh'} // &fileSunhh::openFH($parm{'faFile'}, '<') // &stopErr("[Err] No valid file handle in get_fasta_seq()\n"); 
	$parm{'has_head'} //= 1; 
	$parm{'has_head'} =~ m/^0+$/ or $parm{'has_head'} = 1; 
	ref($fh) eq 'GLOB' or ref($fh) eq '' or &stopErr("[Err] File handle wrong.\n"); 
	
	my %backH; 
	my $seqOrderNum = 0; 
	# For the first sequence; 
	{
		my $relHR = $self->get_fasta_seq( 'faFh'=>$fh, 'has_head'=>$parm{'has_head'} ); 
		if ( $relHR->{'has_get'} == 1 ) {
			$backH{$relHR->{'key'}} = $relHR; 
			++$seqOrderNum; 
			$backH{$relHR->{'key'}}{'Order'} = $seqOrderNum; 
		}
	}
	for (my $relHR = $self->get_fasta_seq( 'faFh'=>$fh ); $relHR->{'has_get'} == 1; $relHR = $self->get_fasta_seq( 'faFh'=>$fh ) ) {
		if ( defined $backH{ $relHR->{'key'} } ) {
			&tsmsg("[Err] Key [$relHR->{'key'}] repeated, and I choose the first one.\n"); 
			next; 
		}
		++$seqOrderNum; 
		$relHR->{'Order'} = $seqOrderNum; 
		$backH{ $relHR->{'key'} } = $relHR; 
	}
	
	return (\%backH); 
}# sub save_seq_to_hash () 

=head2 get_fasta_seq('faFh'=>$faFileHandle, 'faFile'=>$faFileName, 'has_head'=>[1/0])

Required     : 'faFh' , or 'faFile' when 'faFh' missed. 

Function     : Input a fasta file's handle and a signal showing whether it has a head line, 
                 then it returns a hash storing one sequence record once. 
Return       : (\%seq_record_hash)
 \%seq_record_hash : {key} => value 
   {'head'}       => the line heading with '>', but with heading '>' removed; 
   {'key'}        => sequence ID ; 
   {'definition'} => Same to 'head' except that 'key' is removed. 
   {'seq'}        => sequence information in the same format of input. 
   {'has_get'}   => Tell if there is another sequence in the next. 
   {'FH'}         => File handle being used. 

Example      : 
 
my $fs = fastaSunhh->new(); 
my $inFa = "tt.fa"; 
my $inFh = &openFH($inFa, "<"); # Input file handle 
# for (my $relHR = $fs->get_fasta_seq('faFh'=>$inFh); $relHR->{'has_get'} == 1; $relHR = $fs->get_fasta_seq('faFh'=>$inFh) ) {
for (my $relHR = $fs->get_fasta_seq('faFile'=>$inFa); $relHR->{'has_get'} == 1; $relHR = $fs->get_fasta_seq('faFh'=>$relHR->{'FH'}) ) {
	my $seq_id  = $relHR->{'key'}; 
	my $seq_seq = $relHR->{'seq'}; $seq_seq =~ s!\s!!g; 
	my $seq_head = $relHR->{'head'}; 
	my $seq_def  = $relHR->{'definition'}; 
}
# All sequences in $inFh have been processed. 
=cut
sub get_fasta_seq {
	my $self = shift; 
	my %parm = $mathObj->_setHashFromArr(@_); 
	my $fh = $parm{'faFh'} // &fileSunhh::openFH($parm{'faFile'}, '<') // &stopErr("[Err] No valid file handle in get_fasta_seq()\n"); 
	$parm{'has_head'} //= 1; 
	$parm{'has_head'} =~ m/^0+$/ or $parm{'has_head'} = 1; 
	ref($fh) eq 'GLOB' or ref($fh) eq '' or &stopErr("[Err] File handle wrong.\n"); 
	
	my %backH; 
	$backH{'has_get'} = 0; 
	$backH{'FH'} = $fh; 
	# Retrieve head information. 
	if ( $parm{'has_head'} == 1 ) {
		defined ( $backH{'head'} = readline($fh) ) or return (\%backH); 
		$backH{'head'} =~ s!^>!!g; chomp($backH{'head'}); 
		( $backH{'definition'} = $backH{'head'} ) =~ s!^(\S+)!!; 
		$backH{'key'} = $1; 
		$backH{'has_get'} = 1; 
	}
	# Retrieve sequence information. 
	my $r = $/; local $/ = "$r>"; 
	unless ( defined ( $backH{'seq'} = readline($fh) ) ) {
		$parm{'has_head'} and &tsmsg("[Wrn] The last sequence [$backH{'head'}] is empty, so it is not calculated!\n");
		return (\%backH); 
	}
	( chomp($backH{'seq'}) ) > length($r) and $backH{'has_get'} = 1; 
	local $/ = $r; chomp($backH{'seq'}); 
	
	# Check if this sequence is a NULL one. 
	# When the sequence is a NULL one, the $backH{'seq'} will read the next whole sequence including the header line. 
	# So I want to fix this error. 
	while ($backH{'seq'} =~ s!^>!!gs) {
		&tsmsg("[Wrn] Sequence [$backH{'head'}] is empty, and it is skipped!\n"); 
		# if ( $backH{'seq'} =~ s!^([^$r]+)(?:$r|$)!!s ) {
		if ( $backH{'seq'} =~ s!^(.+)(?:$r|$)!!s ) { # I prefer to edit this line because I am afraid of Nt/Nr database fasta. 
			$backH{'definition'} = $backH{'head'} = $1; 
			chomp( $backH{'definition'} ); chomp( $backH{'head'} ); 
			$backH{'definition'} =~ s!^(\S+)!!; $backH{'key'} = $1; 
		}
	}
	
	return (\%backH); 
}# sub get_fasta_seq ()

##########################################################
#########  Sub-functions.
##########################################################

=head1 rcSeq( \$seq_str, $deal_tag ) 

Input       : 
 \$seq_string   : A string storing sequence without '\s' (continuous)
 $deal_tag      : What to do with this $seq_string
                  'r' => reverse , 
                  'c' => complemented , 
                  'rc' => reverse and complemented [Default] 
Function    : Reverse and complementing a DNA string. 
Return      : No return, just edit $seq_string . 

=cut
sub rcSeq {
	my $seq_r = shift; 
	my $tag = shift; $tag //= 'rc'; 
	my ($Is_r, $Is_c) = (0, 0); 
	$tag =~ m/r/i and $Is_r = 1; 
	$tag =~ m/c/i and $Is_c = 1; 
	!$Is_r and !$Is_c and &stopErr("[Err] Why do you call me rcSeq() as deal_tag is [$tag]?\n"); 
	$Is_r and $$seq_r = reverse( $$seq_r ); 
	$Is_c and $$seq_r =~ tr/acgturykmbvdhACGTURYKMBVDHwWsSnN/tgcaayrmkvbhdTGCAAYRMKVBHDwWsSnN/; 
	return; 
}# rcSeq() 


=head1 siteList( \$expr_pattern, \$sequence_to_search, $check_mode )

Input       : 
 \$expr_pattern : Something like '[Nn]+'
 \$sequence     : 
 $check_mode    : 'min' - give the least number of matched patterns without overlap. 
                  'max' - give the most number of matched patterns including overlap. 
Function    : Search $expr_pattern in $sequence. 
Return      : ( [match_Start_1, match_End_1, match_Seq_1], [match_Start_2, match_End_2, match_Seq_2], ... )
=cut
sub siteList {
	my $self = shift; 
	my $siteR = ( ref($self) eq 'fastaSunhh' ) ? shift : $self ; 
	my ($refR, $modeChk) = @_; 
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

=head1 chop_seq( 'seq'=>$seq_str, 'len'=>$length_of_piece, 'step'=>$NextStart_minus_CurrStart, 'min'=>$MinLenOfPiece )

Return       : ( \@piece_strings )
 @piece_strings = ( [$piece_string, $piece_start, $piece_end], [$piece_string, $piece_start, $piece_end], ... )

=cut
sub chop_seq {
	my %parm = $mathObj->_setHashFromArr(@_); 
	defined $parm{'seq'} or &stopErr( "[Err] 'seq' in chop_seq() is not defined.\n" ); 
	$parm{'len'}  //= 100; 
	$parm{'step'} //= $parm{'len'}; 
	$parm{'min'}  //= 0; 

	$parm{'seq'} =~ s!\s!!g; 
	my $seqLen = length( $parm{'seq'} ); 
	
	my @back; 
	for (my $i=1; ($i-1)*$parm{'step'} + 1 < $seqLen ; $i++) {
		my $s = ($i-1) * $parm{'step'} + 1; 
		my $e = $s + $parm{'len'} - 1; 
		$e > $seqLen and $e = $seqLen; 
		$e-$s+1 >= $parm{'min'} or next; 
		my $sub_seq = substr( $parm{'seq'}, $s-1, $e-$s+1 ); 
		push(@back, [$sub_seq, $s, $e]); 
		$e >= $seqLen and last; 
	}
	return (\@back); 
}# chop_seq() 

=head1 _possible_bbb( 'ATB' )

Return    : (@possible_AAAs)

  _possible_bbb( 'ABT' ) returns ( 'ACT', 'AGT', 'ATT' ); 

=cut
sub _possible_bbb {
	my ($bbb) = @_; 
	my @back; 
	my @bbb_arr = split(//, $bbb); 
	my @bN; 
	for (my $i=0; $i<3; $i++) {
		$bbb_arr[$i] = uc($bbb_arr[$i]); 
		if ( defined $IUPAC_d2b{$bbb_arr[$i]} ) {
			$bN[$i] = [ @{$IUPAC_d2b{$bbb_arr[$i]}} ]; 
		} else {
			&stopErr("[Err] Unsupported base [$bbb_arr[$i]]\n"); 
		}
	}
	for my $b1 (@{$bN[0]}) {
		for my $b2 (@{$bN[1]}) {
			for my $b3 (@{$bN[2]}) {
				push(@back, "$b1$b2$b3"); 
			}
		}
	}
	return(@back); 
}# _possible_bbb() 

=head1 setup_codon_tbl() 

Return : ( \%codon_tbl )

Setup %codon_tbl which is used only in PM; 

$codon_tbl{'src'}{'transl_tbl_\d\d'} : text from NCBI http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi 
$codon_tbl{'id2num'}{$some_ID} = $tbl_num_wo0 
$codon_tbl{'tbl2inf'}{$tbl_num_wo0}{'bbb2aa'}    = { 'TTT' => F , ... }
$codon_tbl{'tbl2inf'}{$tbl_num_wo0}{'bbb2start'} = { 'ATG' => 1 , ... }
$codon_tbl{'has_setup'} = 1; 

=cut
sub setup_codon_tbl {
	defined $codon_tbl{'has_setup'} and $codon_tbl{'has_setup'} == 1 and return (\%codon_tbl); 
# &tsmsg("[Msg] Setting codon tables.\n"); 
$codon_tbl{'src'}{'transl_tbl_01'} = <<'TTT'; 
    AAs  = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  Starts = ---M---------------M---------------M----------------------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
TTT

$codon_tbl{'src'}{'transl_tbl_02'} = <<'TTT'; 
    AAs  = FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG
  Starts = --------------------------------MMMM---------------M------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
TTT
$codon_tbl{'src'}{'transl_tbl_03'} = <<'TTT'; 
    AAs  = FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  Starts = ----------------------------------MM----------------------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
TTT
$codon_tbl{'src'}{'transl_tbl_04'} = <<'TTT'; 
    AAs  = FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  Starts = --MM---------------M------------MMMM---------------M------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
TTT
$codon_tbl{'src'}{'transl_tbl_05'} = <<'TTT'; 
    AAs  = FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG
  Starts = ---M----------------------------MMMM---------------M------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
TTT
$codon_tbl{'src'}{'transl_tbl_06'} = <<'TTT'; 
    AAs  = FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  Starts = -----------------------------------M----------------------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
TTT
$codon_tbl{'src'}{'transl_tbl_07'} = <<'TTT'; 
TTT
$codon_tbl{'src'}{'transl_tbl_08'} = <<'TTT'; 
TTT
$codon_tbl{'src'}{'transl_tbl_09'} = <<'TTT'; 
    AAs  = FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG
  Starts = -----------------------------------M---------------M------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
TTT
$codon_tbl{'src'}{'transl_tbl_10'} = <<'TTT'; 
    AAs  = FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  Starts = -----------------------------------M----------------------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
TTT
$codon_tbl{'src'}{'transl_tbl_11'} = <<'TTT'; 
    AAs  = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  Starts = ---M---------------M------------MMMM---------------M------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
TTT
$codon_tbl{'src'}{'transl_tbl_12'} = <<'TTT'; 
    AAs  = FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  Starts = -------------------M---------------M----------------------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
TTT
$codon_tbl{'src'}{'transl_tbl_13'} = <<'TTT'; 
    AAs  = FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG
  Starts = ---M------------------------------MM---------------M------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
TTT
$codon_tbl{'src'}{'transl_tbl_14'} = <<'TTT'; 
    AAs  = FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG
  Starts = -----------------------------------M----------------------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
TTT
$codon_tbl{'src'}{'transl_tbl_15'} = <<'TTT'; 
TTT
$codon_tbl{'src'}{'transl_tbl_16'} = <<'TTT'; 
    AAs  = FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  Starts = -----------------------------------M----------------------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
TTT
$codon_tbl{'src'}{'transl_tbl_21'} = <<'TTT'; 
    AAs  = FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG
  Starts = -----------------------------------M---------------M------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
TTT
$codon_tbl{'src'}{'transl_tbl_22'} = <<'TTT'; 
    AAs  = FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  Starts = -----------------------------------M----------------------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
TTT
$codon_tbl{'src'}{'transl_tbl_23'} = <<'TTT'; 
    AAs  = FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  Starts = --------------------------------M--M---------------M------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
TTT
$codon_tbl{'src'}{'transl_tbl_24'} = <<'TTT'; 
    AAs  = FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG
  Starts = ---M---------------M---------------M---------------M------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
TTT
$codon_tbl{'src'}{'transl_tbl_25'} = <<'TTT'; 
    AAs  = FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  Starts = ---M-------------------------------M---------------M------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
TTT
$codon_tbl{'src'}{'transl_tbl_26'} = <<'TTT'; 
    AAs  = FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  Starts = -------------------M---------------M----------------------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
TTT

	for my $tbl_ID (sort keys %{$codon_tbl{'src'}}) {
		# &tsmsg("[Msg] Setting codon table '$tbl_ID'\n"); 
		( defined $codon_tbl{'src'}{$tbl_ID} and $codon_tbl{'src'}{$tbl_ID} !~ m!^\s*$! ) or next; 
		$tbl_ID =~ m/^transl_tbl_(\d+)$/ or &stopErr("[Err] Bad transl_tbl ID [$tbl_ID]\n"); 
		my $tbl_num = $1; 
		my $tbl_num_wo0 = $tbl_num; $tbl_num_wo0 =~ s!^0+!!; 
		$codon_tbl{'id2num'}{$tbl_ID}      = $tbl_num_wo0; 
		$codon_tbl{'id2num'}{$tbl_num}     = $tbl_num_wo0; 
		$codon_tbl{'id2num'}{$tbl_num_wo0} = $tbl_num_wo0; 
		# $codon_tbl{'tbl2inf'}{$tbl_num_wo0}{'bbb2aa'}
		my @ll = split(/\n/, $codon_tbl{'src'}{$tbl_ID}); 
		@ll == 5 or &stopErr("[Err] Bad transl_tbl text [$codon_tbl{'src'}{$tbl_ID}]\n"); 
		$ll[0] =~ s!^\s*AAs\s*=\s*!!i    or &stopErr("[Err] Bad line [$ll[0]]\n"); $ll[0] = [split(//, $ll[0])]; 
		$ll[1] =~ s!^\s*Starts\s*=\s*!!i or &stopErr("[Err] Bad line [$ll[1]]\n"); $ll[1] = [split(//, $ll[1])]; 
		$ll[2] =~ s!^\s*Base1\s*=\s*!!i  or &stopErr("[Err] Bad line [$ll[2]]\n"); $ll[2] = [split(//, $ll[2])]; 
		$ll[3] =~ s!^\s*Base2\s*=\s*!!i  or &stopErr("[Err] Bad line [$ll[3]]\n"); $ll[3] = [split(//, $ll[3])]; 
		$ll[4] =~ s!^\s*Base3\s*=\s*!!i  or &stopErr("[Err] Bad line [$ll[4]]\n"); $ll[4] = [split(//, $ll[4])]; 
		for (my $i=0; $i<@{$ll[0]}; $i++) {
			my $bbb = "$ll[2][$i]$ll[3][$i]$ll[4][$i]"; 
			$codon_tbl{'tbl2inf'}{$tbl_num_wo0}{'bbb2aa'}{$bbb} = $ll[0][$i]; 
			$ll[1][$i] eq 'M' and $codon_tbl{'tbl2inf'}{$tbl_num_wo0}{'bbb2start'}{$bbb} = 1; 
		}
		# Setup with one 'ATGC' codon 
		for my $i ( 0 .. 2 ) {
			for my $tb_i (qw/A T G C/) {
				my @bbb_arr = ('N','N','N'); 
				$bbb_arr[$i] = $tb_i; 
				my $bbb = join('', @bbb_arr); 
				my @all_bbb = &_possible_bbb($bbb); 
				my %aa; 
				my @aa_k; 
				for my $curr_bbb ( @all_bbb ) {
					defined $codon_tbl{'tbl2inf'}{$tbl_num_wo0}{'bbb2aa'}{$curr_bbb} or &stopErr("[Err] Unknown bbb [$curr_bbb]\n"); 
					defined $aa{ $codon_tbl{'tbl2inf'}{$tbl_num_wo0}{'bbb2aa'}{$curr_bbb} } or push(@aa_k, $codon_tbl{'tbl2inf'}{$tbl_num_wo0}{'bbb2aa'}{$curr_bbb}); 
					$aa{ $codon_tbl{'tbl2inf'}{$tbl_num_wo0}{'bbb2aa'}{$curr_bbb} } ++; 
					scalar(@aa_k) > 1 and last; 
				}
				$codon_tbl{'tbl2inf'}{$tbl_num_wo0}{'bbb2aa'}{$bbb} = (scalar(@aa_k) == 1) ? $aa_k[0] : 'X' ; 
			}
		}
		# Setup with two 'ATGC' codon 
		for (my $i=0; $i<3; $i++) {
			my $bbb = 'NNN'; 
			substr($bbb, $i, 1) = 'A'; 
			my @all_bbb = map { substr($_, $i, 1) = 'N'; $_; } &_possible_bbb($bbb); 
			for my $curr_chk (@all_bbb) {
				my @all_bbb_chk = &_possible_bbb($curr_chk); 
				my (%aa, @aa_k); 
				for my $curr_bbb ( @all_bbb_chk ) {
					defined $codon_tbl{'tbl2inf'}{$tbl_num_wo0}{'bbb2aa'}{$curr_bbb} or &stopErr("[Err] Unknown bbb [$curr_bbb]\n"); 
					defined $aa{ $codon_tbl{'tbl2inf'}{$tbl_num_wo0}{'bbb2aa'}{$curr_bbb} } or push(@aa_k, $codon_tbl{'tbl2inf'}{$tbl_num_wo0}{'bbb2aa'}{$curr_bbb}); 
					$aa{ $codon_tbl{'tbl2inf'}{$tbl_num_wo0}{'bbb2aa'}{$curr_bbb} } ++; 
					scalar(@aa_k) > 1 and last; 
				}
				$codon_tbl{'tbl2inf'}{$tbl_num_wo0}{'bbb2aa'}{$curr_chk} = (scalar(@aa_k) == 1) ? $aa_k[0] : 'X' ; 
			}
		}
		# Setup with no 'ATGC' codon 
		$codon_tbl{'tbl2inf'}{$tbl_num_wo0}{'bbb2aa'}{'NNN'} = 'X'; 
	}
	$codon_tbl{'has_setup'} = 1; 
	# &tsmsg("[Msg] Codon table setup finished.\n"); 

	return( \%codon_tbl ); 
}# setup_codon_tbl () 

=head1 bbb2aa ( 'ATG', $tbl_num_wo0 )

Return : ( $aa, $is_start )

=cut
sub bbb2aa {
	my ( $bbb, $tbl_num ) = @_; 
	$tbl_num //= 1; 
	&setup_codon_tbl(); 
	defined $codon_tbl{'id2num'}{$tbl_num} or &stopErr("[Err] Bad table number [$tbl_num]\n"); 
	$tbl_num = $codon_tbl{'id2num'}{$tbl_num}; 
	defined $codon_tbl{'tbl2inf'}{$tbl_num} or &stopErr("[Err] Failed to find table for [$tbl_num]\n"); 
	$bbb = uc($bbb); 
	my $b_start = 0; 
	if ( $bbb =~ m/^\-\-\-$/ ) {
		return( '-', 0 ); 
	}
	if ( $bbb =~ m/^[ATGC]{3}$/ ) {
		defined $codon_tbl{'tbl2inf'}{$tbl_num}{'bbb2aa'}{$bbb} or &stopErr("[Err] Failed to find aa for bbb [$bbb]\n"); 
		defined $codon_tbl{'tbl2inf'}{$tbl_num}{'bbb2start'}{$bbb} and $codon_tbl{'tbl2inf'}{$tbl_num}{'bbb2start'}{$bbb} and $b_start = 1; 
		return( $codon_tbl{'tbl2inf'}{$tbl_num}{'bbb2aa'}{$bbb} , $b_start ); 
	}
	if ( defined $codon_tbl{'tbl2inf'}{$tbl_num}{'bbb2aa'}{$bbb} ) {
		defined $codon_tbl{'tbl2inf'}{$tbl_num}{'bbb2start'}{$bbb} and $codon_tbl{'tbl2inf'}{$tbl_num}{'bbb2start'}{$bbb} and $b_start = 1; 
		return( $codon_tbl{'tbl2inf'}{$tbl_num}{'bbb2aa'}{$bbb} , $b_start ); 
	}
	
	&tsmsg("[Wrn] Not able to deal with [$bbb] yet.\n"); 
	return( 'X' , 0 ); 
	return; 
}# bbb2aa

sub aa2cds {
	my ( $cds_fn, $prot_fn, $trans_tbl ) = @_; 
	$trans_tbl //= 1; 
	&setup_codon_tbl(); 
	my %stop_codon; 
	for my $tc (qw/TAA TGA TAG/) {
		$stop_codon{$tc} = 1; 
	}
	my %cds_seq; 
	my $fh_cds = &openFH( $cds_fn , '<' ); 
	while (<$fh_cds>) {
		if ( m/^\s*\>\s*(\S+)/ ) {
			$cds_seq{'tk'} = $1; 
			$cds_seq{ 'seq' }{ $cds_seq{'tk'} } = ''; 
		} else {
			$cds_seq{ 'seq' }{ $cds_seq{'tk'} } .= $_; 
		}
	}
	close($fh_cds); 
	delete($cds_seq{'tk'}); 
	for my $tk (keys %{$cds_seq{'seq'}}) {
		$cds_seq{'seq'}{$tk} =~ s![\s\-]!!g; 
		$cds_seq{'seq'}{$tk} = uc($cds_seq{'seq'}{$tk}); 
		$cds_seq{'len'}{$tk} = length($cds_seq{'seq'}{$tk}); 
	}
	my %prot_seq; 
	$prot_seq{'nn'} = 0; 
	my $fh_prot = &openFH( $prot_fn , '<' ); 
	while (<$fh_prot>) {
		if ( m/^\s*\>\s*(\S+)/ ) {
			$prot_seq{'tk'} = $1; 
			$prot_seq{'seq'}{$prot_seq{'tk'}} = ''; 
			$prot_seq{'nn'} ++; 
			$prot_seq{'ord'}{ $prot_seq{'tk'} } = $prot_seq{'nn'}; 
		} else {
			$prot_seq{'seq'}{$prot_seq{'tk'}} .= $_; 
		}
	}
	delete($prot_seq{'tk'}); 
	delete($prot_seq{'nn'}); 

	for my $tk (sort { $prot_seq{'ord'}{$a} <=> $prot_seq{'ord'}{$b} } keys %{$prot_seq{'seq'}}) {
		unless ( defined $cds_seq{'seq'}{ $tk } ) {
			&tsmsg("[Err] Skip prot_seq [$tk] because of lack of cds seq.\n"); 
			next; 
		}
		$prot_seq{'seq'}{$tk} =~ s!\s!!g; 
		$prot_seq{'seq'}{$tk} = uc($prot_seq{'seq'}{$tk}); 

		my $new_cds_seq = &aa2cds_1seq( $cds_seq{'seq'}{$tk}, $prot_seq{'seq'}{$tk}, $trans_tbl ); 
		unless ( defined $new_cds_seq ) {
			&tsmsg("[Err] Skip prot_seq [$tk] because of incompatibility between prot and cds sequences.\n"); 
			next; 
		}

		$new_cds_seq =~ s!(.{60})!$1\n!g; 
		chomp($new_cds_seq); 
		print STDOUT ">$tk\n$new_cds_seq\n"; 
	}
}# aa2cds() 

=head1 aa2cds_1seq( $cds_seq, $prot_seq, $translation_table_number )

Return    : ( $cds_seq_gapped_by_protein_positions )

  If returning "undef()", there must be some error in it. 

=cut
sub aa2cds_1seq {
	# Input sequences should be upper case, and no blank or return symbol! 
	my ( $cseq, $pseq, $tblN ) = @_; 
	$tblN //= 1; 
	&setup_codon_tbl(); 
	my $aa2cds_cseq; 

	my @aa = split(//, $pseq); 
	my $aa_pos = -1; 
	for (my $j=0; $j<@aa; $j++) {
		if ( $aa[$j] eq '-' or $aa[$j] eq ' ' ) {
			$aa2cds_cseq .= $aa[$j] x 3; 
		} else {
			$aa_pos ++; 
			my $bbb = substr( $cseq, $aa_pos*3, 3 ); 
			my $bbb_len = length($bbb); 
			$bbb_len > 0 or do { &tsmsg("[Err] Failed to extract bbb at [" . ($aa_pos*3+1) . "]: $cseq\n"); return(undef()); }; 
			if ( $bbb_len < 3 ) {
				my $addN = 'N' x (3-$bbb_len); 
				$bbb .= $addN; 
			}
			my ($curr_bbb2aa) = &bbb2aa( $bbb, $tblN ); 
			$curr_bbb2aa eq $aa[$j] or do { &tsmsg("[Err] Inconsistency between AA and CDS at seq_pos=$j [aa_pos=$aa_pos] [$aa[$j] vs $curr_bbb2aa vs $bbb]\n"); return(undef()); }; 
			$aa2cds_cseq .= "$bbb"; 
		}
	}
	my $clen = length($cseq); 
	if ( ($aa_pos+1)*3 < $clen ) {
		my $res = $clen % 3; 
		my $t_seq = $cseq; 
		$res > 0 and $t_seq .= ( 'N' x (3-$res) ); 
		my $t_pos = length($t_seq)/3-1; 
		for (my $ip = $t_pos; $ip > $aa_pos; $ip--) {
			my $bbb = substr( $t_seq, $ip*3, 3 ); 
			my ($curr_bbb2aa) = &bbb2aa( $bbb, $tblN ); 
			$curr_bbb2aa eq '*' or $curr_bbb2aa eq 'X' or do { &tsmsg("[Err] CDS seq is longer than prot seq: $cseq\n"); return(undef()); }; 
		}
	}
	return($aa2cds_cseq); 
}# aa2cds_1seq() 

=head1 get_4d_codon( $translation_table_number )

Return     : ( \%fourD_codon )

$fourD_codon{'ACA'} = [ 'T', 3 ]; # [ AA_symbol, frame ]
$fourD_codon{'TCA'} = [ 'S', 3 ]; 

=cut
sub get_4d_codon {
	my ($tblN) = @_; 
	my %back; 
	for my $b1 (qw/A T G C/) {
		for my $b2 (qw/A T G C/) {
			for my $sN (0 .. 2) {
				my $aa; 
				my $is_good = 1; 
				for my $b3 (qw/A T G C/) {
					my @tb = ($b1, $b2); 
					splice(@tb, $sN, 0, $b3); 
					my $bbb = join('', @tb); 
					my ($aa_t) = &bbb2aa($bbb); 
					$aa //= $aa_t; 
					$aa eq $aa_t or do { $is_good = 0; last; }; 
				}
				$is_good == 1 or next; 
				for my $b3 (qw/A T G C/) {
					my @tb = ($b1, $b2); 
					splice(@tb, $sN, 0, $b3); 
					my $bbb = join('', @tb); 
					$back{$bbb} = [ $aa, $sN+1 ]; 
				}
			}
		}
	}
	return(\%back); 
}# get_4d_codon() 

1; 
