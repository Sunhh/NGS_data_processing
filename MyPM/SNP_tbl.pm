package SNP_tbl; 
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
# 2014-12-23 Time flies fast! It's time to step on higher perl skills. I am coming. 

##########################################################
# Methods included: 
#  1. new()         : method to construct object. 
#  2. _initialize() : Initialization. 
#  3. newSubObj()   : Use 'RowN' and 'ColN' to get a new subset of an old object. 
#  4. rm_lmiss()    : Similar to rm_Nmiss, which removes SNP sites by missing rate per line (missing_indv%). 
#  5. rm_noVar()    : Remove SNP sites without accepted variations. (Only A/T/G/C alleles variance accepted). 
#  6. get_tv()      : Get Transversion SNP sites obj when 'type' is set as 'tv' (default). (type can be 'ts/tv/un')
#  7. readTbl()     : Shoud be used at first. Generate 'data_arr' Read in the whole SNP table. Maybe edit to read a subset in the future. 
#  8. is_tstv()     : Generate 'is_tstv' indicating if this site is 'ts'-transition, 'tv'-transversion, or 'un'-unknown(maybe more than two alleles)
#  9. SingleChar()  : Return a Single character for given DNA_base(s). &SingleChar( 'onlyATGC'=>0, 'maxAlleleN'=>0 ) 
#  10. SingleCharData() : Reformat the whole SNP table ('data_arr') to SingleChar format. &SingleCharData( 'onlyATGC'=>0, 'maxAlleleN'=>0 )
#  11. writeTbl()   : Output SNP table within current 'data_arr'. 'fmt' can be SNP/illumina. 
#  12. max2Allele() : Edit 'data_arr' to allow only Top2 counted allele per SNP site (line). Use 'force_single'=>1 to force reformat singleChar. 
#  13. tbl2illumina() : Output a table in illumina format. 'wrapTitle' could control if we wrap title with '"'; 
#  14. tbl2nex()    : Output a table in nexus format, which is used for phylip software. 
#  15. tbl2meg()    : Output a table in mega format, which is used by Mega. 
#  16. tbl2fas()    : Output a table in fasta foramt. 
#  17. tbl2structure() : Output an input table for Structure software. 
#  18. tbl2ped()    : Output input files ('ofile'=>"oo.ped", 'omapfile'=>"oo.map") for plink software. 
#  19. tbl2seq()    : A mixed function of many tbl2xx() functions. 
#  20. guessChrNum(): Given ($chr_txt, \%chr2num), it guess the best number of $chr_txt; 
#  21. cnt_maf()    : Generate 'cnt_maf'([]) for each SNP site. 
#  22. cnt_genotype(): Generate qw/cnt_alleleCnt cnt_alleleSum cnt_alleleTypeN cnt_alleleTypeBase cnt_lmiss cnt_imiss/ from 'data_arr'; 
#  23. skipLine()   : Used by readTbl(). Skip lines according to $self->{'skip'} number ; 
#  24. geno2num()   : Change DNA_bases to numbers (@); Check %allele2num ; 
# Internal functions: 
#  25. rmClose_idx() : Given a small-to-large sorted numbers' array, and 'within_dist'=>Value, it returns a indics array of numbers not close to any other. 
#  26. permutations(): Given (\@list_of_ele, $n_in_class), return all permutations by array. Return ([@perm1_of_ele], [@perm2], ...)
#  27. combinations(): Given (\@list_of_ele, $n_in_class), return all combinations by array. Return ([@comb1_of_ele], [@comb2_of_ele], ...)
##########################################################


##########################################################
#########  Basic settings. 
##########################################################

# IUPAC settings. 
my %IUPAC_b2d; # bases to degenerate; 
my %IUPAC_d2b; # degenerate to bases; 
{
my @ta = qw(
AA   CC   GG   TT   UU
WAT  SCG  MAC  KGT  RAG   YCT
BCGT DAGT HACT VACG NACGT
); 
	for my $tb (@ta) {
		my @bb = split(//, $tb); 
		my $dbase = shift(@bb); 
		$IUPAC_b2d{$dbase} = $dbase; 
		for (&permutations(\@bb, scalar(@bb))) {
			$IUPAC_b2d{join('', @$_)} = $dbase; 
		}
		for (&combinations(\@bb, scalar(@bb))) {
			$IUPAC_d2b{$dbase} = [sort @$_]; 
		}
	}
}# End of IUPAC_xxx; 

# For is_tstv test. 
my %tvtsV = qw(
A  1 G  1
T -1 C -1
); 

# For tbl2structure () test. 
my %allele2num = qw(
	A  1
	C  2
	G  3
	T  4
	N -9
); 

##########################################################
#########  Methods. 
##########################################################

=head1 new()

Function     : Create an SNP_tbl object. 
Input        : 
 'filename'  : input_SNP.tbl filename 
 'FH'        : File handle of filename
 'skip'      : [0] Number of head lines to skip (ignore) when reading in filename
 'header'    : [1] Tell if the input table has header line. 
 'chrColN'   : [0] Column No. of chromosome ID. 
 'posColN'   : [1] Column No. of position . 
 'refColN'   : [2] Column No. of reference base informaton. This is not used yet. 
 'chrColID'  : ['chrCol'] The title of column for chromosome ID. 
 'posColID'  : ['posCol'] The title of column for position information. 
 'skipColN'  : [] Columns that are not genotypes and will not exist in 'data_arr' array. 
                The 'chrCol' and 'posCol' will be skipped automatically. 
=cut
sub new {
	my $class = shift;     # Get the request class name. Use "class->new()" to call this function! 
	
	my $self = {}; 
	bless $self, $class; 

	$self->_initialize(@_); 
	
	return $self; 
}#sub new

sub _initialize {
	my $self = shift; 
	my %parm = @_; 
	$self->{'filename'} = $parm{'filename'}; 
	$self->{'FH'}   = $parm{'FH'}; 
	$self->{'skip'} = $parm{'skip'} // 0; 
	$self->{'header'} = $parm{'header'} // 1; 
	$self->{'chrColN'} = $parm{'chrColN'} // 0; 
	$self->{'posColN'} = $parm{'posColN'} // 1; 
	$self->{'refColN'} = $parm{'refColN'} // 2; 
	$self->{'chrColID'} = 'chrCol'; 
	$self->{'posColID'} = 'posCol'; 
	$self->{'skipColN'} = $parm{'skipColN'}; 
	my %th = %$self; 
	for my $tn (@th{qw/chrColN posColN/}) {
		$self->{'skipColN'}{$tn} = $self->{'skipColN'}{$tn} // 1; 
	}
}#End sub _initialize() 

=head2 newSubObj( 'RowN'=>[-1], 'ColN'=>[-1] )

Function     : Creat a new object which is a subset of old_current object, according to 'RowN' and 'ColN' information. 
 'RowN'   : [-1], should be a \@arr_indice of row numbers. 0-based. 
 'ColN'   : [-1], should be a \@arr_indice of column numbers. 0-based. 
Return       : A new object with 'data_arr' is a subset of old object. 

#  3. newSubObj()   : Use 'RowN' and 'ColN' to get a new subset of an old object. 
=cut
sub newSubObj {
	my $self = shift; 
	my %parm = @_; 
	$parm{'RowN'} = $parm{'RowN'} // [-1]; 
	$parm{'ColN'} = $parm{'ColN'} // [-1]; 
	my $new_obj = {}; 
	bless $new_obj, ref($self); # I don't know other better metod. 
	
	# Setting data_arr
	defined $parm{'RowN'}[0] and $parm{'RowN'}[0] == -1 and $parm{'RowN'} = [ 0 .. $#{$self->{'data_arr'}}    ]; 
	defined $parm{'ColN'}[0] and $parm{'ColN'}[0] == -1 and $parm{'ColN'} = [ 0 .. $#{$self->{'data_arr'}[0]} ]; 
	$new_obj->{'data_arr'} = []; 
	for my $idx1 ( @{$parm{'RowN'}} ) {
		push(@{$new_obj->{'data_arr'}}, [ @{$self->{'data_arr'}[$idx1]}[ @{$parm{'ColN'}} ] ]); 
	}
	
	# Setting title / header / chrCol / posCol 
	# By Key 
	for my $kk (qw/header filename FH skip header chrColN posColN refColN chrColID posColID skipColN chrIDrel/) {
		defined $self->{$kk} and $new_obj->{$kk} = $self->{$kk}; 
	}
	# By Row 
	for my $rowKey (qw/chrColV posColV/) {
		if (defined $self->{$rowKey}) {
			$new_obj->{$rowKey} = [ @{$self->{$rowKey}}[ @{$parm{'RowN'}} ] ]; 
		}
	}
	# By Col
	for my $colKey (qw/title/) {
		if (defined $self->{$colKey}) {
			$new_obj->{$colKey} = [ @{$self->{$colKey}}[ @{$parm{'ColN'}} ] ]; 
		}
	}
	
	return $new_obj; 
}# sub newSubObj () 

# Remove SNP sites by genotyped_missing_indv% (missing per line). 
sub rm_lmiss {
	my $self = shift; 
	my %parm = @_; 
	$parm{'maxAllow'} = $parm{'maxAllow'} // 0.10; # Default 10% 
	$self->cnt_genotype(); 
	$parm{'ttl_indvN'} = $parm{'ttl_indvN'} // $#{$self->{'data_arr'}[0]}+1 ; 
	my $min_indvN = $parm{'ttl_indvN'} * $parm{'maxAllow'} ; 
	my $new_obj; 
	my @goodRows; 
	for (my $i=0; $i<@{$self->{'cnt_lmiss'}}; $i++) {
		$self->{'cnt_lmiss'}[$i] <= $min_indvN and push(@goodRows, $i); 
	}
	$new_obj = $self->newSubObj('RowN'=>\@goodRows); 
	return $new_obj; 
}# sub rm_lmiss

sub rm_noVar {
	my $self = shift; 
	$self->cnt_genotype(); 
	my @goodRows; 
	for (my $i=0; $i<@{$self->{'cnt_alleleTypeN'}}; $i++) {
		$self->{'cnt_alleleTypeN'}[$i] > 1 and push(@goodRows, $i); 
	}
	my $new_obj = $self->newSubObj('RowN'=>\@goodRows); 
	return $new_obj; 
}# sub rm_noVar () 

=head2 get_tv( 'type'=>'tv|ts|un' ) 

Required      : None. 
 'type' => 'tv' -- retrieve tranversion mutation. 
           'ts' -- retrieve transition mutation. 
           'un' -- retrieve non-tv and non-ts sites, which means more than two alleles in this site
Function      : Return a new object with sites in required type. 
                This function will use is_tstv() method to determine the type of each site. 
Return        : SNP_tbl object. 
=cut
sub get_tv {
	my $self = shift; 
	$self->is_tstv(); 
	my %parm = @_; 
	$parm{'type'} = $parm{'type'} // 'tv'; # Could be tv/ts/un
	my @goodRows; 
	for (my $i=0; $i<@{$self->{'is_tstv'}}; $i++) {
		$self->{'is_tstv'}[$i] eq $parm{'type'} and push(@goodRows, $i); 
	}
	my $new_obj = $self->newSubObj('RowN'=>\@goodRows); 
	return $new_obj; 
}# sub get_tv ()  

sub readTbl {
	my $self = shift; 
	for my $tk (qw/FH title data_arr/) {
		undef $self->{$tk}; 
	}
	my %parm = @_; 
	$self->{'filename'} = $parm{'filename'} // $self->{'filename'}; 
	$self->{'FH'} = $parm{'FH'} // &openFH( $self->{'filename'}, '<' ); 
	$self->skipLine(); 
	
	# Read in header. 
	if ( $self->{'header'} ) {
		my $tl = readline($self->{'FH'}); 
		chomp($tl); 
		my @ta = split(/\t/, $tl); 
		for (my $i=0; $i<@ta; $i++) {
			defined $self->{'skipColN'}{$i} or push(@{$self->{'title'}}, $ta[$i]); 
		}
		for my $tpref (qw/chrCol posCol/) {
			defined $self->{"${tpref}N"} and $self->{"${tpref}ID"} = $ta[ $self->{"${tpref}N"} ]; 
		}
	}# if ( header )
	
	# Clean data_idx
	for my $tpref (qw/chrCol posCol/) {
		$self->{"${tpref}V"} = []; 
	}
	# Read in data_arr 
	while (my $tl = readline($self->{'FH'})) {
		chomp($tl); 
		my @ta = split(/\t/, $tl); 
		for my $tpref (qw/chrCol posCol/) {
			defined $self->{"${tpref}N"} and push(@{$self->{"${tpref}V"}}, $ta[ $self->{"${tpref}N"} ]); 
		}
		my @tb; 
		for (my $i=0; $i<@ta; $i++) {
			defined $self->{'skipColN'}{$i} or push(@tb, $ta[$i]); 
		}
		push(@{$self->{'data_arr'}}, [@tb]); 
	}
	close ($self->{'FH'}); 
	$self->{'FH'} = undef(); 
	return 0; 
}# sub readTbl()

# construct an array in 'is_tstv' : 'ts/tv/un' for each SNP site. 
#   'un' means there are more than two alleles in this site. 
sub is_tstv {
	my $self = shift; 
	for my $tk (qw/is_tstv/) {
		undef $self->{$tk}; 
	}
	( defined $self->{'cnt_alleleCnt'} and scalar(@{$self->{'cnt_alleleCnt'}}) == scalar(@{$self->{'data_arr'}}) ) or $self->cnt_genotype(); 
	$self->{'is_tstv'} = []; 
	for (my $i=0; $i<@{$self->{'cnt_alleleCnt'}}; $i++) {
		$_ = $self->{'cnt_alleleCnt'}[$i]; 
		my @ta = keys %$_; 
		if (scalar(@ta) == 2) {
			push(@{$self->{'is_tstv'}}, ( $tvtsV{$ta[0]} * $tvtsV{$ta[1]} == -1 ) ? 'tv' : 'ts'); 
		}else{
			push(@{$self->{'is_tstv'}}, 'un'); 
		}
	}
	return 0; 
}# sub is_tstv()

sub SingleChar {
	my $self = shift; 
	my $char = shift; 
	$char = uc($char); 
	my %parm = @_; 
	$parm{'onlyATGC'} = $parm{'onlyATGC'} // 0; 
	$parm{'maxAlleleN'} = $parm{'maxAlleleN'} // 0; 
	if ( $parm{'onlyATGC'} ) {
		return ($char =~ m/^[ATGC]$/i) ? $char : 'N' ; 
	} else {
		while ( $char =~ s!((.).*)\2!$1!g ) { 1; } 
		if (defined $IUPAC_b2d{$char}) {
			if ($parm{'maxAlleleN'} > 0) {
				return ( @{$IUPAC_d2b{$IUPAC_b2d{$char}}} <= $parm{'maxAlleleN'} ) ? $IUPAC_b2d{$char} : 'N' ; 
			} else {
				return $IUPAC_b2d{$char}; 
			}
		} else {
			return 'N'; 
		}
	}
	return 1; 
}# sub SingleChar () 

=head2 SingleCharData( 'onlyATGC'=>0, 'maxAlleleN'=>0 )

Required : Need 'data_arr' available in current object. 
Function : Change genotypes in 'data_arr' to single character.
 'onlyATGC'    : If == 1, only accept 'A|T|G|C' characters in raw input genotypes. 
 'maxAlleleN'  : If > 0, any genotypes with more than 'maxAlleleN' 'A|T|G|C' probable bases will be changed to 'N'
Return   : None. $this_object->{'data_arr'} has been changed. 
=cut
sub SingleCharData {
	my $self = shift; 
	for my $tk (qw/is_single/) {
		undef $self->{$tk}; 
	}
	my %parm = @_; 
	my %pass_parm = %parm; 
	$parm{'onlyATGC'} = $parm{'onlyATGC'} // 0; 
	# Not to pass some parameters. 
	for my $toRM (qw//) {
		exists $pass_parm{$toRM} and delete $pass_parm{$toRM}; 
	}
	for my $t1 (@{$self->{'data_arr'}}) {
		for my $t2 (@$t1) {
			$t2 = uc($t2); # Repeated but I want to use it. 
			$t2 = $self->SingleChar($t2, %pass_parm); 
		}
	}
	$self->{'is_single'} = 1; 
	return 0; 
}# SingleCharData

sub writeTbl {
	my $self = shift; 
	my %parm = @_; 
	my $ofh = ( defined $parm{'ofile'} ) ? &openFH($parm{'ofile'}, '>') : \*STDOUT ; 
	$parm{'fmt'} = $parm{'fmt'} // 'SNP'; 
	$parm{'wrapTitle'} = $parm{'wrapTitle'} // 0; 
	
	# Output header
	if ( $self->{header} ) {
		my @new_hh = @{$self->{'title'}}; 
		if ( $parm{'wrapTitle'} ) {
			for my $thh (@new_hh) {
				$thh = "\"$thh\""; 
			}
		}
		if ($parm{'fmt'} =~ m/^SNP$/i) {
			my @to_ar; 
			for my $tpref ( qw/chrCol posCol/ ) {
				defined $self->{"${tpref}N"} and push(@to_ar, $self->{"${tpref}ID"}); 
			}
			print $ofh join("\t", @to_ar, @new_hh)."\n"; 
		} elsif ( $parm{'fmt'} =~ m/^illu(?:mina)$/i ) {
			my @to_ar = (qw/name chr pos/); 
			print $ofh join("\t", @to_ar, @new_hh)."\n"; 
		} else {
			&stopErr("[Err] Unknown fmt value [$parm{'fmt'}] for subroutine " . (caller(0))[3] . "\n"); 
		}
	}# if (header) 
	# Output data_arr
	if ( $parm{'fmt'} =~ m/^SNP$/i ) {
		for (my $i=0; $i<@{$self->{'data_arr'}}; $i++) {
			my @to_ar; 
			for my $tpref ( qw/chrCol posCol/ ) {
				defined $self->{"${tpref}N"} and push(@to_ar, $self->{"${tpref}V"}[$i]); 
			}
			print $ofh join("\t", @to_ar, @{$self->{'data_arr'}[$i]})."\n"; 
		}# for (data_arr)
	} elsif ( $parm{'fmt'} =~ m/^illu(?:mina)$/i ) {
		$self->{header} or &stopErr("[Err] Illumina format needs [header] line to assign [chr] and [pos] as well as [name] columns.\n"); 
		# ( defined $self->{'is_single'} and $self->{'is_single'} == 1 ) or $self->SingleCharData('maxAlleleN'=>2); 
		$self->SingleCharData('maxAlleleN'=>2); 
		$parm{'max2Allele'} = $parm{'max2Allele'} // 0; 
		$parm{'max2Allele'} and $self->max2Allele(); 
		defined $self->{'chrIDrel'} or $self->{'chrIDrel'} = {}; 
		for (my $i=0; $i<@{$self->{'data_arr'}}; $i++) {
			my @to_ar; 
			my $chrNum ; 
			($chrNum, $self->{'chrIDrel'}) = &guessChrNum( $self->{'chrColV'}[$i] ,$self->{'chrIDrel'} ); 
			my $posNum = $self->{'posColV'}[$i]; 
			@to_ar = ("s${chrNum}_${posNum}", $chrNum, $posNum); 
			print $ofh join("\t", @to_ar); 
			for my $tb (@{$self->{'data_arr'}[$i]}) {
				if ( $tb eq 'N') { 
					print {$ofh} "\t"."00"; 
				} elsif ( defined $IUPAC_d2b{$tb} ) {
					my $tc = join('', @{$IUPAC_d2b{$tb}}); 
					print {$ofh} "\t".substr("$tc$tc",0,2); 
				} else {
					print {$ofh} "\t"."00"; 
				}
			}
			print $ofh "\n"; 
		}# for (data_arr)
	} else {
		&stopErr("[Err] Unknown fmt value [$parm{'fmt'}] for subroutine " . (caller(0))[3] . "\n"); 
	}
	return 0; 
}# sub writeTbl 

### Keep only two top counted alleles in genotypes. 
### Genotypes with other alleles are changed to "N"
sub max2Allele {
	my $self = shift; 
	my %parm = @_; 
	$parm{'force_single'} = $parm{'force_single'} // 0; 
	if ( $parm{'force_single'} ) {
		$self->SingleCharData('maxAlleleN'=>2); 
	} elsif ( defined $self->{'is_single'} and $self->{'is_single'} == 1 ) {
		; 
	} else {
		$self->SingleCharData('maxAlleleN'=>2); 
	}
	$self->cnt_genotype(); 
	for (my $i=0; $i<@{$self->{'data_arr'}}; $i++) {
		# Check each [$i] SNP site. 
		# Get top 2 alleles 
		my %good; 
		for (my $j=0; $j<@{$self->{'cnt_alleleTypeBase'}[$i]} and $j <= 1; $j++) {
			$good{ $self->{'cnt_alleleTypeBase'}[$i][$j] } = 1; 
		}
		# Convert genotype to good, final alleles in SingleChar format. 
		# Filter genotypes including non-good alleles by changing them to 'N'. 
		for my $tb (@{$self->{'data_arr'}[$i]}) {
			my $is_toN = 0; 
			for my $tc (@{$IUPAC_d2b{$tb}}) {
				defined $good{$tc} or do { $is_toN = 1; last; }; 
			}
			$is_toN == 1 and $tb = 'N'; 
		}
	}
	return 0; 
}# sub max2Allele 

sub tbl2illumina {
	my $self = shift; 
	return $self->writeTbl('fmt'=>'illumina', 'wrapTitle'=>1, @_); 
}#sub tbl2illumina

sub tbl2nex {
	my $self = shift; 
	return $self->tbl2seq('fmt'=>'nex', @_); 
}#sub tbl2nex ()

sub tbl2meg {
	my $self = shift; 
	return $self->tbl2seq('fmt'=>'meg', @_); 
}#sub tbl2meg

sub tbl2fas {
	my $self = shift; 
	return $self->tbl2seq('fmt'=>'fasta', @_); 
}#sub tbl2fas

sub tbl2structure {
	my $self = shift; 
	return $self->tbl2seq('fmt'=>'structure', @_); 
}#sub tbl2structure

sub tbl2ped {
	my $self = shift; 
	return $self->tbl2seq('fmt'=>'ped', @_); 
}

sub tbl2seq {
	my $self = shift; 
	my %parm = @_; 
	&tsmsg("[Msg] Preparing data.\n"); 
	defined $self->{'data_arr'} or $self->readTbl(); 
	$parm{'onlyATGC'} = $parm{'onlyATGC'} // 0; 
	&tsmsg("[Msg] Formatting data.\n"); 
	( defined $self->{'is_single'} and $self->{'is_single'} == 1 ) or $self->SingleCharData('onlyATGC'=>$parm{'onlyATGC'}); 
	my $ofh = ( defined $parm{'ofile'} ) ? &openFH($parm{'ofile'}, '>') : \*STDOUT ; 
	$parm{'rm_noVar'} = $parm{'rm_noVar'} // 0; 
	$parm{'seqWid'} = $parm{'seqWid'} // 100; 
	$parm{'more_infor'} = $parm{'more_infor'} // 0; 
	
	$parm{'fmt'} = $parm{'fmt'} // 'fasta'; 
	defined $self->{'chrIDrel'} and $self->{'chrIDrel'} = {}; 
	if ( $parm{'fmt'} =~ m/^ped$/i ) {
		&tsmsg("[Msg] Limiting only two alleles per site.\n"); 
		$self->max2Allele('force_single'=>1); # 
		&tsmsg("[Msg] Removing no variation sites.\n"); 
		$self = $self->rm_noVar(); 
		$parm{'omapfile'} = $parm{'omapfile'} // 'oo.map' ; 
	}
	$parm{'rm_noVar'} and $self->rm_noVar(); 
	
	&tsmsg("[Msg] Generating sequences.\n"); 
	my @hh = @{$self->{'title'}}; # Sequence's IDs 
	my @seqs;                     # Sequences string. 
	my ($n_tax, $n_char) = (0,0); 
	$n_tax  = scalar(@hh); 
	$n_char = scalar(@{$self->{'data_arr'}}); 
	for (my $i=0; $i<@{$self->{'data_arr'}}; $i++) {
		for (my $k = 0; $k<@{$self->{'data_arr'}[$i]}; $k++) {
			if ( $parm{'onlyATGC'} ) {
				$seqs[$k] .= ($self->{'data_arr'}[$i][$k] =~ m/^([ATGC])$/i ) ? $1 : 'N' ; 
			} else {
				$seqs[$k] .= $self->{'data_arr'}[$i][$k]; 
			}
		}
	}
	
	### Different format to output. 
	&tsmsg("[Msg] Output files.\n"); 
	if ($parm{'fmt'} =~ m/^fa(?:s(?:ta))?$/i) {
		for (my $i=0; $i<@hh; $i++) {
			my $ts = $seqs[$i]; 
			$ts =~ s/(.{$parm{'seqWid'}})/$1\n/g; chomp($ts); 
			my $addHead = ( defined $parm{'more_infor'} ) ? " [Len:$n_char]" : '' ; 
			print {$ofh} ">$hh[$i]$addHead\n$ts\n"; 
		}
	} elsif ($parm{'fmt'} =~ m/^meg(?:a)?$/i) {
		print {$ofh} <<MEGA; 
#mega
!Title tbl2meg; 

MEGA
		for (my $i=0; $i<@hh; $i++) {
			my $ts = $seqs[$i]; 
			$ts =~ s/(.{$parm{'seqWid'}})/$1\n/g; chomp($ts); 
			print {$ofh} "#$hh[$i]\n$ts\n"; 
		}
	} elsif ($parm{'fmt'} =~ m/^nex(?:us)?$/i) {
		print {$ofh} <<HH;
#NEXUS

BEGIN DATA;
 DIMENSIONS NTAX=$n_tax NCHAR=$n_char;
 FORMAT DATATYPE=DNA INTERLEAVE MISSING=N;
HH
		for ( my $i=0; $i<@hh; $i++ ) {
			print {$ofh} "[Name: $hh[$i] Len:$n_char  Check:    0]\n"; 
		}
		print {$ofh} "\nMATRIX\n"; 
		my $len_perLine = 100; 
		for (my $i=0; $i<$n_char; $i+=$len_perLine) {
			my $j=$i+$len_perLine-1; 
			$j >= $n_char and $j=$n_char-1; 
			for (my $k = 0; $k<@hh; $k++) {
				print {$ofh} "$hh[$k]    " . substr($seqs[$k], $i, $j-$i+1) . "\n"; 
			}
			print {$ofh} "\n"; 
		}
		print {$ofh} ";\nEND;\n"; 
	} elsif ( $parm{'fmt'} =~ m/^struc(?:ture)?$/i ) { 
		for (my $i=0; $i<@hh; $i++) {
			my ($s1, $s2) = ($hh[$i], $hh[$i]); 
			for my $tb (split(//, $seqs[$i])) {
				my @nn = &geno2num($tb); 
				$s1 .= "\t$nn[0]"; 
				$s2 .= "\t$nn[1]"; 
			}
			print {$ofh} "$s1\n"; 
			print {$ofh} "$s2\n"; 
		}
	} elsif ( $parm{'fmt'} =~ m/^ped$/i ) { 
		# http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml
		#  It says all markers should be biallelic. 
		
		
		# Write .map file. 
		my $oMfh; 
		$oMfh = &openFH( $parm{'omapfile'}, '>' ); 
		for (my $i=0; $i<@{$self->{chrColV}}; $i++) {
			my $chrNum; 
			($chrNum, $self->{'chrIDrel'}) = &guessChrNum( $self->{'chrColV'}[$i] ,$self->{'chrIDrel'} ); 
			my $posNum = $self->{'posColV'}[$i]; 
			print {$oMfh} join("\t", $chrNum, "s${chrNum}_${posNum}", 0, $posNum)."\n"; 
		}
		close ($oMfh); 
		&tsmsg("[Msg] omapfile=$parm{'omapfile'} written.\n"); 
		
		for (my $i=0; $i<@hh; $i++) {
			my $pedID = $hh[$i]; 
			my $idvID = $hh[$i]; 
			my $fatID = 0; 
			my $monID = 0; 
			my $sexID = 1; 
			my $case  = 0; 
			print {$ofh} join("\t", $pedID, $idvID, $fatID, $monID, $sexID, $case); 
			my @geno_dd; 
			for my $tb (split(//, $seqs[$i])) {
				my @nn = &geno2num($tb); 
				for my $tn (@nn) { $tn == -9 and $tn = 0; }
				push(@geno_dd, "$nn[0] $nn[1]"); 
			}
			print {$ofh} join("\t", '', @geno_dd)."\n"; 
		}
	} else {
		&stopErr("[Err] Unknown fmt value [$parm{'fmt'}] for subroutine " . (caller(0))[3] . "\n"); 
	}
	
	return 0; 
}#sub tbl2seq ()

sub guessChrNum {
	my ($raw_id, $raw_hash) = @_; 
	my ($back_id, %back_hash); 
	$raw_hash = $raw_hash // \%back_hash; 
	%back_hash = %{$raw_hash}; 
	if ( defined $back_hash{c2n}{$raw_id} ) {
		$back_id = $back_hash{c2n}{$raw_id}; 
	} else {
		if ( $raw_id =~ m/^c(?:hr(?:om(?:osome)?)?)?(\d+)$/i ) {
			$back_id = $1; 
		} elsif ( $raw_id =~ m/^(\d+)$/i ) {
			$back_id = $1; 
		} elsif ( $raw_id =~ m/^(WM97_)?c(?:hr(?:om(?:osome)?)?)?(\d+)$/ ) {
			$back_id = $1; 
		} else {
			$back_id = 1; 
			while (1) {
				defined $back_hash{n2c}{$back_id} or last; 
				$back_id ++; 
			}
		}
		if ( defined $raw_hash->{n2c}{$back_id} ) {
			&tsmsg("[Err] Already defined relationship between Num:[$back_id] to ID:[$raw_id]!\n"); 
			$back_id = 1; 
			while (1) {
				defined $back_hash{n2c}{$back_id} or last; 
				$back_id ++; 
			}
			&tsmsg("[Err] Now re-assign a new Num:[$back_id] for ID:[$raw_id]\n"); 
		}
		$raw_hash->{c2n}{$raw_id} = $back_id; 
		$raw_hash->{n2c}{$back_id} = $raw_id; 
	}#End if else 
	return($back_id, $raw_hash); 
}# sub guessChrNum ()


# Calculate MAF for each SNP site. 
sub cnt_maf {
	my $self = shift; 
	( defined $self->{'cnt_alleleCnt'} and scalar(@{$self->{'cnt_alleleCnt'}}) == scalar(@{$self->{'data_arr'}}) ) or $self->cnt_genotype(); 
	$self->{'cnt_maf'} = []; 
	for (my $i=0; $i<@{$self->{'cnt_alleleCnt'}}; $i++) {
		$_ = $self->{'cnt_alleleCnt'}[$i]; 
		if ($self->{'cnt_alleleTypeN'}[$i] <= 1 or $self->{'cnt_alleleSum'}[$i] <= 0) {
			push(@{$self->{'cnt_maf'}}, -1); 
		} else {
			my @ta = sort { $_->{$b} <=> $_->{$a} } keys %{$_}; 
			push(@{$self->{'cnt_maf'}}, $_->{$ta[1]} / $self->{'cnt_alleleSum'}[$i] * 100); # {'cnt_maf'} = [ MAF_ratio1, MAF_ratio2, ... ]; 
		}
	}
	return 0; 
}# sub cnt_maf ()

=head2 cnt_genotype()

Required     : object->{'data_arr'} filled. 
Function     : Count genotype information for each site (line). 
Return       : obj->{'tag'}, most of which are for related sites, when tag == 
 {'cnt_alleleCnt'}      : [\%cnt], {[ATGC]=>base_count}
 {'cnt_alleleSum'}      : [$sum], Sum of allele numbers (AA gives 2 count)
 {'cnt_alleleTypeN'}    : [scalar(keys %cnt)] , Number of allele types in this site. 
 {'cnt_alleleTypeBase'} : [ sorted keys %cnt ], Allele bases sorted by count number large to small 
 {'cnt_lmiss'}          : [$lmiss], Number of individuals missing ("N") genotype in the current site. This value doesn't include number of Indels! 
 {'cnt_imiss'}          : \@imiss, Number of sites missing with in individual. 
 {'cnt_hete'}           : [$hete_indv_number] Only genotypes with 2 alleles ('A|T|G|C') counted. 
 {'cnt_homo'}           : [$homo_indv_number] Only genotypes with 'A|T|G|C' counted. 
=cut
# Only count genotypes in homo/biHete SNP. 
sub cnt_genotype {
	my $self = shift; 
	defined $self->{'data_arr'} or return 0; 
	for my $tkey (qw/cnt_alleleCnt cnt_alleleSum cnt_alleleTypeN cnt_alleleTypeBase cnt_lmiss cnt_imiss cnt_hete cnt_homo/) {
		undef $self->{$tkey}; 
	}
	my @imiss; # Same to definition in vcftools. the missingness on a per-individual basis
	for my $tp (@{$self->{'data_arr'}}) {
		my %cnt; 
		my $sum = 0; 
		my $lmiss = 0; 
		my $i_idx = -1; 
		my $cnt_hete = 0; 
		my $cnt_homo = 0; 
		for my $tb (@$tp) {
			$i_idx ++; 
			$tb = uc($tb); 
			$tb eq 'N' and do { $lmiss++; $imiss[$i_idx]++; next; }; 
			( $tb =~ m/\*/ or $tb =~ m/\+/) and next; 
			length($tb) <= 2 or next; 
			defined $IUPAC_b2d{$tb} or do { $lmiss++; $imiss[$i_idx]++; next; }; 
			my @tc = @{$IUPAC_d2b{ $IUPAC_b2d{$tb} }}; 
			if ( scalar(@tc) == 1 ) {
				$cnt{$tc[0]} += 2; 
				$sum += 2; 
				$cnt_homo ++; 
			} elsif ( scalar(@tc) == 2 ) {
				$cnt{$tc[0]} ++; 
				$cnt{$tc[1]} ++; 
				$sum += 2; 
				$cnt_hete ++; 
			} else {
				$lmiss++; $imiss[$i_idx]++; next; 
			}
		}
		push(@{$self->{'cnt_alleleCnt'}}, \%cnt);   # {[ATGC]=>value} recording counts of allele, 
		push(@{$self->{'cnt_alleleSum'}}, $sum); # Sum of allele numbers (AA gives 2)
		push(@{$self->{'cnt_alleleTypeN'}}, scalar(keys %cnt)); # Number of allele types in this site. 
		push(@{$self->{'cnt_alleleTypeBase'}}, [sort { $cnt{$b} <=> $cnt{$a} || $a cmp $b } keys %cnt]); # Allele bases sorted by count number. 
#		push(@{$self->{'cnt_alleleTypeBase'}}, [sort { $cnt{$b} <=> $cnt{$a} } keys %cnt]); # Allele bases sorted by count number. 
		push(@{$self->{'cnt_hete'}}, $cnt_hete); 
		push(@{$self->{'cnt_homo'}}, $cnt_homo); 
		push(@{$self->{'cnt_lmiss'}}, $lmiss); # Number of individuals missing ("N") genotype in the current site. This value doesn't include number of Indels. 
	}
	$self->{'cnt_imiss'} = \@imiss; 
	return 0; 
}# sub cnt_genotype () 


sub skipLine {
	my $self = shift; 
	int($self->{skip}) > 0 or return 0; 
	for (1 .. int($self->{skip})) { readdir($self->{'FH'}); } 
	return 0; 
}# sub skipLine

#######################################################
# Independent methods. 
#######################################################

#  Check if there is any variant in given \@array, and return 0-noVariant 1-HasVariant; 
#  Input ( 'arr'  =>["AA", "a", "N", "*", "A+C", ...], 
#          'mode' =>'single/skipN'/'skipIndel'/'skipHete'/
#        )
sub chk_VarInArray {
	my $self = shift; 
	my %parm = @_; 
	defined $parm{'arr'} or &stopErr("[Err] No array information asigned.\n"); 
	my @arr = @{$parm{'arr'}}; 
	$parm{'mode'} = $parm{'mode'} // 'simple'; 
	$parm{'mode'} = lc($parm{'mode'}); 
	if ( $parm{'mode'} eq 'simple' or $parm{'mode'} eq 'skipn' ) {
		my $baseRef = ''; 
		my $is_var = 0; 
		for (@arr) {
			$_ =~ m/n/i and next; 
			$_ = uc(@_); 
			$baseRef eq '' and $baseRef = $_; 
			$baseRef ne $_ and do { $is_var = 1; last; }; 
		}
		return $is_var; 
	} elsif ( $parm{'mode'} eq 'skipindel' ) {
		my $baseRef = ''; 
		my $is_var = 0; 
		for (@arr) {
			$_ =~ m/n/i and next; 
			$_ =~ m/\*|\+/i and next; 
			$_ = uc(@_); 
			$baseRef eq '' and $baseRef = $_; 
			$baseRef ne $_ and do { $is_var = 1; last; }; 
		}
		return $is_var; 
	} elsif ( $parm{'mode'} eq 'skiphete' ) {
		my $baseRef = ''; 
		my $is_var = 0; 
		for (@arr) {
			$_ =~ m/n/i and next; 
			$_ =~ m/^([ATGC])\1*$/i or next; 
			$_ = uc(@_); 
			$baseRef eq '' and $baseRef = $_; 
			$baseRef ne $_ and do { $is_var = 1; last; }; 
		}
		return $is_var; 
	} elsif ( $parm{'mode'} eq 'onlyhomo' or $parm{'mode'} eq 'skipindelhete' or $parm{'mode'} eq 'skipheteindel' ) {
		my $baseRef = ''; 
		my $is_var = 0; 
		for (@arr) {
			$_ =~ m/n/i and next; 
			$_ =~ m/\*|\+/i and next; 
			$_ =~ m/^([ATGC])\1*$/i or next; 
			$_ = uc(@_); 
			$baseRef eq '' and $baseRef = $_; 
			$baseRef ne $_ and do { $is_var = 1; last; }; 
		}
		return $is_var; 
	} else {
		&stopErr("[Err] Bad mode type [$parm{'mode'}].\n"); 
	}
	return undef(); 
}# sub chk_VarInArray () 


#######################################################
# Internal sub-routine. 
#######################################################

# Function for tbl2structure and tbl2ped
sub geno2num {
	my $tb = shift; 
	my @bb = @{$IUPAC_d2b{$tb}}; 
	my $bC = scalar(@bb); 
	if ( $bC == 1 ) {
		return ( $allele2num{$bb[0]}, $allele2num{$bb[0]}); 
	} elsif ( $bC == 2 ) { 
		return ( $allele2num{$bb[0]}, $allele2num{$bb[1]} ); 
	} else {
		return ( $allele2num{'N'}, $allele2num{'N'} ); 
	}
}# sub geno2num ()


sub rmClose_idx {
	my $ar = shift; 
	my %parm = @_; 
	$parm{'within_dist'} = $parm{'within_dist'} // 5; 
	
	my @backIdx; 
	my @t2 = @$ar; 
	my %prev; 
	$prev{'idx'} = 0; 
	$prev{'pos'} = $t2[0]; 
	$prev{'is_good'} = 1; 
	for (my $i=1; $i<@t2; $i++) {
		my $curr_good = 1; 
		my $dist2prev = $t2[$i]-$prev{'pos'}+1; 
		if ( $dist2prev <= $parm{'within_dist'} ) {
			$curr_good = 0; 
		} else {
			$prev{'is_good'} == 1 and push(@backIdx, $prev{'idx'}); 
		}
		$prev{'is_good'} = $curr_good; 
		$prev{'idx'} = $i; 
		$prev{'pos'} = $t2[$i]; 
	}

	return \@backIdx; 
}# sub rmClose_idx () 

sub permutations {
	my ($list, $n) = @_; 
	$n = $n // scalar(@$list); 
	$n > @$list and return ($list); 
	$n <= 1 and return(map {[$_]} @$list); 
	my @perm; 
	for my $i (0 .. $#$list) {
		my @rest = @$list; 
		my $val = splice(@rest, $i, 1); 
		for ( &permutations(\@rest, $n-1) ) {
			push(@perm, [$val, @$_]); 
		}
	}
	return @perm; 
}#sub permutations() 

sub combinations {
	my ($list, $n) = @_; 
	$n = $n // scalar(@$list); 
	$n > @$list and return ($list); 
	$n <= 1 and return(map {[$_]} @$list); 
	my @comb; 
	for (my $i=0; $i+$n<=@$list; $i++) {
		my $val = $list->[$i]; 
		my @rest = @$list[$i+1 .. $#$list]; 
		for (&combinations(\@rest, $n-1)) {
			push(@comb, [$val, @$_]); 
		}
	}
	return @comb; 
}#sub combinations


1; # Terminate the package with the required 1; 
