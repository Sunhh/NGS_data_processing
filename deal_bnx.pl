#!/usr/bin/perl 
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"in_bnx:s", # Input.bnx file. 
	"rand_slct:i", # 5000; 
	"min_len:f",   # 100000; 
	"out:s", # 
); 

$opts{'outFH'} = \*STDOUT; 
defined $opts{'out'} and $opts{'outFH'} = &openFH($opts{'out'}, '>'); 

my $help_txt = <<HH; 

perl $0 -in_bnx in.bnx 

-help 

-out               [STDOUT] Output. 

-rand_slct         [5000] No replacement. 

-min_len           [100000] Remove molecule map shorter than ... 

HH
$opts{'help'} and &LogInforSunhh::usage($help_txt); 

if (defined $opts{'rand_slct'}) {
	&rand_slct(); 
} elsif (defined $opts{'min_len'}) {
	&rm_shrt(); 
}

####################################################################################################
#      Function sub-routines. 
####################################################################################################
sub rm_shrt {
	defined $opts{'in_bnx'} or &LogInforSunhh::usage($help_txt); 
	&tsmsg("[Msg] Loading bnx file [$opts{'in_bnx'}]\n"); 
	my $bnx_data = &_load_whole_bnx($opts{'in_bnx'}); 
	&_write_bnx_header( $opts{'outFH'}, $bnx_data->{'header'} ); 
	$bnx_data->{'records'} = &_filter_mol_length( $bnx_data->{'records'}, $opts{'min_len'} ); 
	&_write_bnx_records( $opts{'outFH'}, $bnx_data->{'records'} ); 
	return; 
}# rm_shrt() 

sub rand_slct {
	$opts{'rand_slct'} > 0 or $opts{'rand_slct'} = 5000; 
	defined $opts{'in_bnx'} or &LogInforSunhh::usage($help_txt); 
	&tsmsg("[Msg] Loading bnx file [$opts{'in_bnx'}]\n"); 
	my $bnx_data = &_load_whole_bnx($opts{'in_bnx'}) ; 
	&_write_bnx_header( $opts{'outFH'}, $bnx_data->{'header'} ); 
	# chomp( $bnx_data->{'header'} ); print {$opts{'outFH'}} $bnx_data->{'header'} . "\n"; 
	my $max_cnt = scalar( @{$bnx_data->{'records'}} ); 
	&tsmsg("[Msg]   Total $max_cnt records.\n"); 
	if ( $max_cnt <= $opts{'rand_slct'} ) {
		&_write_bnx_records( $opts{'outFH'}, $bnx_data->{'records'} ); 
		#for my $lr ( @{$bnx_data->{'records'}} ) {
		#	chomp($lr); print {$opts{'outFH'}} "$lr\n"; 
		#}
		return; 
	}
	my @use_idx; 
	if ( $opts{'rand_slct'} < $max_cnt/2 ) {
		@use_idx = @{ &_rand_number( $max_cnt, $opts{'rand_slct'} ) }; 
	} else {
		my %t = map { $_ => 1 } @{ &_rand_number( $max_cnt, $max_cnt - $opts{'rand_slct'} ) }; 
		@use_idx = grep { !defined $t{$_} } ( 0 .. ($max_cnt-1)); 
	}
	&_write_bnx_records( $opts{'outFH'}, [ @{$bnx_data->{'records'}}[@use_idx] ] ); 
	#for my $ti ( @use_idx ) {
	#	chomp( $bnx_data->{'records'}[$ti] ); 
	#	print {$opts{'outFH'}} $bnx_data->{'records'}[$ti]."\n"; 
	#}
	return; 
}# rand_slct() 


####################################################################################################
#      Inner sub-routines. 
####################################################################################################
sub _filter_mol_length {
	# $_[0] : $bnx_data->{'records'}
	# $_[1] : min_length_allowed. 
	$_[1] //= 0; 
	my @back; 
	for my $lr (@{$_[0]}) {
		my @l_a = split(/\n/, $lr); 
		my @ta0 = split(/\t/, $l_a[0]); 
		$ta0[2] >= $_[1] or next; 
		push(@back, $lr); 
	}
	return(\@back); 
}
sub _write_bnx_header {
	# $_[0] : Out_FH
	# $_[1] : $bnx_data->{'header'}
	my $a=$_[1]; 
	chomp($a); 
	print {$_[0]} "$a\n"; 
	return; 
}
sub _write_bnx_records {
	# $_[0] : Out_FH
	# $_[1] : $bnx_data->{'records'}
	my $b; 
	for my $a (@{$_[1]}) {
		$b = $a; 
		chomp($b); 
		print {$_[0]} "$b\n"; 
	}
	return; 
}

sub _rand_number {
	# $_[0] : Select numbers from [0-($_[0]-1)] 
	# $_[1] : Number of numbers to select 
	$_[1] //= 1; 
	my %h; 
	my $cnt = 0; 
	while ($cnt < $_[1]) {
		my $n = int(rand($_[0])); 
		defined $h{$n} and next; 
		$h{$n} = 1; 
		$cnt ++; 
	}
	return([sort { $a <=> $b } keys %h]); 
}

sub _load_whole_bnx {
	# $_[0] : file name 
	my $fh = &openFH($_[0], '<'); 
	my %back; 
	while (<$fh>) {
		if (m/^#/) {
			$back{'header'} .= $_; 
			next; 
		}
		$_ =~ m/^0\s/ or &stopErr("[Err] Bad first line of molecule record: $_\n"); 
		push(@{$back{'records'}}, $_); 
		$back{'records'}[-1] .= <$fh>; 
		$back{'records'}[-1] .= <$fh>; 
		$back{'records'}[-1] .= <$fh>; 
	}
	close($fh); 
	return(\%back); 
}# _load_whole_bnx() 

