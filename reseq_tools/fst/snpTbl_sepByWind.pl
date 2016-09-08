#!/usr/bin/perl
# 2015-10-19 
#   Workflow: 
#    1. Set windows according to [Win_Chr, Win_Start, Win_End, Win_Len, Win_Step]; 
#    2. Divide the whole input table to small tables, each table stands for a window; 
#   should setup a basic window list. 
# 2016-09-08
#   Add directory to store all window.snp instead of temporary dir. 
use strict; 
use warnings; 

# An example of how to use this. 
use Parallel::ForkManager; 
#my $MAX_PROCESSES = 10; 
#my $pm = new Parallel::ForkManager($MAX_PROCESSES); 
#for my (my $i=0; $i<100; $i++) {
#	my $pid = $pm->start and next; 
#	$pm->finish; 
#}
#$pm->wait_all_children; 


use Getopt::Long; 

use fileSunhh; 
use LogInforSunhh; 
use mathSunhh; 
my $ms_obj = mathSunhh->new(); 

############################################################
# Basic settings 
############################################################
my %opts; 
my (%wind, @chrIDs); 


############################################################
# Main 
############################################################
&tsmsg("[Rec] Start.\n"); 

&set_opts(); 
&load_snp_tbl();   # Setup $opts{'_inner'}{'tbl_lines'} : recording [ $line1, $line2, ... ]; 
&reform_snp_tbl(); # Reformat $opts{'_inner'}{'tbl_lines'} to format : [ [chrID, posVal, raw_line] ] ; Record $opts{'_inner'}{'max_pos'}{$cur_chr} values ; sort if needed. 
&setup_windows(); # Setup ($opts{'_inner'}{'geno_cols'}, %wind); 
my $dvd_file_list = &dvd_snp_tbl(); 
&write_list(); 
# &del_tmp(); 

&tsmsg("[Rec] Done.\n"); 

############################################################
# Sub-routines 
############################################################

sub usage {
	print <<HH; 
############################################################
# perl $0 -snp_tbl input_SNP.tbl 
# 
# -help 
# 
# -snp_tbl          [filename] 
# -skipSort         [Boolean] Skip sorting table if given. 
#                       In case 1, 'A' or 'W' will be kept, 'D' will be changed to 'N'. 
#                       This affects the markers' genotype, as well as the population size. 
# -ncpu             [1] 
#
# -out              [*STDOUT] Output file storing the values per window. 
#                   
# -wind_start       [1]
# -wind_end         [99999999]
# -wind_end_useMax  [Boolean] If given, the max positions of each chr are recorded, 
#                     and the windows larger than this value are deleted. 
# -wind_length      [10000]
# -wind_step        [1000]
# -chr_colN         [0] column of chromosome. 
# -pos_colN         [1] column of position on chromosome. 
# -skipHN           [0] Skip header lines number. 
# -skipHead         [Boolean] Will add lines in -skipHN if not given. 
# -geno_col         [''] By default I will use the columns in the first reading line except chr_colN and pos_colN; 
#
# -tmp_dir          [''] Default is empty. If it is not empty, this directory is used to store all windows. 
# -rm_oldDir        [Boolean]
############################################################
# Output : 
#  -out   Format: 
#    window1_file_name \\t wind1_chrID \\t wind1_Start \\t wind1_End \\t wind1_length
#    window2_file_name \\t wind2_chrID \\t wind2_Start \\t wind2_End \\t wind2_length
#    .... 
#  Format of window1_file_name : (-skipHN is required)
#    chr \\t pos \\t base \\t indv1 \\t indv2 ... 
#
#
############################################################
HH
	exit 1; 
}

# Set options from the input; 
sub set_opts {
	&tsmsg("[Msg] Setup options.\n"); 
	GetOptions(\%opts, 
		"help!", 
		"snp_tbl:s", 
		"skipSort!", 
		"out:s", 
		"ncpu:i", 
		"wind_start:i", "wind_end:i", "wind_end_useMax!", "wind_length:i", "wind_step:i", 
		"chr_colN:i", "pos_colN:i", 
		"skipHN:i", "skipHead!", 
		"geno_col:s", 
		"tmp_dir:s", # ''
		"rm_oldDir!", 
	); 
	$opts{'help'} and &usage(); 
	-t and !(defined $opts{'snp_tbl'}) and &usage(); 
	&prepare_input_fh(); # Setup $opts{'_inner'}{'inFh'} and skip $opts{'skipHN'} lines;
	&prepare_output_fh(); # Setup $opts{'_inner'}{'outFh'} 
	
	$opts{'ncpu'} //= 1; 
	$opts{'wind_start'} //= 1; 
	$opts{'wind_end'}   //= 99999999; 
	$opts{'wind_length'} //= 10000; 
	$opts{'wind_step'} //= 1000; 
	$opts{'wind_end_useMax'} //= undef(); 
	$opts{'skipSort'} //= undef(); 
	$opts{'chr_colN'}  //= 0; 
	$opts{'pos_colN'}  //= 1; 
	$opts{'geno_col'} //= ''; 
	$opts{'tmp_dir'}  //= ''; 
	$opts{'_inner'}{'geno_cols'} = [ &mathSunhh::_parseCol( $opts{'geno_col'} ) ]; 
	# Usable keys: 
	#  '_inner':'inFh'  => file_handle of input. 
	#  'wind_start|end|length|step' => global parameters to setup windows. 
	#  'wind_end_useMax' => False. If True, will trim windows by CHR's maximum position; 
	#  'chr|pos_colN'    => columns for indexing. 
	#  '_inner':'geno_cols' => [@cols_to_use]
}# sub set_opts() 

# Setup $opts{'_inner'}{'inFh'} and skip $opts{'skipHN'} lines; 
sub prepare_input_fh {
	&tsmsg("[Msg] Prepare input file.\n"); 
	$opts{'_inner'}{'inFh'} = \*STDIN; 
	defined $opts{'snp_tbl'} and $opts{'_inner'}{'inFh'} = &openFH($opts{'snp_tbl'}, '<'); 
	$opts{'skipHN'} //= 0; 
	for ( 1 .. $opts{'skipHN'} ) {
		$opts{'_inner'}{'header'} .= readline($opts{'_inner'}{'inFh'}); 
	}
	return 0; 
	#  'skipHN'   => number of lines to skip from head of file. 
}# sub prepare_input_fh () 

# Setup $opts{'_inner'}{'outFh'} 
sub prepare_output_fh {
	$opts{'_inner'}{'outFh'} = \*STDOUT; 
	defined $opts{'out'} and $opts{'_inner'}{'outFh'} = &openFH($opts{'out'}, '>'); 
	return 0; 
}# prepare_output_fh() 

# Setup $opts{'_inner'}{'tbl_lines'} : recording [ $line1, $line2, ... ]; 
sub load_snp_tbl {
	&tsmsg("[Msg] Loading data.\n"); 
	@{ $opts{'_inner'}{'tbl_lines'} } = readline( $opts{'_inner'}{'inFh'} ) ; 
	chomp(@{$opts{'_inner'}{'tbl_lines'}}); 
	my $nn = scalar( @{$opts{'_inner'}{'tbl_lines'}} ); 
	&tsmsg("[Msg] $nn lines data loaded.\n"); 
	return 0; 
}# load_snp_tbl () 

# Reformat $opts{'_inner'}{'tbl_lines'} to format : [ [chrID, posVal, raw_line] ]. 
# Record $opts{'_inner'}{'max_pos'}{$cur_chr} values ; 
# Setup @chrIDs ; 
#  and sort it if needed. 
sub reform_snp_tbl {
	&tsmsg("[Msg] Reform snp table data.\n"); 
	for ( @{$opts{'_inner'}{'tbl_lines'}} ) {
		my @ta = split(/\t/, $_); 
		my $cur_chr = $ta[ $opts{'chr_colN'} ]; 
		my $cur_pos = $ta[ $opts{'pos_colN'} ]; 
		$_ = [ $cur_chr, $cur_pos, $_ ]; 
		defined $opts{'_inner'}{'max_pos'}{$cur_chr} or do { push(@chrIDs, $cur_chr); $opts{'_inner'}{'max_pos'}{$cur_chr} = $cur_pos; }; 
		$opts{'_inner'}{'max_pos'}{$cur_chr} < $cur_pos and $opts{'_inner'}{'max_pos'}{$cur_chr} = $cur_pos; 
	}
	unless ( $opts{'skipSort'} ) {
		@{$opts{'_inner'}{'tbl_lines'}} = sort my_sort_01 @{$opts{'_inner'}{'tbl_lines'}} ; 
	}
	return 0; 
}# sub reform_snp_tbl () 

sub my_sort_01 {
	no strict; 
	no warnings; 
	my $result = ($a->[0] <=> $b->[0]) || ($a->[0] cmp $b->[0]) || ($a->[1] <=> $b->[1]) ; 
	return $result; 
}

# Setup ($opts{'_inner'}{'geno_cols'}, %wind); 
sub setup_windows {
	&tsmsg("[Msg] Setting up windows.\n"); 
	if ( scalar(@{$opts{'_inner'}{'geno_cols'}}) == 0 ) { 
		# There is no 'geno_cols' defined. 
		my $ln = $opts{'_inner'}{'tbl_lines'}[0][2]; 
		my @ta = split(/\t/, $ln); 
		for (my $i=0; $i<@ta; $i++) {
			$i == $opts{'chr_colN'} and next; 
			$i == $opts{'pos_colN'} and next; 
			push(@{$opts{'_inner'}{'geno_cols'}}, $i); 
		}
	}
	scalar(@{$opts{'_inner'}{'geno_cols'}}) == 0 and &stopErr("[Err] No good geno_cols available. Try use -geno_col\n"); 
	
	# Setup %wind; {chr}
	#for (@{$opts{'_inner'}{'tbl_lines'}}) {
	#	my $cur_chr = $_->[0]; 
	#	my $cur_pos = $_->[1]; 
	#	defined $opts{'_inner'}{'max_pos'}{$cur_chr} or do { push(@chrIDs, $cur_chr); $opts{'_inner'}{'max_pos'}{$cur_chr} = $cur_pos; }; 
	#	$opts{'_inner'}{'max_pos'}{$cur_chr} < $cur_pos and $opts{'_inner'}{'max_pos'}{$cur_chr} = $cur_pos; 
	#}
	for my $cur_chr ( @chrIDs ) {
		my $max_len = ( $opts{'wind_end_useMax'} ) ? $opts{'_inner'}{'max_pos'}{$cur_chr} : $opts{'wind_end'} ; 
		&tsmsg("[Msg]   Setting windows for $cur_chr [ $opts{'wind_start'} - $max_len ]\n"); 
		$wind{$cur_chr} = $ms_obj->setup_windows(
		  'ttl_start'   =>  $opts{'wind_start'}, 
		  'ttl_end'     =>  $max_len, 
		  'wind_size'   =>  $opts{'wind_length'}, 
		  'wind_step'   =>  $opts{'wind_step'}, 
		  'minRatio'    =>  0, 
		); 
	}
	return 0; 
}# setup_windows() 

# Setup smaller snp_tbl files ( @{$opts{'_inner'}{'tmp_wind_file'}} ). 
#  $opts{'_inner'}{'tmp_dir'}
#  $opts{'_inner'}{'windFN2windTI'} : {wind_filename}=>[chrID, chr_wind_idx]
#  $opts{'_inner'}{'chrIdx2fIdx'}   : {chrID}{wind_idx} => file_idx 
sub dvd_snp_tbl {
	&tsmsg("[Msg] Dividing windows.\n"); 
	my $tmpDir; 
	if ( defined $opts{'tmp_dir'} and $opts{'tmp_dir'} ne '' ) {
		$opts{'_inner'}{'tmp_dir'} = $opts{'tmp_dir'}; 
	} else {
		$opts{'_inner'}{'tmp_dir'} = &fileSunhh::new_tmp_dir(); 
		defined $opts{'_inner'}{'tmp_dir'} or &stopErr("[Err] failed to find a temporary directory.\n"); 
	}
	$tmpDir = $opts{'_inner'}{'tmp_dir'}; 
	$opts{'rm_oldDir'} and -d $tmpDir and &fileSunhh::_rmtree($tmpDir); 
	-d $tmpDir or mkdir($tmpDir); 
	my %used; 
	for ( my $ln=0; $ln<@{$opts{'_inner'}{'tbl_lines'}}; $ln++ ) { 
		$ln % 500e3 == 1 and &tsmsg("[Msg] Processed $ln line.\n"); 
		my $cur_chr = $opts{'_inner'}{'tbl_lines'}[$ln]->[0]; 
		my $cur_pos = $opts{'_inner'}{'tbl_lines'}[$ln]->[1]; 
		my (@wind_i) = @{ $ms_obj->map_windows( 'posi'=>$cur_pos, 'wind_hash'=>$wind{$cur_chr} ) }; 
		for my $ti ( @wind_i ) {
			my $file_idx; 
			if ( defined $opts{'_inner'}{'chrIdx2fIdx'}{$cur_chr}{$ti} ) {
				$file_idx = $opts{'_inner'}{'chrIdx2fIdx'}{$cur_chr}{$ti}; 
			} else {
				$file_idx = $ms_obj->newNumber(); 
				$opts{'_inner'}{'chrIdx2fIdx'}{$cur_chr}{$ti} = $file_idx; 
			}
			my $wind_fname = "$tmpDir/wind_${file_idx}"; 
			$opts{'_inner'}{'windFN2windTI'}{$wind_fname} //= [$cur_chr, $ti]; 
			$opts{'_inner'}{'windFN2windTI'}{$wind_fname}[1] == $ti or &stopErr("[Err] bug here.\n"); 
			push( @{$opts{'_inner'}{'windFN2LineN'}{$wind_fname}}, $ln ); 
			defined $used{$wind_fname} or do { push(@{$opts{'_inner'}{'tmp_wind_file'}}, $wind_fname); $used{$wind_fname} = 1; }; 
		}
	}
	my $nn = scalar( @{$opts{'_inner'}{'tmp_wind_file'}} ); 
	&tsmsg("[Msg] Total $nn windows to process.\n"); 

	my $batch_grp = &mathSunhh::dvd_array( $opts{'_inner'}{'tmp_wind_file'}, $opts{'ncpu'}, 1 ); 
	my $pm = new Parallel::ForkManager( $opts{'ncpu'} ); 
	for my $grp ( @$batch_grp ) {
		my $pid = $pm->start and next; 
		for my $inFname (@$grp) {
			open O,'>',"$inFname" or &stopErr("[Err] Failed to write to [$inFname]\n"); 
			$opts{'skipHead'} or do { defined $opts{'_inner'}{'header'} and print O $opts{'_inner'}{'header'}; }; 
			for my $tr ( @{ $opts{'_inner'}{'tbl_lines'} }[ @{ $opts{'_inner'}{'windFN2LineN'}{$inFname} } ] ) { 
				print O $tr->[2]."\n"; 
			} 
			close O; 
		}
		$pm->finish; 
	}
	$pm->wait_all_children; 
	return $opts{'_inner'}{'tmp_wind_file'}; 
}# dvd_snp_tbl () 

sub del_tmp {
	&tsmsg("[Msg] Delete temporary Dir: $opts{'_inner'}{'tmp_dir'}\n"); 
	&fileSunhh::_rmtree( $opts{'_inner'}{'tmp_dir'} ); 
	return 0; 
}

sub write_list {
	for (@{$opts{'_inner'}{'tmp_wind_file'}}) {
		my ($cur_chr, $cur_si) = @{$opts{'_inner'}{'windFN2windTI'}{$_}}; 
		print { $opts{'_inner'}{'outFh'} } join("\t", $_, $cur_chr, @{ $wind{$cur_chr}{'loci'}{$cur_si} } )."\n"; 
	}
	return 0; 
}

############################################################
############################################################

=head1 _write_fst_R( $outR_filename ) 
=cut
sub _write_fst_R {
	my $fh = &openFH($_[0], '>'); 
	print {$fh} <<'HH';

argvs <- commandArgs( trailingOnly=TRUE ) ;
flist <- read.table( file = argvs[1], stringsAsFactors=F, colClasses=c('character'), header=F ) ;
# The first column is input file name, the second column is output file prefix.
library(hierfstat);

for ( i in 1:nrow(flist) ) {
        # input file    : flist[i,1]
        # output prefix : flist[i,2]
        aa <- read.table( file= flist[i,1], header=T, colClasses="numeric", stringsAsFactors=F )
        aa.stats <- basic.stats( aa )
        write.table( aa.stats$perloc,  file=paste0(flist[i, 2], ".perSite", sep=""), append=F, row.names=T, col.names=NA, quote=F, sep="\t" )
        write.table( aa.stats$overall, file=paste0(flist[i, 2], ".perWind", sep=""), append=F, row.names=T, col.names=NA, quote=F, sep="\t" )
}

HH
	close($fh); 
	return ($_[0]); 
}# _write_fst_R() 

