#!/bin/env perl
# Author email : hs738@cornell.edu or biosunhh@gmail.com

# 2006-11-17 13:18 ׼������һ��best_uniq��������������*.psl��ʽ�ļ�����δ��ɣ�
# 2006-11-24 13:33 ����һ��С����,UniqColLine,finished.
# 2006-11-24 15:16 ����best_uniq
# 2007-1-17 17:16 ��΢������best_uniqCol, ���best_rep����;
# 2007-06029 add a para -label_mark
# 2007-0711  add a para -col_repCount
# 2007-9-20 10:58:10 add a output for -col_stat; 
# 2007-11-02 09:34 correct a bug in best_uniq function. 
# 2012-09-05 17:21 Add a cbind function like in R language. 
# 2012-11-07 11:15 Change '-column' option to fit a input type like "0-10". 
# 2013-02-08 10:54 Add '-best_disOrdCol' option to compare 'A\tB' line with 'B\tA' line. 
# 2013-02-22 09:58 Writing a new sub-function to cluster given names from a list. 
# 2013-04-08 Edit cluster_lines, and add cluster_cols. 
# 2013-07-31 Add a function to trim end blanks from Windows. 
# 2013-08-19 Add a function to fill NULL cells with assigned character. 
# 2013-08-21 Add a function to transpose a matrix table. Also fix a bug for function &fillNull() which does not fill the first column of line. 
# 2013-09-11 Edit -cbind log infor. 
# 2013-11-26 Add in uniqComb.pl function. Add time tag for message output. 
# 2015-04-21 Add -colByTbl to select columns according to an input table. 
# 2015-07-14 Add -spec_loci_1from2 to find region in -spec_loci_f1 (Format: ID\tstart\tEnd) that is absent in -spec_loci_f2 (Format: ID\tStart\tEnd). 
# 2016-03-09 Add -col_sort_rule for more flexible. 
# 2016-03-11 Finally, I accept to add my own perl modules into this frequently used script. Change all split to &splitL(); 
# 2016-10-31 Add -colByTbl_cut to use 'cut & paste' method. 
# 2016-12-01 Add -cpuN for selected, time consuming tasks. 
# 2016-12-02 'cut & paste' method for -colByTbl_cut is not fast enough. trying to use -cpuN
# 2016-12-03 -colByTbl_cut is removed. 
# 2020-08-23 Fix -col_sort for mixed character and number column; 

use strict;
use warnings; 
use fileSunhh; 
use LogInforSunhh; 

use Getopt::Long;
my %opts;
GetOptions(\%opts,
	"combine!","column:s","kick_col!", 
	"colByTbl:s", "colByTbl_also:s", 
	"max_col:s","min_col:s",
	"skip:i",
	"reverse!",
	"col_sort:s", "col_sort_rule:s", 
	"col_uniq:s","col_reps:s","col_repCount:s", # uniq/repeat columns pick.
	"col_stat:i", "col_stat_asINS!", 
	"symbol:s",
	"best_uniqCol:s","select_col:s","select_rule:s","best_rep:s","best_disOrdCol!", 
	"UniqColLine:s", "topN:i", 
	"label_mark:s",
	"cbind!", "chkColNum!", 
	"cluster_lines!", "cluster_cols:s", 
	"trimEndReturn!", 
	"fillNull:s", 
	"transpose!", "skip_null_line!", "beMatrix!", "fill_new:s",      # transpose a matrix. 
	"kSrch_idx:s","kSrch_idxCol:s","kSrch_srcCol:s","kSrch_drop!", "kSrch_line!", # Similar to linux command join, without joining and with more index columns. Combined from uniqComb.pl 
	"dR2dN!", 
	"col_head!", 
	"chID_RefLis:s", "chID_Row!", "chID_OldColN:i", "chID_NewColN:i", "chID_skipH:i", "chID_RowColN:i", # Change Column/Row names according to reference. 
	"spec_loci_1from2!", "spec_loci_f1:s", "spec_loci_f2:s", "spec_loci_minLen:i", 
	"self_agp_from_keyLen!", # output of "deal_fasta.pl -attribute key:len"
	"join_sameRow!", 
	"cpuN:i", 
	"log_ln:i", 
	"help!");

sub usage {

	my $info=<<INFO;
##################################################
command:perl $0 <STDIN|parameters>
# 2016-10-31

  -help           help infomation;
  -combine        combine files;
  -column<int>    cols "num1,num2,num3..." will be picked out and joined to a new line;
   -kick_col      [Boolean] If given, the columns defined by -column will be removed instead of extracted. 
  -colByTbl<Str>  [index_file]. The 1st column of index_file is used to match the first line of input.table. Output matching columns from input.table. 
   -colByTbl_also [Str] Cols "num1,num2,num3" will be output before checking the matching between input.table and index_file. 
  -max_col        Similar to -column, compare by order;
  -min_col        Similar to -column, compare by order;
  -skip<int>      Skip first <int> lines.Usually for skipping head lines;
  -col_sort       [num,num,num,...]sort by col (0-based).
  -col_sort_rule  [-1,-1,-1,...] -1 for small to max, 1 for max to small. 

  -col_uniq       pick lines with its columns uniq
  -col_reps       pick lines with its columns repeat
  -col_repCount   Show col repeat times.

  -col_stat<int>  Get Sum, Mean, Median, Min, Max of the Column;
  -col_stat_asINS Add ins_calc() functioin to calculate Mean like bwa. 
  -reverse        reverse the file lines order.

  -best_uniqCol   1st column is 0; col_1[,col_2...]
  -select_col     1st column is 0; col_1[,col_2...]
  -select_rule    1/-1[,1/-1...]; 1 for Max, -1 for Min;
  -best_rep       out a best uniq cols' repeat records if given;
  -best_disOrdCol Do not care the order of uniqCols, so "A\\tB" line will be compared with "B\\tA" line. 

  -UniqColLine    output a line with its col=UniqColLine one time.<*,*,*...>
  -topN           [1] Get the first topN number of lines. 

  -label_mark     Marks to be added to label_col. in format \" 001..090::a..z::1[::...]\"

  -cbind          Bind input files by columns. 
  -chkColNum      Check column number for each line if given. 

  -cluster_lines  Cluster genes in the same line into a same group. 
  -cluster_cols   Set columns to be clustered. 

  -trimEndReturn  Trim end return from windows lines. 

  -fillNull       Fill NULL cells with character assigned here. [String]. 

  -transpose      Transpose a matrix table. 
  -skip_null_line Will not read in empty lines (with only blanks) if assigned. 
  -beMatrix       The output table will be filled with NULL blank (empty) to be an intact (M x N) matrix. 
  -fill_new       Fill new added cells in output matrix with character assigned here. [String]. 
                    Be aware that only the new added cells will be filled with this character! 
                    So the raw NULL cells will remain unchanged! 

  '-kSrch_idx:s','-kSrch_idxCol:s','-kSrch_srcCol:s','-kSrch_drop!', '-kSrch_line!', # Similar to linux command join, without joining and with more index columns. Combined from uniqComb.pl 

  # Change Column/Row names according to reference. Need -chID_RefLis . 
  -chID_RefLis      [filename] With OldCol and NewCol for name conversion; 
    -chID_OldColN   [ColN] Column number of OldCol in RefLis
    -chID_NewColN   [ColN] Column number of NewCol in RefLis
  -chID_Row         [Boolean]  Together with -chID_RowColN for replacing IDs in specified column_N of each row. 
    -chID_skipH     [Number]  If >= 1, the heading Number of lines will not be changed. 
    -chID_RowColN   [Number]  The column_N that needs to be changed. 
  

  '-spec_loci_1from2!', '-spec_loci_f1:s', '-spec_loci_f2:s' , # Find region in -spec_loci_f1 that is absent in -spec_loci_f2. 
                  # Format of '-spec_loci_f1/2' : ID \\t Start \\t End \\n
    '-spec_loci_minLen:i', # Minimum length of specific region. Default is 0; 

  -self_agp_from_keyLen   [Boolean] Input should be output of "deal_fasta.pl -attr key:len"
                            This could be in format "Key\\tLen\\n" or format "Key\\tLen\\tNew_ID\\n"

  -symbol         Defining the symbol to divide data.Default is "\\t";
  -dR2dN          [Boolean] Change \\r to \\n in files. 

  -join_sameRow   [Boolean] Join col-1 accoding to col-0; 

  -log_ln         [0]

  -cpuN           [1] Multi-threading. Available for : 
                    -column
                    -chID_RefLis
                    -fillNull
                    -colByTbl
                    -kSrch_idx
##################################################
INFO

	print STDOUT "$info"; 
	&stop("[Err] Exit $0\n"); 
	exit(1); 
}

if ($opts{help} or (-t and !@ARGV) and !$opts{spec_loci_1from2}) {
	&usage; 
}

#****************************************************************#
#--------------Main-----Function-----Start-----------------------#
#****************************************************************#

my %glob; 

# Making File handles for reading;
our @InFp = () ;  
if ( !@ARGV )
{
	-t or @InFp = (\*STDIN);
}
else
{
	for (@ARGV) {
		push( @InFp, &openFH($_,'<') );
	}
}


my $symbol = "\t";
&goodVar($opts{symbol}) and $symbol = $opts{symbol};

$opts{'log_ln'} //= 0; 
$opts{'cpuN'}   //= 1; 
$opts{'cpuN'} = int( $opts{'cpuN'} ); 
$opts{'cpuN'} < 1 and $opts{'cpuN'} = 1; 


#****************************************************************#
#--------------Main-----Use multi-processes----------------------#
#****************************************************************#
my $pm; 
if ( $opts{'cpuN'} > 1 ) {
	use Parallel::ForkManager; 
	$pm = new Parallel::ForkManager($opts{'cpuN'}); 
}


#****************************************************************#
#--------------Main-----Invoke functions-t-----------------------#
#****************************************************************#



&reverse_lines() if($opts{reverse});
&skip() if (&goodVar($opts{skip}) and $opts{skip}>0);
&combine() if ($opts{combine});
&column() if ( &goodVar($opts{column}) ); # Should have; 
&extreme() if ( &goodVar($opts{max_col}) || &goodVar($opts{min_col}) );
&colStat() if ( &goodVar($opts{col_stat}) );
if ( &goodVar($opts{col_sort}) ) {

	&_prep_cSort_ruls(); 

	my @file;
	for my $fh (@InFp) {
		while (<$fh>) {
			chomp; 
			push(@file, $_); 
		}
	}# End $fh @InFp
	my @output = sort col_sort @file; 
	foreach my $line (@output) {
		print STDOUT "$line\n"; 
	}
}# end if

&uniq_rep() if ( &goodVar($opts{col_uniq}) || &goodVar($opts{col_reps}) || &goodVar($opts{col_repCount}) );

&best_uniqCol() if ( &goodVar($opts{best_uniqCol}) );

&UniqColLineF() if ( &goodVar($opts{UniqColLine}) );

&LabelTbl($opts{label_mark}) if ( &goodVar($opts{label_mark}) );

&cbind() if ($opts{cbind}); 

&cluster_group() if ($opts{cluster_lines}); 

&trimEndReturn() if ( $opts{trimEndReturn} ); 

&fillNull() if ( defined $opts{fillNull} ); 

&t() if ( $opts{transpose} ); 

&kSrch() if ( &goodVar($opts{kSrch_idx}) ); 

&dR2dN() if ( $opts{dR2dN} ); 

&showHeader() if ( $opts{col_head} ); 

&chRowColName() if ( defined $opts{chID_RefLis} ); 

&colByTbl() if ( defined $opts{'colByTbl'} ); 

&specLoci() if ( defined $opts{'spec_loci_1from2'} ); 

&kl2agp() if ( defined $opts{'self_agp_from_keyLen'} ); 

&joinColByCol() if ( $opts{'join_sameRow'} ); 

for (@InFp) {
	close ($_);
}


######################################################################
## sub-routines for functions. 
######################################################################

sub joinColByCol {
	my %h; 
	my @k_cols = (0); 
	my @v_cols = (1); 
	for my $fh (@InFp) {
		while (&wantLineC($fh)) {
			my @ta = &splitL($symbol, $_); 
			my $kk = join("\t", @ta[@k_cols]); 
			my $vv = join(";;", @ta[@v_cols]); 
			push(@{$h{$kk}}, $vv); 
		}
	}
	for my $kk (sort keys %h) {
		print STDOUT join("\t", $kk, join(', ', @{$h{$kk}}))."\n"; 
	}
	return; 
}# joinColByCol() 

sub kl2agp {
	for my $fh (@InFp) {
		my $has_newID = 0; 
		my $cnt = 0; 
		while (&wantLineC($fh)) {
			my @ta = &splitL($symbol, $_); 
			$ta[0] =~ m/^key$/i and next; 
			my $new_ID = $ta[0]; 
			$cnt ++; 
			$cnt == 1 and defined $ta[2] and $has_newID = 1; 
			if ( $has_newID == 1 ) {
				defined $ta[2] or &stopErr("[Err] Failed to get new_ID at line : $_\n"); 
				$new_ID = $ta[2]; 
			}
			print STDOUT join("\t", $new_ID, 1, $ta[1], 1, "W", $ta[0], 1, $ta[1], '+')."\n"; 
		}
	}
}# kl2agp () 


# 
sub specLoci {
	$opts{'spec_loci_minLen'} //= 0; 
	my $f1h = &openFH($opts{'spec_loci_f1'}, '<'); 
	my $f2h = &openFH($opts{'spec_loci_f2'}, '<'); 
	# ID \\t Start \\t End \\n
	my %loc1 = %{ &loc_tbl2hash( $f1h ) }; 
	my %loc2 = %{ &loc_tbl2hash( $f2h ) }; 

	for my $id1 (sort keys %loc1) {
		my @spec_in_1; 
		for my $ar1 ( @{$loc1{$id1}} ) {
			my ($s1, $e1) = @$ar1; 
			my $cur_p = $s1; 
			for my $ar2 ( @{$loc2{$id1}} ) {
				my ($s2, $e2) = @$ar2; 
				if ( $e2 < $cur_p ) {
					next; 
				} elsif ( $s2 > $e1 ) { 
					last; 
				} elsif ( $s2 <= $cur_p ) {
					$cur_p = $e2+1; 
				} elsif ( $s2 > $cur_p ) {
					push(@spec_in_1, [$cur_p, $s2-1]); 
					$cur_p = $e2+1; 
				}
			}
			if ( $cur_p <= $e1 ) {
				push(@spec_in_1, [$cur_p, $e1]); 
			}
		}
		for my $ar (@spec_in_1) {
			$ar->[1]-$ar->[0]+1 >= $opts{'spec_loci_minLen'} or next; 
			print STDOUT join("\t", $id1, $ar->[0], $ar->[1])."\n"; 
		}
	}
}# specLoci() 

# Extract columns according to $opts{'colByTbl'}; 
sub colByTbl {
	my @col_also = (); 
	defined $opts{'colByTbl_also'} and @col_also = &parseCol($opts{'colByTbl_also'}); 
	my $idxFh = &openFH($opts{'colByTbl'}, '<'); 
	my (%needColID, $cnt); 
	$cnt = 0; 
	while (<$idxFh>) {
		chomp; 
		my @ta = &splitL($symbol, $_); 
		defined $needColID{$ta[0]} or do { $needColID{$ta[0]} = $cnt; $cnt ++; }; 
	}
	close($idxFh); 


	for my $fh ( @InFp ) {
		# header_txt
		my $header_txt = ''; 
		my @cur_ColNs = @col_also; 
		{
			my $first_line = <$fh>; 
			my %cur_ID2ColN; 
			chomp($first_line); 
			my @ta = &splitL($symbol, $first_line); 
			for (my $i=0; $i<@ta; $i++) {
				defined $needColID{$ta[$i]} and $cur_ID2ColN{$ta[$i]} = $i; 
			}
			for (sort { $needColID{$a} <=> $needColID{$b} } keys %needColID) {
				if ( defined $cur_ID2ColN{$_} ) {
					push(@cur_ColNs, $cur_ID2ColN{$_}); 
				} else {
					push(@cur_ColNs, 'NA'); 
				}
			}
			$header_txt = join("\t", map { ( $_ eq 'NA' ) ? '' : $ta[$_] ; } @cur_ColNs)."\n"; 
		}
		if ( defined $pm ) {
			$opts{'cpuN'} > 1 or &stopErr("[Err] Check 01 in column()\n"); 
			# separate files 
			my $wrk_dir = &fileSunhh::new_tmp_dir( 'create' => 1 ); 
			my @sub_fn = &fileSunhh::dvd_file( $fh, $opts{'cpuN'}, 'keep_order' => 1, 'with_header' => 0, 'sub_pref' => "$wrk_dir/sub_", 'tmpFile' => "$wrk_dir/base_0" ); 
			# Process sub-files
			for my $sfn ( @sub_fn ) {
				my $pid = $pm->start and next; 
				open F,'<',"$sfn" or &stopErr("[Err] Failed to open subfile [$sfn]\n"); 
				open O,'>',"$sfn.o" or &stopErr("[Err] Failed to write subfile [$sfn.o]\n"); 
				while (<F>) {
					chomp; 
					my @ta = &splitL($symbol, $_); 
					print O join("\t", map { ( $_ eq 'NA' ) ? '' : $ta[$_] ; } @cur_ColNs)."\n"; 
				}
				close O; 
				close F; 
				$pm->finish; 
			}
			$pm->wait_all_children; 
			# Merge sub-files
			print STDOUT $header_txt; 
			for my $sfn ( @sub_fn ) {
				open F,'<',"$sfn.o" or &stopErr("[Err] Failed to open subfile [$sfn.o]\n"); 
				while (<F>) {
					print STDOUT $_; 
				}
				close F; 
			}
			# Delete temp_dir 
			&fileSunhh::_rmtree($wrk_dir); 
		} else {
			print STDOUT $header_txt; 
			while (<$fh>) {
				chomp; 
				my @ta = &splitL($symbol, $_); 
				print STDOUT join("\t", map { ( $_ eq 'NA' ) ? '' : $ta[$_] ; } @cur_ColNs)."\n"; 
			}
			close ($fh); 
		}
	}

}# sub colByTbl() 

# change table's name according to reference list. 
sub chRowColName {
	$opts{chID_OldColN} = $opts{chID_OldColN} // 0; 
	$opts{chID_NewColN} = $opts{chID_NewColN} // 1; 
	$opts{chID_skipH} = $opts{chID_skipH} // 0; 
	my $is_rowID = (defined $opts{chID_Row}) ? 1 : 0 ; 
	$opts{chID_RowColN} = $opts{chID_RowColN} // 0; 

	my $lisFH = &openFH($opts{chID_RefLis}, '<'); 
	my %old2new; 
	while (<$lisFH>) {
		chomp; 
		my @ta = &splitL($symbol, $_); 
		$old2new{ $ta[ $opts{chID_OldColN} ] } = $ta[ $opts{chID_NewColN} ]; 
	}
	close ($lisFH); 
	
	for my $fh ( @InFp ) {
		my $cur_n = -1; 

		if ( defined $pm ) {
			$opts{'cpuN'} > 1 or &stopErr("[Err] Check 02 in chRowColName()\n"); 
			# Keep header 
			my $header_txt = ''; 
			if ( $opts{'chID_skipH'} > 0 ) {
				my $cnt = 0; 
				while (<$fh>) {
					$cnt ++; 
					$header_txt .= $_; 
					$cnt == $opts{'chID_skipH'} and last; 
				}
			}
			unless ( $is_rowID ) {
				my $h_line = <$fh>; 
				chomp($h_line); 
				my @ta = &splitL( $symbol, $h_line ); 
				for my $tb (@ta) {
					$tb = $old2new{$tb} // $tb; 
				}
				$header_txt .= join("\t", @ta)."\n"; 
			}

			my $wrk_dir = &fileSunhh::new_tmp_dir(); 
			mkdir($wrk_dir) or &stopErr("[Err] Failed to create dir [$wrk_dir]\n"); 
			my @sub_fn = &fileSunhh::dvd_file( $fh, $opts{'cpuN'}, 'keep_order' => 1, 'with_header' => 0, 'sub_pref' => "$wrk_dir/sub_", 'tmpFile' =>      "$wrk_dir/base_0" );
			for my $sfn (@sub_fn) {
				my $pid = $pm->start and next; 
				if ( $is_rowID ) {
					open F,'<',"$sfn" or &stopErr("[Err] Failed to open file [$sfn]\n");
					open O,'>',"$sfn.o" or &stopErr("[Err] Failed to open file [$sfn.o]\n"); 
					while (<F>) {
						chomp; 
						my @ta = &splitL( $symbol, $_ ); 
						$ta[$opts{chID_RowColN}] = $old2new{ $ta[$opts{chID_RowColN}] } // $ta[$opts{chID_RowColN}] ; 
						print O join("\t", @ta)."\n"; 
					}
					close O; 
					close F; 
				} else {
					&fileSunhh::_move( $sfn, "$sfn.o" ); 
				}
				$pm->finish;
			}
			$pm->wait_all_children; 

			# Add $header_txt and sub_fn to final one. 
			print STDOUT $header_txt; 
			for my $sfn ( @sub_fn ) {
				open F,'<',"$sfn.o" or &stopErr("[Err] Failed to reopen file [$sfn.o]\n"); 
				while (<F>) {
					print STDOUT $_; 
				}
				close F; 
			}
			&fileSunhh::_rmtree($wrk_dir); 
		} else {
			while (<$fh>) {
				chomp; 
				if ( $. <= $opts{chID_skipH} ) {
					print STDOUT "$_\n"; 
					next; 
				} else {
					$cur_n ++; 
				}
				if ( m!^\s*(##|$)! ) {
					print STDOUT "$_\n"; 
					next; 
				}
				if ( $is_rowID ) {
					my @ta = &splitL($symbol, $_); 
					$ta[$opts{chID_RowColN}] = $old2new{ $ta[$opts{chID_RowColN}] } // $ta[$opts{chID_RowColN}] ; 
					print STDOUT join("\t", @ta)."\n"; 
				} elsif ( $cur_n == 0 ) {
					my @ta = &splitL($symbol, $_); 
					for my $tb (@ta) {
						$tb = $old2new{$tb} // $tb; 
					}
					print STDOUT join("\t", @ta)."\n"; 
				} else {
					print STDOUT "$_\n"; 
				}
			}# End while 
		}

	}#End for 
	return 0; 
}# sub chRowColName () 

# Show header column-number of a file. 
sub showHeader {
	for my $fh ( @InFp ) {
		while (<$fh>) {
			chomp; 
			my @ta = &splitL($symbol, $_); 
			for (my $i=0; $i<@ta; $i++) {
				print STDOUT join("\t", $i, $ta[$i])."\n"; 
			}
			last; 
		}
	}
}#End sub showHeader

# Change \r to \n
sub dR2dN {
	for my $fh (@InFp) {
		while (<$fh>) {
			s/\r/\n/g; 
			print STDOUT $_; 
		}
	}
}#sub dR2dN

# Similar to linux_join and from perl script uniqComb.pl 
## "kSrch_idx:s","kSrch_idxCol:s","kSrch_srcCol:s","kSrch_drop!", # Similar to linux command join, without joining and with more index columns. Combined from uniqComb.pl 
sub kSrch {
	my @src_Cols = ( defined $opts{kSrch_srcCol} ) ? ( &parseCol( $opts{ kSrch_srcCol } ) ) : ( 0 ) ; 
	my @idx_Cols = ( defined $opts{kSrch_idxCol} ) ? ( &parseCol( $opts{ kSrch_idxCol } ) ) : ( 0 ) ; 
	my %idx_key; 
	for my $idxF ( &splitL(",", $opts{'kSrch_idx'}) ) {
		my $fh = &openFH( $idxF, '<' ); 
		my %tmp_cnt = ( 'cntN_base'=>0 , 'cntN_step'=>$opts{'log_ln'} ); 
		if ( $opts{kSrch_line} ) {
			while (<$fh>) {
				$opts{'log_ln'} > 0 and &fileSunhh::log_section( $. , \%tmp_cnt ) and &tsmsg("[Msg] Reading [$.] line in kSrch_idx\n"); 
				$idx_key{$_} = 1; 
			}
		}else{
			while (<$fh>) {
				$opts{'log_ln'} > 0 and &fileSunhh::log_section( $. , \%tmp_cnt ) and &tsmsg("[Msg] Reading [$.] line in kSrch_idx\n"); 
				chomp; 
				my @temp = &splitL($symbol, $_); 
				my $tkey = join("\t", @temp[@idx_Cols]); 
				$idx_key{$tkey} = 1; 
			}
		}
		close ($fh); 
	}# for my $idxF 

	for my $fh (@InFp) {
		my %tmp_cnt = ( 'cntN_base'=>0 , 'cntN_step'=>$opts{'log_ln'} ); 
		if ( $opts{kSrch_line} ) {
			# Using the whole line as a key. 
			if ( defined $pm ) {
				$opts{'cpuN'} > 1 or &stopErr("[Err] cpuN not > 1\n"); 
				# header_txt
				# separate file
				my $wrk_dir = &fileSunhh::new_tmp_dir( 'create' => 1 ); 
				my @sub_fn = &fileSunhh::dvd_file( $fh, $opts{'cpuN'}, 'keep_order' => 1, 'with_header' => 0, 'sub_pref' => "$wrk_dir/sub_", 'tmpFile' => "$wrk_dir/base_0" ); 
				# process sub-file
				for my $sfn ( @sub_fn ) { 
					my $pid = $pm->start and next; 
					open F,'<',"$sfn" or die; 
					open O,'>',"$sfn.o" or die; 
					while (<F>) {
						if ($opts{kSrch_drop}) {
							$idx_key{$_} or print O $_; 
						}else{
							$idx_key{$_} and print O $_; 
						}
					}#End while 
					close O; 
					close F; 
					$pm->finish; 
				}
				$pm->wait_all_children; 
				# merge sub-file
				for my $sfn ( @sub_fn ) {
					open F,'<',"$sfn.o" or die; 
					while (<F>) { print STDOUT $_; }
					close F; 
				}
				# delete temporary
				&fileSunhh::_rmtree($wrk_dir); 
			} else {
				while (<$fh>) {
					$opts{'log_ln'} > 0 and &fileSunhh::log_section( $. , \%tmp_cnt ) and &tsmsg("[Msg] Reading [$.] line in kSrch_src\n"); 
					if ($opts{kSrch_drop}) {
						$idx_key{$_} or print; 
					}else{
						$idx_key{$_} and print; 
					}
				}#End while 
			}
		}else{
			# Using the selected columns as key. 
			if ( defined $pm ) {
				$opts{'cpuN'} > 1 or &stopErr("[Err] cpuN not > 1\n"); 
				# header_txt
				# separate file
				my $wrk_dir = &fileSunhh::new_tmp_dir( 'create' => 1 ); 
				my @sub_fn = &fileSunhh::dvd_file( $fh, $opts{'cpuN'}, 'keep_order' => 1, 'with_header' => 0, 'sub_pref' => "$wrk_dir/sub_", 'tmpFile' => "$wrk_dir/base_0" ); 
				# process sub-file
				for my $sfn ( @sub_fn ) { 
					my $pid = $pm->start and next; 
					open F,'<',"$sfn" or die; 
					open O,'>',"$sfn.o" or die; 
					while (<F>) {
						chomp; 
						my @temp = &splitL($symbol, $_); 
						my $tkey = join("\t", @temp[@src_Cols]); 
						if ($opts{kSrch_drop}) {
							$idx_key{$tkey} or print O "$_\n"; 
						}else{
							$idx_key{$tkey} and print O "$_\n"; 
						}
					}#End while 
					close O; 
					close F; 
					$pm->finish; 
				}
				$pm->wait_all_children; 
				# merge sub-file
				for my $sfn ( @sub_fn ) {
					open F,'<',"$sfn.o" or die; 
					while (<F>) { print STDOUT $_; }
					close F; 
				}
				# delete temporary
				&fileSunhh::_rmtree($wrk_dir); 
			} else {
				while (<$fh>) {
					$opts{'log_ln'} > 0 and &fileSunhh::log_section( $. , \%tmp_cnt ) and &tsmsg("[Msg] Reading [$.] line in kSrch_src\n"); 
					chomp; 
					my @temp = &splitL($symbol, $_); 
					my $tkey = join("\t", @temp[@src_Cols]); 
					if ($opts{kSrch_drop}) {
						$idx_key{$tkey} or print STDOUT "$_\n"; 
					}else{
						$idx_key{$tkey} and print STDOUT "$_\n"; 
					}
				}# End while 
			}
		}# End if kSrch_line 
	}# for my $fh (@InFp)
}# End kSrch 

# Transpose a matrix . Like the t() function in R language. 
sub t {
	my @result_lines; 
	my $line_idx = -1; 
	for my $fh (@InFp) {
		while (<$fh>) {
			$opts{skip_null_line} and /^\s*$/ and next; # Here empty line should not be skipped until demanded. 
			chomp; 
			my @ta = &splitL($symbol, $_); 
			$line_idx ++; 
			for (my $i=0; $i<@ta; $i++) {
				$result_lines[$i][$line_idx] = $ta[$i]; 
			}
		}# End while 
	}# End $fh in @inFp
	my $fill_new = (defined $opts{fill_new}) ? $opts{fill_new} : '' ; 
	for my $r1 (@result_lines) {
		if ($opts{beMatrix}) {
			for (my $i=0; $i <= $line_idx; $i++) { defined $r1->[$i] or $r1->[$i] = $fill_new; } 
		} else {
			for my $ta (@$r1) { defined $ta or $ta = $fill_new; } 
		}
		print STDOUT join("\t", @$r1)."\n"; 
	}
}#End sub t() 


# Fill NULL cells in table with character assigned. 
sub fillNull {
	my $tc = 1; 
	$opts{fillNull} ne '' and $tc = $opts{fillNull}; 
	for my $fh (@InFp) {
		if ( defined $pm ) {
			$opts{'cpuN'} > 1 or &stopErr("[Err] Check 01 in column()\n"); 
			my $wrk_dir = &fileSunhh::new_tmp_dir( 'create' => 1 ); 
			my @sub_fn = &fileSunhh::dvd_file( $fh, $opts{'cpuN'}, 'keep_order' => 1, 'with_header' => 0, 'sub_pref' => "$wrk_dir/sub_", 'tmpFile' => "$wrk_dir/base_0" ); 
			for my $sfn (@sub_fn) {
				my $pid = $pm->start and next; 
				open F,'<',"$sfn" or &stopErr("[Err] Failed to open file [$sfn]\n"); 
				open O,'>',"$sfn.o" or &stopErr("[Err] Failed to open file [$sfn.o]\n"); 
				while (<F>) {
					while (s/$symbol$symbol/$symbol$tc$symbol/og) { 1; } 
					s/^$symbol/$tc$symbol/o; 
					s/$symbol$/$symbol$tc/o; 
					print O $_; 
				}
				close O; 
				close F; 
				$pm->finish; 
			}
			$pm->wait_all_children; 
			for my $sfn ( @sub_fn ) {
				open F,'<',"$sfn.o" or &stopErr("[Err] Failed to open file [$sfn.o]\n"); 
				while (<F>) {
					print STDOUT $_; 
				}
				close F; 
			}
			&fileSunhh::_rmtree($wrk_dir); 
		} else {
			while (<$fh>) {
				while (s/$symbol$symbol/$symbol$tc$symbol/og) { 1; } 
				s/^$symbol/$tc$symbol/o; 
				s/$symbol$/$symbol$tc/o; 
				print; 
			}# End while 
		}
	}# End $fh in @InFp
}# End sub fillNull

# trim end returns from windows. 
sub trimEndReturn {
	for my $fh (@InFp) {
		while (<$fh>) {
			s/[\b\r\n]+$//; 
			print "$_\n"; 
		}# End while 
	}# End $fh in @InFp
}# End sub trimEndReturn 

#cluster_lines , cluster_group
sub cluster_group {
	my (%groupID, $gID); 
	my %id2genes; 
	$gID = 0; 
	my (@tcols, $is_col); 
	$is_col = 0; 
	if (defined $opts{cluster_cols}) {
		@tcols = &parseCol($opts{cluster_cols}); 
	}
	@tcols > 0 and $is_col = 1; 
	for my $fh (@InFp) {
		while (<$fh>) {
			$. % 1000 == 1 and warn "$. lines.\n"; 
			/^\s*$/ and next; 
			chomp; 
			my @ta = &splitL($symbol, $_); 
	
			$is_col == 1 and @ta = @ta[@tcols]; 
			my @ta1; 
			for my $tid (@ta) {
				(defined $tid and $tid ne '') or next; 
				push(@ta1, $tid); 
			}
			@ta = @ta1; 
			@ta > 0 or next; 
	
			my ($minV, $maxV) = ('X', 'X'); 
			my %usedV; 
			for my $tid (@ta) {
				if (defined $groupID{$tid}) {
					$usedV{ $groupID{$tid} } = 1; 
					$minV = ($minV eq 'X') 
						? $groupID{$tid}
						: ($minV > $groupID{$tid}) 
							? $groupID{$tid}
							: $minV
					; 
					$maxV = ($maxV eq 'X') 
						? $groupID{$tid}
						: ($maxV < $groupID{$tid}) 
							? $groupID{$tid}
							: $maxV
					; 
				}
			}#End for my $tid 
			if ($minV eq 'X') {
				for my $tid (@ta) {
					$groupID{$tid} = $gID; $id2genes{$gID}{$tid} = 1; 
				}
				$gID ++; 
			} elsif ($minV < $maxV) {
				delete $usedV{$minV}; 
				for my $tv (keys %usedV) {
					for my $tk (keys %{$id2genes{$tv}}) {
						$groupID{$tk} = $minV; $id2genes{$minV}{$tk} = 1; 
					}
					delete $id2genes{$tv}; 
				}
				for my $tid (@ta) {
					exists $groupID{$tid} or do { $groupID{$tid} = $minV; $id2genes{$minV}{$tid} = 1; }; 
				}
			} else {
				for my $tid (@ta) {
					exists $groupID{$tid} or do { $groupID{$tid} = $minV; $id2genes{$minV}{$tid} = 1; }; 
				}
			}#End if ( $minV eq 'X' 
		}#End while 
	}# End for $fh in @InFp 
	my (%olines, %linev); 
	while ( my ($tk, $tv) = each %groupID ) {
		push(@{$olines{$tv}}, $tk);
		$linev{$tv} ++;
	}
	my $o_idx = 0; 
	print STDOUT join("\t", qw/ClusterID ElementNum Elements/)."\n"; 
	for my $tv (sort { $linev{$b}<=>$linev{$a} || $a<=>$b } keys %olines) {
		print STDOUT join("\t", $o_idx, $linev{$tv}, sort @{$olines{$tv}})."\n"; 
		$o_idx ++; 
	}
}

#cbind: bind files by columns. 
# Similar to linux-paste command. 
sub cbind {
	@ARGV > 1 or &stop( "[Err]There should be no less than two files for -cbind function.\n" ); 
	my @col_len = (-1) x scalar(@InFp); 
	while (!eof($InFp[0])) {
		$. % 1e6 == 1 and &tmsg( "[Msg]$. lines.\n" ); 
		my @line; 
		for (my $i=0; $i<@InFp; $i++) {
			my $tfh = $InFp[$i]; 
			&tmsg("[Rec]Dealing with [0-Idx $i] file [$ARGV[$i]]\n"); 
			if (eof($tfh)) {
				&stop("[Err][$i] file [$ARGV[$i]] ended before others!\n"); 
			}
			my $tl = <$tfh>; 
			chomp($tl); $tl =~ s/[\r\b\n]+$//; 
			if ($opts{chkColNum}) {
				my @cc = ($tl =~ m/($symbol)/og); 
				my $coln = $#cc+2; 
				if ($col_len[$i] == -1) {
					$col_len[$i] = $coln; 
				}elsif ($col_len[$i] != $coln) {
					&tmsg("[Wrn]Column number of file $ARGV[$i] changed at line ($col_len[$i] to $coln): '$tl'\n"); 
					&tmsg("[Wrn]I modify columns to $col_len[$i].\n"); 
					$tl = join("$symbol", (&splitL($symbol, $tl))[0..($col_len[$i]-1)]); 
				}
			}# End chkColNum 
			push(@line, $tl); 
		}# End for $i<@InFp 
		print STDOUT join("\t", @line)."\n"; 
	}# End while !eof($InFp[0])
	for (my $i=1; $i<@InFp; $i++) {
		!eof($InFp[$i]) and &tmsg("[Wrn][0-Idx $i] file [$ARGV[$i]] rests lines.\n"); 
	}
	&tmsg( "[Rec] -cbind Over.\n" ); 
}#End of sub cbind. 

#Labeltbl 2007-06-29
sub LabelTbl {
	my ($mark) = @_;
	$mark =~ s/^\s+//; $mark =~ s/\s+$//g;
	my @Marks = &splitL("::", $mark); 
	my @markT = ();
	for (my $i=0; $i<@Marks; $i++) {
		if ($Marks[$i] =~ /^([^.].*)\.{2}/ ) {
			$markT[$i]  = 'step';
			$Marks[$i]  = $1;
		}else{
			$markT[$i] = 'same';
		}
	}

	for my $fh (@InFp) {
		while (<$fh>) {
			# Make individual markers. 
			my @OutPre = ();
			foreach my $i (0..$#Marks) {
				push(@OutPre, $Marks[$i]);
				$markT[$i] eq 'same' or $Marks[$i]++;
			}
			# output. 
			print STDOUT join($symbol,@OutPre,$_); # Only in this function I choose "$symbol" to delimiting output columns. 
		}# End while $fh 
	}# End $fh in @InFp
}# End sub LabelTbl 

# UniqColLine
sub UniqColLineF{
	my %uniqLine;
	my @Cols = &parseCol($opts{UniqColLine});
	$opts{'topN'} //= 1; 
	for my $fh (@InFp) {
		while (<$fh>) {
			chomp; 
			s/[^\S$symbol]+$//; 
			my @temp = &splitL($symbol, $_); 
			$temp[-1] =~ s/[^\S$symbol]+$//; 
			my $key = join($symbol,@temp[@Cols]);
			if (defined $uniqLine{$key}) {
				$uniqLine{$key} ++; 
				$uniqLine{$key} > $opts{'topN'} and next; 
				print STDOUT "$_\n"; 
			}else{
				$uniqLine{$key} = 1;
				print STDOUT "$_\n";
			}
		}# End while $fh 
	}# End $fh in @InFp
}# 2006-11-24 13:41


### "best_uniqCol:s","select_col:s","select_rule:s","best_rep:s",2006-11-17 13:24
sub best_uniqCol{
	my $rep = 0;
	defined $opts{'best_rep'} and $rep = 1;
	$rep == 1 and do { open(REP,'>',$opts{best_rep}) or &stop( "[Err]open file $opts{best_rep} failed!\n" ) };
	my @BestCol   = &parseCol($opts{best_uniqCol});
	my @SelctCol  = &parseCol($opts{select_col});
#	my @SelctRule = &parseCol($opts{select_rule});
	my @SelctRule = &splitL(',', $opts{'select_rule'}); 
	foreach my $test (@SelctRule) {
		$test =~ s/(^\s+|\s+$)//; 
		($test == 1 or $test == -1) or &stop( "[Err] select_rule mustbe 1/-1!\n" );
	}
	($#BestCol > -1 and $#SelctCol > -1 and $#SelctRule > -1) or &stop("[Err] $#BestCol > -1 and $#SelctCol > -1 and $#SelctRule > -1\n");
	$#SelctCol == $#SelctRule or &stop("[Err] SelctCol != SelctRule\n");
	my %lines;
	my %slctLines;
	my %repLines;
	my %repCount;
	my @Keys;
	for my $fh (@InFp) {
		READING:while (<$fh>) {
			chomp; 
			my @temp = &splitL($symbol, $_); 
			my $key = ($opts{best_disOrdCol}) ? join($symbol, (sort @temp[@BestCol])) : join($symbol,@temp[@BestCol]) ; # 2013-02-08 
			my $newSlct = join($symbol,@temp[@SelctCol]);
			if (defined $slctLines{$key}) {
				my @OriginSelct = &splitL($symbol, $slctLines{$key}); 
				if ( $newSlct eq $slctLines{$key} ) {
					if ( $rep == 1 ) {
						$repLines{$key} .= "\n$_"; 
						$repCount{$key} ++;
					}
					$lines{$key} .= "\n$_"; 
					next READING; 
				}
				no warnings; 
				COMP:for (my $i=0; $i<@SelctCol; $i++) {
					my ($slctCol,$slctRule,$originSlct) = ($SelctCol[$i],$SelctRule[$i],$OriginSelct[$i]);
					my $sym_cmp = ($temp[$slctCol] <=> $originSlct || $temp[$slctCol] cmp $originSlct); 
					$sym_cmp == 0 and next COMP; 
					if ($sym_cmp == $slctRule) {
						$slctLines{$key} = $newSlct;
						$lines{$key} = $_;
						$repLines{$key} = $_;
						$repCount{$key} = 1;
					}
					last COMP; 
				}
			}else{
				$slctLines{$key} = $newSlct;
				$lines{$key} = $_;
				$repLines{$key} = $_;
				$repCount{$key} = 1;
				push(@Keys,$key);
			}
		} # end while
	}# End for $fh in @InFp
	foreach my $ttKey (@Keys) {
		print STDOUT "$lines{$ttKey}\n";
	}
	if ($rep == 1) {
		foreach my $ttKey (@Keys) {
			$repCount{$ttKey} > 1 and print REP "$repLines{$ttKey}\n";
		}
	}
	$rep == 1 and close REP;
#	exit 1;
}# 2007-11-02 09:34

#### reverse lines
sub reverse_lines{
	my @file;
	for my $fh (@InFp) {
		while (<$fh>) {
			chomp;
			push(@file,$_);
		}
	}# End $fh in @InFp 
	my @output=reverse(@file);
	@file=();
	foreach my $line (@output) {
		print STDOUT "$line\n";
	}
}# End sub reverse_lines 


#### skip head lines
sub skip{
	my $skip=$opts{skip}+1;
	for my $fh (@InFp) {
		while (<$fh>) {
			$. < $skip and next;
			print STDOUT $_;
		}
	}# End for $fh in @InFp 
}# end skip
#### combine files to STDOUT
sub combine{
	for my $fh (@InFp) {
		while (<$fh>) {
			print STDOUT $_;
		}
	}# End for $fh in @InFp 
}# end combine

#### select special cols to generate a new line
sub column{
	my @cols=&parseCol($opts{column});
	my $need_kick = 0; 
	$opts{'kick_col'} and $need_kick = 1; 
	for my $fh (@InFp) {
		if ( defined $pm ) {
			$opts{'cpuN'} > 1 or &stopErr("[Err] Check 01 in column()\n"); 
			my $wrk_dir = &fileSunhh::new_tmp_dir(); 
			mkdir($wrk_dir) or &stopErr("[Err] Failed to create dir [$wrk_dir]\n"); 
			my @sub_fn = &fileSunhh::dvd_file( $fh, $opts{'cpuN'}, 'keep_order' => 1, 'with_header' => 0, 'sub_pref' => "$wrk_dir/sub_", 'tmpFile' => "$wrk_dir/base_0" ); 
			for my $sfn (@sub_fn) {
				my $pid = $pm->start and next; 
				open F,'<',"$sfn" or &stopErr("[Err] Failed to open file [$sfn]\n"); 
				open O,'>',"$sfn.o" or &stopErr("[Err] Failed to open file [$sfn.o]\n"); 
				while (<F>) {
					chomp; s/[^\S$symbol]+$//; 
					my @temp = &splitL($symbol, $_); 
					if ($need_kick == 1) {
						my %t1 = map { $_ => 1 } @cols; 
						@cols = (); 
						for (my $i=0; $i<@temp; $i++) {
							defined $t1{$i} and next; 
							push(@cols, $i); 
						}
						$need_kick = 0; 
					}
					print O join("\t",@temp[@cols])."\n"; 
				}
				close O; 
				close F; 
				$pm->finish; 
			}
			$pm->wait_all_children; 
			for my $sfn ( @sub_fn ) {
				open F,'<',"$sfn.o" or &stopErr("[Err] Failed to open file [$sfn.o]\n"); 
				while (<F>) {
					print STDOUT $_; 
				}
				close F; 
			}
			&fileSunhh::_rmtree( $wrk_dir ); 
		} else {
			while (<$fh>) {
				chomp; s/[^\S$symbol]+$//; 
				my @temp = &splitL($symbol, $_); 
				if ($need_kick == 1) {
					my %t1 = map { $_ => 1 } @cols; 
					@cols = (); 
					for (my $i=0; $i<@temp; $i++) {
						defined $t1{$i} and next; 
						push(@cols, $i); 
					}
					$need_kick = 0; 
				}
				print STDOUT join("\t",@temp[@cols])."\n";
			}
		}
	}# End for $fh in @InFp 
}# end column


#### max_col/min_col give columns to search.
sub extreme{
	my @cols;
	my $mode;
	if ($opts{max_col} ne '') {
		@cols = &parseCol($opts{max_col});
		$mode=1;			# 1 for max mode;
	}else{
		@cols = &parseCol($opts{min_col});
		$mode=0;			# 0 for min mode;
	}
	my @file;
	for my $fh (@InFp) {
		while (<$fh>) {
			chomp;
			push(@file,$_);
		}
	}
	foreach my $col (@cols) {
		$#file>0 or last; # Should not be ">=". 
		my @rest;
		my @temp_extrm = &splitL($symbol, $file[0]); 
		my $extreme=$temp_extrm[$col];
		foreach my $line (@file) {
			my @temp_line = &splitL($symbol, $line); 
			if ($mode) {
				($extreme<$temp_line[$col])?($extreme=$temp_line[$col]):();
			}else{
				($extreme>$temp_line[$col])?($extreme=$temp_line[$col]):();
			}
		}
		foreach my $line (@file) {
			my @temp_line = &splitL($symbol, $line); 
			($extreme==$temp_line[$col])?(push(@rest,$line)):();
		}
		@file=@rest;
	}
	foreach my $output (@file) {
		print STDOUT "$output\n";
	}
	@file=();
}# end extreme

sub col_sort {
	no strict; 
	my @tempa = &splitL($symbol, $a); 
	my @tempb = &splitL($symbol, $b); 
	no warnings; 
	foreach my $rul (@{$opts{'col_sort_ruls'}}) {
		my $result = 0; 
		if ( &is_digital($tempa[$rul->[0]]) and &is_digital($tempb[$rul->[0]]) ) {
			$result = $tempa[$rul->[0]] <=> $tempb[$rul->[0]] || $tempa[$rul->[0]] cmp $tempb[$rul->[0]]; 
		} elsif ( &is_digital($tempa[$rul->[0]]) ) {
			$result = -1; 
		} elsif ( &is_digital($tempb[$rul->[0]]) ) {
			$result = 1; 
		} else {
			$result = $tempa[$rul->[0]] cmp $tempb[$rul->[0]]; 
		}
		$result != 0 and return( $result * $rul->[1] ); 
	}
	0;
}#end col_sort

sub is_digital {
	if      ($_[0] =~ m!^[+-]?\d+(\.?\d*)?([eE][+-]?\d+)?$!) {
		return 1; 
	} elsif ($_[0] =~ m!^[+-]?(\.\d+)?([eE][+-]?\d+)?$!) {
		return 1; 
	}
	return 0; 
}# is_digital()

sub _prep_cSort_ruls {
	defined $opts{'col_sort_ruls'} and @{$opts{'col_sort_ruls'}} > 0 and return; 
	$opts{col_sort} //= 0; 
	$opts{'col_sort_cols'} = [ &parseCol($opts{col_sort}) ]; 
	$opts{'col_sort_ruls'} = []; 
	if (defined $opts{'col_sort_rule'}) {
		my @rr = map { 
		  $_ =~ s/\s//g; 
		  $_ eq '' and $_ = -1; 
		  ($_ eq '1' or $_ eq '-1') or die "[Err] col_sort_rule should be like ' -1,1,-1,...'\n"; 
		  $_ * -1; 
		} &splitL(',', $opts{'col_sort_rule'}); 
		for (my $i=0; $i<@{$opts{'col_sort_cols'}}; $i++) {
			$rr[$i] //= 1; 
			push(@{$opts{'col_sort_ruls'}}, [$opts{'col_sort_cols'}[$i], $rr[$i]]); 
		}
	} else {
		for ( my $i=0; $i<@{$opts{'col_sort_cols'}}; $i++ ) {
			push(@{$opts{'col_sort_ruls'}}, [$opts{'col_sort_cols'}[$i], 1]); 
		}
	}
	return; 
}# sub _prep_cSort_ruls() 



#### cols select.with the form "num,num,num".., 0 for first column.
sub uniq_rep{
	my @cols;
	my %line;
	my %rep_count;
	my @keys;
	if (&goodVar($opts{col_uniq})) {
		@cols = &parseCol($opts{col_uniq});
	}elsif (&goodVar($opts{col_reps})) {
		@cols = &parseCol($opts{col_reps});
	}elsif ( &goodVar($opts{col_repCount}) ) {
		@cols = &parseCol($opts{col_repCount});
	}else{
		&stop( "[Err]No Cols Found!\n" );
	}
	for my $fh (@InFp) {
		while (<$fh>) {
			chomp; 
			my @temp = &splitL($symbol, $_); 
			my $key = join("\t",@temp[@cols]);
			if (defined $line{$key}) {
				$line{$key} .= "\n$_";
				$rep_count{$key}++;
			}else{
				push(@keys,$key);
				$line{$key} = $_;
				$rep_count{$key} = 1;
			}
		}# End while $fh
	}# End for $fh in @InFp
	if (&goodVar($opts{col_uniq})) {
		foreach my $key (@keys) {
			$rep_count{$key} == 1 and print STDOUT "$line{$key}\n";
		}
	}elsif ( &goodVar($opts{col_reps}) ) {
		foreach my $key (@keys) {
			$rep_count{$key} > 1 and print STDOUT "$line{$key}\n";
		}
	}elsif ( &goodVar($opts{col_repCount}) ) {
		print STDOUT "Number\tKeyCols\n";
		foreach my $key (@keys) {
			print STDOUT "$rep_count{$key}\t$key\n";
		}
	}
}#end uniq_rep


### 
sub colStat{
	my @Data;
	my $col = $opts{col_stat};
	my $total = 0; 
	my $useful_count = 0; 
	for my $fh (@InFp) {
		while (<$fh>) {
			chomp; 
			my @temp = &splitL($symbol, $_); 
			chomp; chomp($temp[-1]); 
			(defined $temp[$col] and $temp[$col] =~ m/^[\d.\-eE]+$/) or next; 
			$useful_count ++; 
			push(@Data,$temp[$col]); 
			$total += $temp[$col]; 
		}
	}# End for $fh in @InFp
	$#Data == -1 and &stop( "[Err]No Data found.\n" );
	my @SortData = sort {$a<=>$b;} @Data;
	my ($min,$max,$mean) = ($SortData[0],$SortData[$#SortData],$total/($#SortData+1));
	my $median;
	if ($#SortData%2==0) {
		my $i = $#SortData/2;
		$median = $SortData[$i];
	}else{
		my $i1 = ($#SortData+1)/2;
		my $i2 = $i1-1;
		$median = ($SortData[$i1]+$SortData[$i2])/2;
	}# get median
	if ( $opts{'col_stat_asINS'} ) {
		use mathSunhh; 
		my $bh = &ins_calc(\@Data); 
		my @add_tk = qw/interval_mean interval_median interval_stdev interval_low interval_high limit_low limit_high/; 
		my @add_tv = @{$bh}{@add_tk}; 
		print STDOUT join("\t", qw/SUM MEAN MEDIAN MIN MAX Count NoNull/, @add_tk)."\n"; 
		print STDOUT join("\t", $total, $mean, $median, $min, $max, scalar @SortData, $useful_count, @add_tv)."\n"; 
	} else {
		print STDOUT "SUM\tMEAN\tMEDIAN\tMIN\tMAX\tCount\tNoNull\n"; 
		print STDOUT "$total\t$mean\t$median\t$min\t$max\t",scalar @SortData,"\t$useful_count\n";
	}
}

######################################################################
##     Inner sub-routines. 
######################################################################
sub parseCol {
	my @cols = &splitL(',', $_[0]); 
	my @ncols; 
	for my $tc (@cols) {
		$tc =~ s/(^\s+|\s+$)//g; 
		if ($tc =~ m/^\d+$/) {
			push(@ncols, $tc); 
		} elsif ($tc =~ m/^(\d+)\-(\d+)$/) {
			my ($s, $e) = ($1, $2); 
			if ($s <= $e) {
				push(@ncols, ($s .. $e)); 
			}else{
				push(@ncols, reverse($e .. $s)); 
			}
		} else {
			&stop("[Err]Unparsable column tag for [$tc]\n");
		}
	}
	return (@ncols); 
}


sub goodVar {
	my ($a, $t) = @_; 
	defined $t and $t ne '' and $t = lc($t); 
	(!defined $t or $t eq '') and $t = 'nonull'; 
	if ($t eq 'nonull') {
		defined $a and $a ne '' and return 1; 
	} elsif ( $t eq 'defined' ) {
		defined $a and return 1; 
	}
	return 0; 
}# End goodVar

sub tmsg {
	my $tt = scalar( localtime() ); 
	print STDERR join('', "[$tt]", @_); 
}# End tmsg 

sub stop {
	my $tt = scalar( localtime() ); 
	print STDERR join('', "[$tt]", @_); 
	exit 1; 
}# End stop 

sub loc_tbl2hash {
	# ID \\t Start \\t End \\n
	my $fh = shift; 
	my %loc; 
	while (<$fh>) {
		chomp; 
		my @ta = &splitL($symbol, $_); 
		push(@{$loc{$ta[0]}}, [$ta[1] , $ta[2]]); 
	}
	for my $id (keys %loc) {
		my @new_blk; 
		for my $ar1 (sort {$a->[0] <=> $b->[0] || $a->[1] <=> $b->[1]} @{$loc{$id}}) {
			if (defined $new_blk[0]) {
				if ( $new_blk[-1][1]+1 >= $ar1->[0] ) {
					$ar1->[1] > $new_blk[-1][1] and $new_blk[-1][1] = $ar1->[1]; 
				} else {
					push(@new_blk, [$ar1->[0], $ar1->[1]]); 
				} 
			} else {
				push(@new_blk, [$ar1->[0], $ar1->[1]]); 
			}
		}
		$loc{$id} = \@new_blk; 
	}
	return \%loc; 
}# loc_tbl2hash() 


