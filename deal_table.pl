#!/bin/env perl
# Author email : hs738@cornell.edu or biosunhh@gmail.com

# 2006-11-17 13:18 准备加入一个best_uniq函数，用来处理*.psl格式文件，尚未完成；
# 2006-11-24 13:33 加入一个小函数,UniqColLine,finished.
# 2006-11-24 15:16 测试best_uniq
# 2007-1-17 17:16 略微测试了best_uniqCol, 完成best_rep参数;
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

use strict;
use warnings; 
use fileSunhh; 
use LogInforSunhh; 

use Getopt::Long;
my %opts;
GetOptions(\%opts,
	"combine!","column:s",
	"colByTbl:s", "colByTbl_also:s", "colByTbl_cut:i", 
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
  -colByTbl<Str>  [index_file]. The 1st column of index_file is used to match the first line of input.table. Output matching columns from input.table. 
   -colByTbl_also [Str] Cols "num1,num2,num3" will be output before checking the matching between input.table and index_file. 
   -colByTbl_cut  [CPU_number] Use 'cut' and 'paste' command in Linux to try to shorten the processing time. 
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

  -symbol         Defining the symbol to divide data.Default is "\\t";
  -dR2dN          [Boolean] Change \\r to \\n in files. 

  -log_ln         [0]
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


&reverse_lines() if($opts{reverse});
&skip() if (&goodVar($opts{skip}) and $opts{skip}>0);
&combine() if ($opts{combine});
&column() if ( &goodVar($opts{column}) );
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

for (@InFp) {
	close ($_);
}


######################################################################
## sub-routines for functions. 
######################################################################

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
	use Parallel::ForkManager; 
	my @col_also = (); 
	defined $opts{'colByTbl_also'} and @col_also = &parseCol($opts{'colByTbl_also'}); 
	$opts{'colByTbl_cut'} //= 0; 
	if ( $opts{'colByTbl_cut'} > 0 ) {
		use Parallel::ForkManager; 
		$glob{'pm'} = new Parallel::ForkManager($opts{'colByTbl_cut'}); 
	}
	my $idxFh = &openFH($opts{'colByTbl'}, '<'); 
	my (%needColID, $cnt); 
	$cnt = 0; 
	while (<$idxFh>) {
		chomp; 
		my @ta = &splitL($symbol, $_); 
		defined $needColID{$ta[0]} or do { $needColID{$ta[0]} = $cnt; $cnt ++; }; 
	}
	close($idxFh); 

	if ( $opts{'colByTbl_cut'} > 0 ) {
		my $wrk_dir = &fileSunhh::new_tmp_dir(); 
		mkdir($wrk_dir) or &stopErr("[Err] Failed to create tmp_dir [$wrk_dir]\n"); 
		my $argv_i = -1; 
		for my $fh ( @InFp ) {
			$argv_i ++; 
			my @cur_ColNs = @col_also; 
			my $has_empty = 0; 
			{
				my $first_line = <$fh>; chomp($first_line); 
				my @ta = &splitL($symbol, $first_line); 
				my %cur_ID2ColN; 
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
				&tsmsg("[Msg] Copying raw file\n"); 
				my $base0_fh = &openFH("$wrk_dir/base0", '>'); 
				print {$base0_fh} "$first_line\n"; 
				while (<$fh>) {
					print {$base0_fh} $_; 
				}
				close($base0_fh); 
				&tsmsg("[Msg] Finish copying file\n"); 
			}
			for (my $i=0; $i<@cur_ColNs; $i++) {
				my $pid = $glob{'pm'}->start and next; 
				if ( $cur_ColNs[$i] eq 'NA' ) {
					unless ( $has_empty ) {
						my $iC = $#cur_ColNs+2; 
						&exeCmd_1cmd("cut -f $iC $wrk_dir/base0 > $wrk_dir/empty"); 
						$has_empty = 1; 
					}
					&exeCmd_1cmd("ln -s empty $wrk_dir/s$i"); 
				} else {
					my $iC = $cur_ColNs[$i]+1; 
					&exeCmd_1cmd("cut -f $iC $wrk_dir/base0 > $wrk_dir/s$i"); 
				}
				$glob{'pm'}->finish; 
			}
			$glob{'pm'}->wait_all_children; 
			my $cmdLn = join(" ", 'paste', map { "$wrk_dir/s$_" } (0 .. $#cur_ColNs)); 
			open O, '-|', "$cmdLn" or &stopErr("[Err] Failed to exe paste [$cmdLn]\n"); 
			while (<O>) {
				print STDOUT $_; 
			}
			close O; 
		}
		&fileSunhh::_rmtree($wrk_dir); 
	} else {
		for my $fh ( @InFp ) {
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
				print STDOUT join("\t", map { ( $_ eq 'NA' ) ? '' : $ta[$_] ; } @cur_ColNs)."\n"; 
			}
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
		while (<$fh>) {
			chomp; 
			if ( $. <= $opts{chID_skipH} ) {
				print STDOUT "$_\n"; 
				next; 
			} else {
				$cur_n ++; 
			}
			if ( m!^\s*(#|$)! ) {
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
			while (<$fh>) {
				$opts{'log_ln'} > 0 and &fileSunhh::log_section( $. , \%tmp_cnt ) and &tsmsg("[Msg] Reading [$.] line in kSrch_src\n"); 
				if ($opts{kSrch_drop}) {
					$idx_key{$_} or print; 
				}else{
					$idx_key{$_} and print; 
				}
			}#End while 
		}else{
			# Using the selected columns as key. 
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
		while (<$fh>) {
			while (s/$symbol$symbol/$symbol$tc$symbol/og) { 1; } 
			s/^$symbol/$tc$symbol/o; 
			s/$symbol$/$symbol$tc/o; 
			print; 
		}# End while 
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
	for my $fh (@InFp) {
		while (<$fh>) {
			chomp; s/[^\S$symbol]+$//; 
			my @temp = &splitL($symbol, $_); 
			print STDOUT join("\t",@temp[@cols])."\n";
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
	&_prep_cSort_ruls(); 
	my @tempa = &splitL($symbol, $a); 
	my @tempb = &splitL($symbol, $b); 
	no warnings; 
	foreach my $rul (@{$opts{'col_sort_ruls'}}) {
		if ( my $result=( ($tempa[$rul->[0]]<=>$tempb[$rul->[0]]) || ($tempa[$rul->[0]] cmp $tempb[$rul->[0]]) ) * $rul->[1] ) {
			return $result;
		}
	}
	0;
}#end col_sort

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
			(defined $temp[$col] and $temp[$col] =~ m/^[\d.\-e]+$/) or next; 
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


