#!/usr/bin/perl
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

use strict;

use Getopt::Long;
my %opts;
GetOptions(\%opts,
			"combine!","column:s",
			"max_col:s","min_col:s",
			"skip:i",
			"reverse!",
			"col_sort:s",
			"col_uniq:s","col_reps:s","col_repCount:s",		# uniq/repeat columns pick.
			"col_stat:i",
			"symbol:s",
			"best_uniqCol:s","select_col:s","select_rule:s","best_rep:s","best_disOrdCol!", 
			"UniqColLine:s",
			"label_mark:s",
			"cbind!", "chkColNum!", 
			"cluster_lines!", "cluster_cols:s", 
			"trimEndReturn!", 
			"fillNull:s", 
			"transpose!", "skip_null_line!", "beMatrix!", "fill_new:s",      # transpose a matrix. 
			"help!");
my $info=<<INFO;
##################################################
command:perl $0 <STDIN|parameters>
# 2013-08-19 

	-help           help infomation;
	-combine        combine files;
	-column<int>    cols "num1,num2,num3..." will be picked out and joined to a new line;
	-max_col        Similar to -column, compare by order;
	-min_col        Similar to -column, compare by order;
	-skip<int>      Skip first <int> lines.Usually for skipping head lines;
	-col_sort       sort by col.

	-col_uniq       pick lines with its columns uniq
	-col_reps       pick lines with its columns repeat
	-col_repCount   Show col repeat times.

	-col_stat<int>  Get Sum, Mean, Median, Min, Max of the Column;
	-reverse        reverse the file lines order.
	
	-best_uniqCol   1st column is 0; col_1[,col_2...]
	-select_col     1st column is 0; col_1[,col_2...]
	-select_rule    1/-1[,1/-1...]; 1 for Max, -1 for Min;
	-best_rep       out a best uniq cols' repeat records if given;
	-best_disOrdCol Do not care the order of uniqCols, so "A\\tB" line will be compared with "B\\tA" line. 
	
	-UniqColLine    output a line with its col=UniqColLine one time.<*,*,*...>

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

	-symbol         Defining the symbol to divide data.Default is "\\t";
##################################################
INFO

sub usage {
	print STDOUT "$info"; 
	print STDERR "Error in $0\n"; 
	exit(1); 
}
if ($opts{help} or (-t and !@ARGV)) {
	&usage; 
}

my $symbol = "\t";
$opts{symbol} ne '' and $symbol = $opts{symbol};


&reverse_lines() if($opts{reverse});
&skip() if ($opts{skip}>0);
&combine() if ($opts{combine});
&column() if ($opts{column} ne '');
&extreme() if ($opts{max_col} ne '' || $opts{min_col} ne '');
&colStat() if ($opts{col_stat} ne '');
if ($opts{col_sort} ne '') {
	my @file;
	while (<>) {
		chomp;
		push(@file,$_);
	}continue{
		if (eof) {
			my @output=sort col_sort @file;
			foreach my $line (@output) {
				print STDOUT "$line\n";
			}
		}
	}
}# end if

&uniq_rep() if ($opts{col_uniq} ne '' || $opts{col_reps} ne '' || $opts{col_repCount} ne '');

&best_uniqCol() if ($opts{best_uniqCol} ne '');

&UniqColLineF() if ($opts{UniqColLine} ne '');

&LabelTbl($opts{label_mark}) if (defined $opts{label_mark});

&cbind() if ($opts{cbind}); 

&cluster_group() if ($opts{cluster_lines}); 

&trimEndReturn() if ( $opts{trimEndReturn} ); 

&fillNull() if ( defined $opts{fillNull} ); 

&t() if ( $opts{transpose} ); 

######################################################################
######################################################################

# Transpose a matrix . Like the t() function in R language. 
sub t {
	my @result_lines; 
	my $line_idx = -1; 
	while (<>) {
		chomp; 
		$opts{skip_null_line} and /^\s*$/ and next; # Here empty line should not be skipped until demanded. 
		my @ta = split(/$symbol/o, $_); 
		$line_idx ++; 
		for (my $i=0; $i<@ta; $i++) {
			$result_lines[$i][$line_idx] = $ta[$i]; 
		}
	}
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
	while (<>) {
		while (s/$symbol$symbol/$symbol$tc$symbol/og) { 1; } 
		s/^$symbol/$tc$symbol/o; 
		s/$symbol$/$symbol$tc/o; 
		print; 
	}
}# End sub fillNull

# trim end returns from windows. 
sub trimEndReturn {
	while (<>) {
		s/[\b\r\n]+$//; 
		print "$_\n"; 
	}
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
	while (<>) {
		$. % 1000 == 1 and warn "$. lines.\n"; 
		chomp; /^\s*$/ and next; 
		my @ta=split(/$symbol/o, $_); 

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
sub cbind {
	@ARGV > 1 or die "There should be no less than two files for -cbind function.\n"; 
	my @fileH; 
	my @col_len; 
	for (@ARGV) {
		my $tfh; 
		if ($_ =~ m/\.gz$/) {
			open ($tfh, '-|', "gzip -cd $_") or die; 
		}elsif ($_ =~ m/\.bz2$/) {
			open ($tfh, '-|', "bzip2 -cd $_") or die; 
		}else{
			open ($tfh, '<', "$_") or die; 
		}
		push(@fileH, $tfh); 
		push(@col_len, -1); 
		warn "[Record]File handle [$tfh] for [0-Idx $#fileH] file [$_] " . scalar(localtime()) . "\n";
	}
	while (!eof($fileH[0])) {
		$. % 1000000 == 1 and warn "[Msg]$. lines. " . scalar(localtime()) . "\n"; 
		my @line; 
		for (my $i=0; $i<@fileH; $i++) {
			my $tfh = $fileH[$i]; 
			if (eof($tfh)) {
				die "[Error][$i] file [$ARGV[$i]] (file handle [$tfh]) ended before others!\n"; 
			}
			my $tl = <$tfh>; 
			chomp($tl); $tl =~ s/[\b\n]+$//; 
			if ($opts{chkColNum}) {
				my @cc = ($tl =~ m/($symbol)/og); 
				my $coln = $#cc+2; 
				if ($col_len[$i] == -1) {
					$col_len[$i] = $coln; 
				}elsif ($col_len[$i] != $coln) {
					warn "[Warn]Column number of file $ARGV[$i] changed at line ($col_len[$i] to $coln): '$tl'\n[Warn]I modify columns to $col_len[$i]. " . scalar(localtime()) . "\n"; 
					$tl = join("$symbol", (split(/$symbol/o, $tl))[0..($col_len[$i]-1)]); 
				}
			}


			push(@line, $tl); 
		}
		print STDOUT join("\t", @line)."\n"; 
	}
	for (my $i=1; $i<@fileH; $i++) {
		!eof($fileH[$i]) and warn "[Warn][0-Idx $i] file [$ARGV[$i]] (file handle [$fileH[$i]]) rests lines. " . scalar(localtime()) . "\n"; 
	}
	warn "[Record]-cbind Over. " . scalar(localtime()) . "\n"; 
}#End of sub cbind. 

#Labeltbl 2007-06-29
sub LabelTbl {
  my ($mark) = @_;
  $mark =~ s/^\s+//; $mark =~ s/\s+$//g;
  my @Marks = split(/::/,$mark);
  my @markT = ();
  for (my $i=0; $i<@Marks; $i++) {
    if ($Marks[$i] =~ /^([^.].*)\.{2}/ ) {
      $markT[$i]  = 'step';
      $Marks[$i]  = $1;
    }else{
      $markT[$i] = 'same';
    }
  }
  
  while (<>) {
    my @OutPre = ();
    foreach my $i (0..$#Marks) {
      push(@OutPre, $Marks[$i]);
      $markT[$i] eq 'same' or $Marks[$i]++;
    }
    print STDOUT join($symbol,@OutPre,$_); # Only in this function I choose "$symbol" to delimiting output columns. 
  }
}








# UniqColLine
sub UniqColLineF{
	my %uniqLine;
	my @Cols = &parseCol($opts{UniqColLine});
	while (<>) {
		chomp; s/[^\S$symbol]+$//; 
		my @temp = split(/$symbol/o,$_);
		my $key = join($symbol,@temp[@Cols]);
		if (defined $uniqLine{$key}) {
		}else{
			$uniqLine{$key} = 1;
			print STDOUT "$_\n";
		}
	}
}# 2006-11-24 13:41


### "best_uniqCol:s","select_col:s","select_rule:s","best_rep:s",2006-11-17 13:24
sub best_uniqCol{
	my $rep = 0;
	$opts{best_rep} ne '' and $rep = 1;
	$rep == 1 and (open(REP,'>',$opts{best_rep}) or die "open file $opts{best_rep} failed!\n$info");
	my @BestCol   = &parseCol($opts{best_uniqCol});
	my @SelctCol  = &parseCol($opts{select_col});
#	my @SelctRule = &parseCol($opts{select_rule});
	my @SelctRule = split(/,/, $opts{select_rule});
	foreach my $test (@SelctRule) {
		$test =~ s/(^\s+|\s+$)//; 
		($test == 1 or $test == -1) or die "select_rule mustbe 1/-1!\n";
	}
	($#BestCol > -1 and $#SelctCol > -1 and $#SelctRule > -1) or die $info;
	$#SelctCol == $#SelctRule or die $info;
	my %lines;
	my %slctLines;
	my %repLines;
	my %repCount;
	my @Keys;
	READING:while (<>) {
		chomp;
		my @temp = split(/$symbol/o,$_);
		my $key = ($opts{best_disOrdCol}) ? join($symbol, (sort @temp[@BestCol])) : join($symbol,@temp[@BestCol]) ; # 2013-02-08 
		my $newSlct = join($symbol,@temp[@SelctCol]);
		if (defined $slctLines{$key}) {
			my @OriginSelct = split(/$symbol/o,$slctLines{$key});
			if ($rep == 1 and $newSlct eq $slctLines{$key}) {
				$repLines{$key} .= "\n$_";
				$repCount{$key} ++;
				next READING;
			} # 2007-1-17 16:55
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
	while (<>) {
		chomp;
		push(@file,$_);
	}continue{
		if (eof) {
			close ARGV;
			my @output=reverse(@file);
			@file=();
			foreach my $line (@output) {
				print STDOUT "$line\n";
			}
		}
	}
}


#### skip head lines
sub skip{
	my $skip=$opts{skip}+1;
	while (<>) {
		$.<$skip and next;
		print STDOUT $_;
	}continue{
		if (eof) {
			close ARGV;
		}
	}
}# end skip
#### combine files to STDOUT
sub combine{
	while (<>) {
		print STDOUT $_;
	}
}# end combine

sub parseCol {
	my @cols = split(/,/, $_[0]); 
	my @ncols; 
	for my $tc (@cols) {
		$tc =~ s/(^\s+|\s+$)//; 
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
			die "Unparsable column tag for [$tc]\n";
		}
	}
	return (@ncols); 
}

#### select special cols to generate a new line
sub column{
	my @cols=&parseCol($opts{column});
	while (<>) {
		chomp; s/[^\S$symbol]+$//; 
		my @temp=split(/$symbol/o,$_);
		print STDOUT join("\t",@temp[@cols])."\n";
	}
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
	while (<>) {
		chomp;
		push(@file,$_);
	}continue{
		if (eof) {
			close ARGV;
			foreach my $col (@cols) {
				$#file>0 or last; # Should not be ">=". 
				my @rest;
				my @temp_extrm=split(/$symbol/o,$file[0]);
				my $extreme=$temp_extrm[$col];
				foreach my $line (@file) {
					my @temp_line=split(/$symbol/o,$line);
					if ($mode) {
						($extreme<$temp_line[$col])?($extreme=$temp_line[$col]):();
					}else{
						($extreme>$temp_line[$col])?($extreme=$temp_line[$col]):();
					}
				}
				foreach my $line (@file) {
					my @temp_line=split(/$symbol/o,$line);
					($extreme==$temp_line[$col])?(push(@rest,$line)):();
				}
				@file=@rest;
			}
			foreach my $output (@file) {
				print STDOUT "$output\n";
			}
			@file=();
		}
	}
}# end extreme

sub col_sort{
	my @tempa=split(/$symbol/o,$a);
	my @tempb=split(/$symbol/o,$b);
	my @cols = &parseCol($opts{col_sort});
	foreach my $col (@cols) {
		if (my $result=(($tempa[$col]<=>$tempb[$col]) || ($tempa[$col] cmp $tempb[$col]))) {
			return $result;
		}
	}
	0;
}#end col_sort

#### cols select.with the form "num,num,num".., 0 for first column.
sub uniq_rep{
	my @cols;
	my %line;
	my %rep_count;
	my @keys;
	if ($opts{col_uniq} ne '') {
		@cols = &parseCol($opts{col_uniq});
	}elsif ($opts{col_reps} ne '') {
		@cols = &parseCol($opts{col_reps});
	}elsif ($opts{col_repCount} ne '') {
		@cols = &parseCol($opts{col_repCount});
	}else{
		die "No Cols Found!\n";
	}
	while (<>) {
		chomp;
		my @temp=split(/$symbol/o,$_);
		my $key = join("\t",@temp[@cols]);
		if (defined $line{$key}) {
			$line{$key} .= "\n$_";
			$rep_count{$key}++;
		}else{
			push(@keys,$key);
			$line{$key} = $_;
			$rep_count{$key} = 1;
		}
	}
	if ($opts{col_uniq} ne '') {
		foreach my $key (@keys) {
			$rep_count{$key} == 1 and print STDOUT "$line{$key}\n";
		}
	}elsif ($opts{col_reps} ne '') {
		foreach my $key (@keys) {
			$rep_count{$key} > 1 and print STDOUT "$line{$key}\n";
		}
	}elsif ($opts{col_repCount} ne '') {
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
	while (<>) {
		chomp;
		my @temp=split(/$symbol/o,$_); 
		push(@Data,$temp[$col]); 
		$total += $temp[$col]; 
		defined $temp[$col] and $temp[$col] ne '' and $useful_count++; 
	}
	$#Data == -1 and die "No Data found.\n";
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
	print STDOUT "SUM\tMEAN\tMEDIAN\tMIN\tMAX\tCount\tNoNull\n"; 
	print STDOUT "$total\t$mean\t$median\t$min\t$max\t",scalar @SortData,"\t$useful_count\n";
}



