#!/usr/bin/perl

# file 1 for key, file 2 to srch.
# 2007-3-21 16:52 add a -sign
# 2008/05/08 change -f2 to STDIN. 
# 2009/01/09 add a 'link' para; 
# 2012/11/26 Flexible Column assigning. 10-12,35-14
# 2013/01/04 Fix a bug for typo "Col1" to "col1". 
# 2013-08-15 add a parameter 'fill' for joining absent lines. 
# 2013-08-21 Add some recording information because I find this script runs too slow than I expected. 

use strict;
use Getopt::Long;
my %opts;
GetOptions(\%opts,
			"keyC1:s","keyC2:s",
			"Col1:s","Col2:s",
			"f1:s",
			"out:s",
			"add!",
			"sign:s",
			"link!", # used for join two files. 
			"fill:s", 
			"help!");
my $info = <<INFO;
############################################################
perl $0 srchFileToAddTo -keyC1 0 -Col1 0,1,.. -keyC2 0 -Col2 0,1,.. -f1 KeyFile 

-add        if given, -Col2 will no use.
-link       if given, -Col1 will no use.
-sign       \'yes,no\' used to tag existence for Line2 in Line1, instead of combine columns. 
-out        Output file name, default to STDOUT. 
-fill       String used to fill added columns in absent lines. 
############################################################
INFO

$opts{help} and die "$info";
!@ARGV and -t and die "$info"; 

sub parseCol {
  my $cc = shift(@_); 
  my @return_col; 
  for my $tc1 (split(/,/, $cc)) {
    if ($tc1 =~ m/^\d+$/) {
      push(@return_col, $tc1); 
    } elsif ($tc1 =~ m/^(\d+)\-(\d+)$/) {
      if ($1 <= $2) {
        push(@return_col, ($1 .. $2)); 
      }else{
        push(@return_col, reverse($2 .. $1)); 
      }
    } else {
      die "Unknown Column=[$tc1] in [$cc]\n"; 
    }
  }
  return @return_col; 
}#End sub parseCol

my @Col1 = &parseCol($opts{Col1});
my @Col2 = &parseCol($opts{Col2});
my @KC1 = (0); 
my @KC2 = (0); 
defined $opts{keyC1} and @KC1 = &parseCol($opts{keyC1}); 
defined $opts{keyC2} and @KC2 = &parseCol($opts{keyC2}); 
#my $kc1 = $opts{keyC1};
#my $kc2 = $opts{keyC2};
my $file1 = $opts{f1};
my $fh_out = (defined $opts{out}) ? &openFHwrite($opts{out}) : \*STDOUT;
my ($yS,$nS);
if (defined $opts{sign}) {
	($yS,$nS) = split(/,/,$opts{sign});
	$yS =~ s/^\s+//;
}


open (F1,'<',$file1) or die "f1?\n$info";
print STDERR "[Stat] Begin to parse -f1 [$file1]. " . scalar(localtime()) . "\n"; 
my %need;		# 
while (<F1>) {
	chomp;
	my @temp = split(/\t/,$_);
	my $tk = join("\t", @temp[@KC1]); 
	if ($opts{link}) {
		my @tt; 
		F1_COL:
		for my $id1 (0..$#temp) {
			for my $id2 (@KC1) {
				$id1 == $id2 and next F1_COL;
			}
			push(@tt, $id1);
		}# end of F1_COL;
		$need{$tk} = join("\t", @temp[@tt]);
	}else{
		$need{$tk} = join("\t",@temp[@Col1]);
	}
}
close F1;
print STDERR "[Stat] -f1 $file1 parsed. " . scalar(localtime()) . "\n"; 

while (<>) {
	$. % 100000 == 1 and print STDERR "[Stat] $. line adding. " . scalar(localtime()) . "\n"; 
	chomp;
	my @temp = split(/\t/,$_);
	my ($line,$add);
	my $tk = join("\t", @temp[@KC2]); 
	if (defined $need{$tk}) {
		$add = $need{$tk};
	}else{
		if (defined $opts{fill}) {
			$add = join("\t", ($opts{fill}) x ($#Col1+1) ); 
		}else{
			my @tt;
			$#tt = $#Col1;
			$add = join("\t",@tt);
		}
	}
	if ($opts{add}) {
		$line = $_;
	}else{
		$line = join("\t",@temp[@Col2]);
	}
	if (defined $opts{sign}) {
		$add = (defined $need{$tk})?$yS:$nS;
	}
	print $fh_out "$line\t$add\n";
}
#close OUT;
print STDERR "[Stat] All finished. " . scalar(localtime()) . "\n"; 

sub openFHwrite {
	my $file = shift;
	my $fh;
	open $fh, '>', "$file" or die;
	return $fh;
}



