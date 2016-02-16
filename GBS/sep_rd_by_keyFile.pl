#!/usr/bin/perl 
# 2016-02-16 I don't know why, but FileHandle package doen't work well (skipping reads) in this script. So I temporarily skip this method. 
#   ######## The current 'write2file()' solution needs 50% more time than FileHandle method. 
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use Getopt::Long; 
use FileHandle; 
# $| = 1; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"inKey:s", # key.txt 
	"inFqDir:s", # /path/to/fastq/ : . 
	"outFqDir:s", # ./out_sep ; 
	"keep_index!", 
); 

my $usage = <<HH; 

perl $0 

-help 

-inKey       [key.txt] Required. 
-inFqDir     ['./'] Default use current directory. 
-cutEnd      ['CAGC,CTGC'] ',' means there are more than one possibilities. 
-outFqDir    ['./out_sep']

-keep_index  [Boolean] Default False. Will keep index barcodes in the beginning of reads. 

HH

$opts{'help'} and &LogInforSunhh::usage($usage); 

defined $opts{'inKey'} or &LogInforSunhh::usage($usage); 

$opts{'cutEnd'} //= 'CAGC,CTGC'; 
my @ends; 
for (split(/,/, $opts{'cutEnd'})) {
	s/\s//g; 
	push(@ends, $_); 
}

$opts{'outFqDir'} //= './out_sep'; 
-e $opts{'outFqDir'} or mkdir($opts{'outFqDir'}, 0755); 

my %key_hash = %{ &readInKeyFile( $opts{'inKey'} ) }; 

$opts{'inFqDir'} //= '.'; 
opendir ID,"$opts{'inFqDir'}" or die ; 
my @inFile; 
while ( readdir ID ) {
	chomp; 
	# C1AADACXX_rmAdp3low_3_fastq.txt 
	m/^([^_\s]+)(?:_[^_\s]+)?_(\d+)_fastq(?:.txt)?(?:\.gz)?$/ or do { &tsmsg("[Msg] Skip file [$_]\n"); next; }; 
	my $kk = "$1\t$2"; 
	defined $key_hash{'lane2info'}{$kk} or do { &tsmsg("[Msg] Skip file [$_] not in key file\n"); next; }; 
	push(@inFile, [ "$opts{'inFqDir'}/$_", $1, $2 ]); # [ fastq_file, Flowcell, Lane ] 
}
closedir ID; 

unlink( "$opts{'outFqDir'}/rabbish.fq" ); 
for my $tr1 (@inFile) {
	my ($fq, $flowcell, $lane) = @$tr1; 
	my $kk = "$flowcell\t$lane"; 
	my %idx2fname; 
	my %idx2ofileH; 
	my %idx2startP; 
	for my $tr2 ( @{$key_hash{'lane2info'}{$kk}} ) {
		my ( $fname, $barc, $sample ) = @$tr2; 
		my $ofile = $fname; $ofile =~ s![^\w\d]!_!g; 
		$ofile = "$ofile.fq"; 
		for my $te ( @ends ) {
			my $idx = "${barc}$te"; 
			$idx2fname{$idx} = $fname; 
			my $ofile_D = "$opts{'outFqDir'}/$ofile"; 
			unlink($ofile_D); 

			&fileSunhh::write2file($ofile_D, '', '>'); 
			$idx2ofileH{$idx} = $ofile_D; 
			$idx2startP{$idx} = length( $barc ) ; 
			# There should be another way, but I failed to figure out its problem. 
			### Ref : http://perl.plover.com/FAQs/Buffering.html
			### Ref : http://perldoc.perl.org/FileHandle.html
			# $idx2ofileH{$idx} = FileHandle->new("$ofile_D", '>');
			# $idx2ofileH{$idx}->autoflush(1);
			# 
		}
	}
	my %idx_inf = %{ &info_idx([ keys %idx2ofileH ]) }; 
	defined $idx2ofileH{'rabbish'} and &stopErr("[Err] defined rabbish\n"); 

#	$idx2ofileH{'rabbish'} = &openFH("$opts{'outFqDir'}/rabbish.fq", '>>'); 
	$idx2ofileH{'rabbish'} = "$opts{'outFqDir'}/rabbish.fq"; 

	# Ref : http://perldoc.perl.org/FileHandle.html 
	# $idx2ofileH{'rabbish'} = FileHandle->new("$opts{'outFqDir'}/rabbish.fq", '>');
	# $idx2ofileH{'rabbish'}->autoflush(1); 
	#

	&tsmsg("[Msg] Idx_min=$idx_inf{'minL'} Idx_max=$idx_inf{'maxL'}\n"); 
	&tsmsg("[Msg] Separating raw file [$fq]\n"); 
	my $fqFh = &openFH( $fq, '<' ); 
	my ($l1, $l2, $l3, $l4, $sub_seq); 
	my $fix_idx = 'rabbish'; 
	my $sub_idx = ''; 
	RD: 
	while ($l1 = <$fqFh>) {
		$. % 4e7 == 1 and &tsmsg("[Msg] $. lines.\n"); 
		$l2 = <$fqFh>; 
		$l3 = <$fqFh>; $l4 = <$fqFh>; 
		$fix_idx = 'rabbish'; 
		for (my $i=$idx_inf{'maxL'}; $i>=$idx_inf{'minL'}; $i--) {
			$sub_idx = substr($l2, 0, $i); 
#			defined $idx2ofileH{ $sub_idx } and do { $fix_idxo = $sub_idx; $l2 = substr($l2, $idx2startP{$fix_idx}); $l4 = substr($l4, $idx2startP{$fix_idx}); last; }; 
			defined $idx2ofileH{ $sub_idx } and do { $fix_idx = $sub_idx; $fix_idx=$sub_idx; unless ($opts{'keep_index'}) { $l2 = substr($l2, $idx2startP{$fix_idx}); $l4 = substr($l4, $idx2startP{$fix_idx}); } last; }; 
		}

		&fileSunhh::write2file( $idx2ofileH{ $fix_idx }, "$l1$l2$l3$l4", '>>' ); 

		# Ref : http://perldoc.perl.org/FileHandle.html
		# print { $idx2ofileH{ $fix_idx } } "$l1$l2$l3$l4"; 
		#
	}
	close($fqFh); 
}

##### Sub ##### 
sub readInKeyFile {
	my ($file) = @_; 
	my $fh = &openFH($file, '<'); 
	my %back; 
	while (<$fh>) {
		chomp; s/[^\t\S]+$//; 
		my @ta = split(/\t/, $_); 
		$. == 1 and $ta[0] eq 'Flowcell' and next; 
		my $kk = "$ta[0]\t$ta[1]"; 
		my $full_name = join(':', @ta[3,0,1,7]); 
		push(@{$back{'lane2info'}{$kk}}, [$full_name, $ta[2], $ta[3]]); # [ DNASample:Flowcell:Lane:LibraryPrepID, Barcode, DNASample ] 
		$back{'lane_barc2fname'}{$kk}{$ta[2]} = $full_name; # {Flowcell \t Lane}{Barcode} => DNASample:Flowcell:Lane:LibraryPrepID 
		$back{'lane_barc2sample'}{$kk}{$ta[2]} = $ta[3]; # {Flowcell \t Lane}{Barcode} => DNASample 
	}
	close ($fh); 
	return(\%back); 
}# readInKeyFile () 

sub info_idx {
	my ($seq) = @_; 
	my %back; 
	for my $ts (@$seq) {
		my $ll = length($ts); 
		$back{'minL'} //= $ll; 
		$back{'maxL'} //= $ll; 
		$back{'minL'} > $ll and $back{'minL'} = $ll; 
		$back{'maxL'} < $ll and $back{'maxL'} = $ll; 
	}
	return(\%back); 
}# info_idx


