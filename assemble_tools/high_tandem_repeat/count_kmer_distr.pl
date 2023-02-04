#!/usr/bin/perl
### Do not consider complementary. Different strand means different k-mers.
use strict;
use warnings;
use fastaSunhh;
my $fa_obj = fastaSunhh->new();
use fileSunhh;
use LogInforSunhh;
use Getopt::Long;
my %opts;
GetOptions(\%opts,
  "window_len:i", # 10e3
  "window_step:i", # same as window_len
  "kmer_len:i",   # 3
  "mer_store:s",  # 
  "canonical!",
  "out_dir:s",
  "verbose!",
  "help!",
);
sub usage {
  print STDOUT <<HH;
################################################################################
perl $0 in_sampleID.chr.fa -out_dir out_sampleID_dir/

-window_len         [10e3]
-window_step        [value of -window_len]
-kmer_len           [3]
-out_dir            [out]
-canonical          [Boolean] Combine complementary sequence to one.
-mer_sotre          [''] Overlap -kmer_len. Should not be used with -canonical!
-verbose            [Boolean] Print log.
################################################################################
HH
  exit(1);
}
-t and !@ARGV and &usage();

$opts{'window_len'} //= 10e3;
$opts{'kmer_len'}   //= 3;
$opts{'out_dir'}    //= 'out';

-e $opts{'out_dir'} or -l $opts{'out_dir'} or mkdir($opts{'out_dir'}, 0755);

$opts{'kmer_len'} > 0 or &stopErr("[Err] -kmer_len should be larger than 0\n");
$opts{'window_len'} >= $opts{'kmer_len'} or &stopErr("[Err] -window_len should not be smaller than -kmer_len\n");
$opts{'window_step'} //= $opts{'window_len'};

my $iFh = \*STDIN;
scalar(@ARGV) > 0 and $iFh = &openFH(shift(@ARGV));
$opts{'verbose'} and &tsmsg("[Msg] Loading input sequences\n");
my %seq = %{$fa_obj->save_seq_to_hash('faFh' => $iFh)};
for (keys %seq) {
  $seq{$_}{'seq'} =~ s!\s!!g;
  $seq{$_}{'seq'} = uc($seq{$_}{'seq'});
  $seq{$_}{'len'} = length($seq{$_}{'seq'});
}

# Set up possible k-mers
$opts{'verbose'} and &tsmsg("[Msg] Setting up kmers\n");
my @kmer_store;
if (defined $opts{'mer_store'} and $opts{'mer_store'} !~ m!^\s*$!) {
  my %has;
  my $klen;
  for (split(/[\s,.]+/, $opts{'mer_store'})) {
    $_ = uc($_);
    defined $has{$_} and next;
    push(@kmer_store, $_);
    $klen //= length($_);
    $klen == length($_) or &stopErr("[Err] -mer_store must provide kmers with equal length!\n");
  }
  $opts{'kmer_len'} = $klen;
} else {
  @kmer_store = qw/A T G C/;
  for (my $i=1; $i<$opts{'kmer_len'}; $i++) {
    my @new_store;
    for (@kmer_store) {
       for my $b (qw/A T G C/) {
         push(@new_store, "${_}$b");
       }
    }
    @kmer_store = @new_store;
  }
}
### Shrink kmers to canonical ones.
if ($opts{'canonical'}) {
  $opts{'verbose'} and &tsmsg("[Msg]   Getting canonical kmers\n");
  my %k2c_mer; # original kmer => canonical kmer.
  my %has;
  my @new_store;
  for (@kmer_store) {
    my $a = &get_canonical($_);
    $k2c_mer{$_} = $a;
    defined $has{$a} and next;
    push(@new_store, $a);
    $has{$a} = 1;
  }
  @kmer_store = @new_store;
}
my $kstore_cnt = scalar(@kmer_store);

# Index kmers.
$opts{'verbose'} and &tsmsg("[Msg] Indexing kmers\n");
my %kmer_str2idx;
for (my $i=0; $i<@kmer_store; $i++) {
  $kmer_str2idx{$kmer_store[$i]} = $i;
  my $rc = &get_rc($kmer_store[$i]);
  $kmer_str2idx{$rc} = $i;
}

# Count k-mers
# my %kcnt;
for my $k (sort keys %seq) {
  $opts{'verbose'} and &tsmsg("[Msg] Precessing [$k]\n");
  # my %kcntS;
  my $oFh1 = &openFH("$opts{'out_dir'}/$k.tbl", '>');
  print {$oFh1} join("\t", qw/seqID  windowS  windowE ttlKcnt/, @kmer_store)."\n";
  for (my $i=0; $i*$opts{'window_step'} < $seq{$k}{'len'}; $i++) {
    my $startP = $i*$opts{'window_step'}+1;
    my $endP   = $i*$opts{'window_step'}+$opts{'window_len'};
    $endP > $seq{$k}{'len'} and $endP = $seq{$k}{'len'};
    my $t_len = $endP-$startP+1;
    $t_len < 0.9 * $opts{'window_len'} and last;
    my @kcntSW = (0) x $kstore_cnt;
    my $t_seq = substr($seq{$k}{'seq'}, $startP-1, $t_len);
    my $ttlKcnt = 0;
    for (my $j=0; $j<$t_len-$opts{'kmer_len'};$j++) {
      my $t_ele = substr($t_seq, $j, $opts{'kmer_len'});
      $t_ele =~ m!N|X!i and next;
      $ttlKcnt ++;
      defined $kmer_str2idx{$t_ele} or next;
      $kcntSW[ $kmer_str2idx{$t_ele} ] ++;
    }
    print {$oFh1} join("\t", $k, $startP, $endP, $ttlKcnt, @kcntSW)."\n";
  }
  close($oFh1);
}
$opts{'verbose'} and &tsmsg("[Msg] $0 done\n");

sub get_canonical {
  my ($s) = @_; # Must be upper case!
  my $s1 = &get_rc($s);
  my $a = ($s cmp $s1);
  $a <= 0 and return($s);
  return($s1);
}# get_canonical()

sub get_rc {
  my $s = shift;
  $s = reverse($s);
  $s =~ tr/acgturykmbvdhACGTURYKMBVDHwWsSnN/tgcaayrmkvbhdTGCAAYRMKVBHDwWsSnN/;
  return($s);
}# get_rc()

