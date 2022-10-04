#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
  "help!", 
  "shrt_len:i", # [9]
  "out_id_list:s", # [out_id_list]
  "fmt:s", # 0
); 
use fastaSunhh; 
my $fs_obj = fastaSunhh->new(); 
$opts{'shrt_len'} //= 9; 
$opts{'out_id_list'} //= 'out_id_list'; 

my $help_txt = <<HH; 

perl $0 in_aligned.fasta > out_phyml_interleaved.phyml

-help

-shrt_len       [$opts{'shrt_len'}]
-out_id_list    ['out_id_list']
-fmt            [0] '0' for normal format;
                    '1' for mcmctree input;
                    '2' for mcmctree input with three frames.

HH

!@ARGV and &LogInforSunhh::usage($help_txt); 
$opts{'help'} and &LogInforSunhh::usage($help_txt); 
$opts{'fmt'} //= 0;

&tsmsg("[Msg] Reading fasta.\n\n"); 
my $charN_per_line = 60; 

my %fas_hash = %{$fs_obj->save_seq_to_hash('faFile'=>$ARGV[0], 'has_head'=>1)}; 
my @seq_IDs = sort { $fas_hash{$a}{'Order'} <=> $fas_hash{$b}{'Order'} } keys %fas_hash; 
for (@seq_IDs) {
  $fas_hash{$_}{'seq'} =~ s/\s//g; 
  $fas_hash{$_}{'len'} = length($fas_hash{$_}{'seq'}); 
}
my $seq_num = scalar(@seq_IDs); 
my $seq_len = $fas_hash{ $seq_IDs[0] }{'len'}; 

&tsmsg("[Msg] Printint output.\n");
if ($opts{'fmt'} == 0 or $opts{'fmt'} == 1) {
  print STDOUT "$seq_num $seq_len\n";
} elsif ($opts{'fmt'} == 2) {
  ;
} else {
  &stopErr("[Err] Reach here!\n");
}
my %used_ID; 

open O,'>',"$opts{'out_id_list'}" or &stopErr("$!\n");
if ($opts{'fmt'} == 0) {
  for (my $i=0; $i<$seq_len; $i+=$charN_per_line) {
    if ($i==0) {
      for (@seq_IDs) {
        my $new_ID = substr($_, 0, $opts{'shrt_len'}); $new_ID =~ s/[\s();:]/_/g;
        my $tk = $new_ID; 
        my $suff = 'a'; 
        while ( defined $used_ID{$tk} ) {
          $suff eq 'z' and &stopErr("[Err] Too many repeat names [$new_ID]\n"); 
          $suff ++; 
          $tk = "$new_ID$suff"; 
        }
        $used_ID{$tk} = 1;
        print O "$tk\t$_\n"; 
        my $l2 = $opts{'shrt_len'}+1; 
        $new_ID = sprintf("%-${l2}s", $tk); 
        print STDOUT join(' ', $new_ID, @{ &arr_by_seg( substr($fas_hash{$_}{'seq'}, $i, $charN_per_line), 10 ) })."\n"; 
      }
      print STDOUT "\n"; 
      next; 
    }
    for (@seq_IDs) {
      print STDOUT join(' ', ' ' x 10, @{ &arr_by_seg( substr($fas_hash{$_}{'seq'}, $i, $charN_per_line), 10 ) })."\n"; 
    }
    print STDOUT "\n"; 
  }
} elsif ($opts{'fmt'} == 1) {
  for my $k1 (@seq_IDs) {
    my $new_ID = substr($k1, 0, $opts{'shrt_len'}); $new_ID =~ s/[\s();:]/_/g;
    my $tk = $new_ID;
    my $suff = 'a';
    while ( defined $used_ID{$tk} ) {
      $suff eq 'z' and &stopErr("[Err] Too many repeat names [$new_ID]\n"); 
      $suff ++; 
      $tk = "$new_ID$suff"; 
    }
    $used_ID{$tk} = 1;
    print O "$tk\t$k1\n";
    my $l2 = $opts{'shrt_len'} + 1;
    $new_ID = sprintf("%-${l2}s", $tk);
    print STDOUT join(" ", $new_ID, $fas_hash{$k1}{'seq'})."\n";
  }
  print STDOUT "\n";
} elsif ($opts{'fmt'} == 2) {
  my @s2;
  my @id2;
  for (my $i=0; $i<@seq_IDs; $i++) {
    my $k1 = $seq_IDs[$i];
    my $new_ID = substr($k1, 0, $opts{'shrt_len'}); $new_ID =~ s/[\s();:]/_/g;
    my $tk = $new_ID;
    my $suff = 'a';
    while ( defined $used_ID{$tk} ) {
      $suff eq 'z' and &stopErr("[Err] Too many repeat names [$new_ID]\n"); 
      $suff ++; 
      $tk = "$new_ID$suff"; 
    }
    $used_ID{$tk} = 1;
    print O "$tk\t$k1\n";
    my $l2 = $opts{'shrt_len'} + 1;
    $new_ID = sprintf("%-${l2}s", $tk);
    push(@id2, $new_ID);
    for (my $j=0; $j<$fas_hash{$k1}{'len'}; $j+=3) {
      for (my $m=0; $m<3; $m++) {
        my $jm= $j+$m;
        $s2[$m][$i] //= "";
        $s2[$m][$i] .= substr($fas_hash{$k1}{'seq'}, $jm, 1);
      }
    }
  }
  for (my $m=0; $m<3; $m++) {
    my $v1 = scalar(@{$s2[$m]});
    my $v2 = length($s2[$m][0]);
    print STDOUT "$v1 $v2\n";
    for (my $i=0; $i<@id2; $i++) {
      print STDOUT join(" ", $id2[$i], $s2[$m][$i])."\n";
    }
    print STDOUT "\n";
  }
} else {
  &stopErr("[Err] Should not arrive here!\n");
}
close (O); 

sub arr_by_seg {
  # $_[0] whole sequence 
  # $_[1] length per segment. 
  my @back; 
  for (my $i=0; $i<length($_[0]); $i+=$_[1]) {
    my $ss = substr($_[0], $i, $_[1]); 
    push(@back, $ss); 
  }
  return \@back; 
}




