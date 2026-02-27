#!/usr/bin/perl
# 260227: Add '--double-id' to PLINK command.
use strict;
use warnings;
use fileSunhh;
use LogInforSunhh;

!@ARGV and die "perl $0 outPrefix in.peak in-ld.vcf.gz minR2 flankLength\n";

my $opref = shift;
my $peakFn = shift;
my $vcfFn  = shift;
my $minR2 = shift; # 0.8
my $flankLen = shift; # 500e3
my $flankLK = $flankLen/1e3;

my $ofh1 = &openFH("${opref}_LD.mrk.gz", '>');
print {$ofh1} join("\t", qw/peakID chrV proxPosi r2 mrkID dist/)."\n";
my $ofh2 = &openFH("${opref}_LD.blk", '>');
print {$ofh2} join("\t", qw/peakID chrV start end span/)."\n";
my $ofh3 = &openFH("${opref}_LD.tagVAR", '>');
print {$ofh3} join("\t", qw/peakID chrV proxPosi r2 mrkID dist/)."\n";
my %marker_map;
for my $l1 (&fileSunhh::load_tabFile($peakFn)) {
  my @ld_markers = ($l1->[3]); # 3 for fastlmm.res
  my $s=$l1->[3]-$flankLen; # fastlmm
  my $e=$l1->[3]+$flankLen; # fastlmm
  # my $varID = (defined $l1->[8] and $l1->[8] ne "") ? $l1->[8] : $l1->[2] ; # emmax
  my $varID = $l1->[0]; # fastlmm

  # Change SV ID from '10_8886884_v41319259_22' to 'v41319259_22' when GWAS and VCF have different format of IDs.
  # if ($varID =~ m!^\d+_\d+_(v\d+_\d+)$!) {
  #   $marker_map{$1} = $varID;
  #   $varID = $1;
  # }
  &runCmd("bcftools view $vcfFn -r $l1->[1]:$s-$e -Ov -o $opref.t.vcf --threads 20");
  &runCmd("plink --vcf $opref.t.vcf --output-missing-genotype 0 --make-bed --out $opref.t --geno 0.8 --double-id");
  # Test if associate SNP is in the data.
  &fileSunhh::write2file("$opref.t.chk.txt", "$varID\n", '>');
  system("plink --bfile $opref.t --extract $opref.t.chk.txt --write-snplist --out $opref.t.out --double-id\n");
  if (-f "$opref.t.out.snplist") {
    &runCmd("rm -f $opref.t.out* $opref.t.chk.txt");
    &runCmd("plink --bfile $opref.t --r2 --ld-window $flankLen --ld-window-r2 0 --ld-window-kb $flankLK --out ${opref}_LD --ld-snp $varID --double-id");
  } else {
    next;
  }
  my $ifh1 = &openFH("${opref}_LD.ld", '<');
  # CHR_A         BP_A          SNP_A  CHR_B         BP_B          SNP_B           R2 
  #     2     31961541   v9279436_829      2     31461557   2_31461557_1      0.40214 
  #     2     31961541   v9279436_829      2     31461558   2_31461558_1     0.402228 
  while (my $l2 = <$ifh1>) {
    $l2 =~ s!^\s+!!g;
    chomp($l2);
    my @ta=split(/\s+/, $l2);
    $ta[0] eq 'CHR_A' and next;
    # print {$ofh1} join("\t", $l1->[2], @ta[3,4,6,5], $ta[4]-$l1->[1])."\n"; # emmax
    print {$ofh1} join("\t", $l1->[0], @ta[3,4,6,5], $ta[4]-$l1->[3])."\n"; # fastlmm
    $ta[6] >= $minR2 or next;
    push(@ld_markers, $ta[4]);
    print {$ofh3} join("\t", $l1->[0], @ta[3,4,6,5], $ta[4]-$l1->[3])."\n"; # fastlmm
    # if ($ta[5] =~ m!^(\d+_\d+_)?v\d+!) {
    #   # print {$ofh3} join("\t", $l1->[2], @ta[3,4,6,5], $ta[4]-$l1->[1])."\n"; # emmax
    #   print {$ofh3} join("\t", $l1->[0], @ta[3,4,6,5], $ta[4]-$l1->[3])."\n"; # fastlmm
    # }
  }
  close($ifh1);
  @ld_markers = sort {$a<=>$b} @ld_markers;
  # print {$ofh2} join("\t", $l1->[2], $l1->[0], $ld_markers[0], $ld_markers[-1], $ld_markers[-1]-$ld_markers[0]+1)."\n"; # emmax
  print {$ofh2} join("\t", $l1->[0], $l1->[1], $ld_markers[0], $ld_markers[-1], $ld_markers[-1]-$ld_markers[0]+1)."\n"; # fastlmm
}
close($ofh1);
close($ofh2);
close($ofh3);

