#!/usr/bin/perl
use strict; 
use fileSunhh; 
use LogInforSunhh; 
use mathSunhh; 
use wm97Sunhh; 
use SNP_tbl; 

!@ARGV and die "perl $0   out_pref   apple.snp_addGmP   idv_list_1_objPop   idv_list_2_refPop [fn_chrID2num 30]\n"; 

my $opref  = shift; 
my $fn_snp = shift; 
my $fn_lis1 = shift; 
my $fn_lis2 = shift;
my $fn_chrID2num = shift; 
$fn_chrID2num //= ''; 
my $cpuN    = shift; 
$cpuN      //= 30; 


my %glob; 
my $gmCn = 2; 


if ( defined $fn_chrID2num and $fn_chrID2num ne '' ) {
  my $fh = &openFH( $fn_chrID2num , '<' ); 
  while (&wantLineC($fh)) {
    my @ta = &splitL( "\t", $_ ); 
    ( defined $ta[0] and defined $ta[1] ) or next; 
    $glob{'chr_id2num'}{$ta[0]} = $ta[1]; 
  }
  close($fh); 
}
my %indv_1 = &indv_list($fn_lis1); 
my %indv_2 = &indv_list($fn_lis2); 

my $fh_snp = &openFH( $fn_snp, '<' ); 

if ($cpuN <= 1) {
  &proc_file($fh_snp, "$opref", \%indv_1, \%indv_2, 1); 
} else {
  my $wrk_d0 = &fileSunhh::new_tmp_dir( 'create' => 1 );
  # Use parallele threading. 
  my $pm = &LogInforSunhh::get_pm( $cpuN ); 
  # Divide input. 
  my @sub_f0 = &fileSunhh::dvd_file( $fh_snp, $cpuN, 'keep_order' => 1, 'with_header' => 1, 'sub_pref' => "$wrk_d0/sub_", 'tmpFile' => "$wrk_d0/base_0" );
  close($fh_snp); 

  # Make output file. 
  for my $f0 (@sub_f0) {
    my $pid = $pm->start and next; 
    my $fh0 = &openFH($f0, '<');
    my $o_pref = "$f0.o"; 
    # Process $f0 to generate files:
    #   $f0.o_g1.geno
    #   $f0.o_g2.geno
    #   $f0.o.snp
    &proc_file($fh0, $o_pref, \%indv_1, \%indv_2, 0); 
    $pm->finish; 
  }
  $pm->wait_all_children; 

  # Combine the outputs
  my %fho; 
  for my $f0 (@sub_f0) {
    open FG1,'<',"$f0.o_g1.geno" or die; 
    open FG2,'<',"$f0.o_g2.geno" or die; 
    open FGG,'<',"$f0.o.snp" or die; 
    while (my $lg1 = <FG1>) {
      my $lg2 = <FG2>; 
      my $lgg = <FGG>; 
      $lgg =~ m!^(\S+)_\d+\t! or die "bad line: $lgg\n"; 
      my $id = $1; 
      unless (defined $fho{$id}) {
        $fho{$id}{'geno1'} = &openFH("${opref}.${id}_g1.geno", '>'); 
        $fho{$id}{'geno2'} = &openFH("${opref}.${id}_g2.geno", '>'); 
        $fho{$id}{'SNP'}   = &openFH("${opref}.${id}.snp", '>'); 
      }
      print {$fho{$id}{'geno1'}} $lg1; 
      print {$fho{$id}{'geno2'}} $lg2; 
      print {$fho{$id}{'SNP'}}   $lgg; 
    }
    close(FG1); 
    close(FG2); 
    close(FGG); 
  }# for my $f0 (@sub_f0) 
  for my $k1 (keys %fho) {
    for my $k2 (qw/geno1 geno2 SNP/) {
      close($fho{$k1}{$k2}); 
    }
  }# for my $k1 (keys %fho)
  &fileSunhh::_rmtree($wrk_d0); 
}

&tsmsg("[Rec] Done. $0\n"); 

#my $fh_oSite = &openFH("${opref}.sites", '>'); 
#my $fh_oLoci = &openFH("${opref}.locs",  '>'); 
sub indv_list {
  my $fn = shift; 
  my %back; 
  my $fh = &openFH($fn, '<'); 
  while (&wantLineC($fh)) {
    my @ta=&splitL("\t", $_); 
    defined $back{'has'}{$ta[0]} and next; 
    push(@{$back{'arr'}}, $ta[0]); 
    $back{'has'}{$ta[0]} = $#{$back{'arr'}}; 
  }
  close ($fh); 
  return(%back); 
}


sub load_gmP {
  my $fn = shift; 
  my %back; 
  my $fh = &openFH($fn, '<'); 
  while (&wantLineC($fh)) {
    # chr \\t pos \\t cM 
    my @ta = &splitL("\t", $_); 
    $ta[0] eq 'chr' and next; 
    $back{$ta[0]}{$ta[1]} = $ta[2]; 
  }
  close($fh); 
  return(%back); 
}

sub proc_file {
  my ($fh_snp, $o_pref, $indv_1h, $indv_2h, $useID) = @_; 
  $useID //= 0; 
  # If useID == 1, the output files will be "${o_pref}.chrID_g1.geno" and "${o_pref}.chrID.snp"; 
  # If useID == 0, the output files will be "${o_pref}_g1.geno" and "${o_pref}.snp"; 
  my %indv_1 = %$indv_1h; 
  my %indv_2 = %$indv_2h; 
  my %fho; 
  my @h; 
  { my $a=&wantLineC($fh_snp); @h=split(/\t/, $a); } 
  
  {
    my %nn; 
    my @new_h; 
    for (my $i=0; $i<@h; $i++) {
  
      defined $indv_1{'has'}{$h[$i]} and push( @{$indv_1{'goodIdx'}}, $i ); 
      defined $indv_2{'has'}{$h[$i]} and push( @{$indv_2{'goodIdx'}}, $i ); 
      defined $indv_1{'has'}{$h[$i]} and defined $indv_2{'has'}{$h[$i]} and &stopErr("[Err] Indv [$h[$i]] exists in both groups.\n"); 
    }
    @{$indv_1{'goodSample'}} = @h[ @{$indv_1{'goodIdx'}} ]; 
    @{$indv_2{'goodSample'}} = @h[ @{$indv_2{'goodIdx'}} ]; 
  }
  
  my %cc = ( "cntN_base"=>0, "cntN_step"=>1e4 ); 
  
  my %locs; 
  my %sites; 
  my %prev; 
  while (&wantLineC($fh_snp)) { 
    &fileSunhh::log_section($. , \%cc) and &tsmsg("[Msg] $. line.\n"); 
    my @ta=split(/\t/, $_); 
  #  defined $phy2gm_P{$ta[0]} or next; 
  #  defined $phy2gm_P{$ta[0]}{$ta[1]} or next; 
  
    my @ta_1 = @ta[ @{$indv_1{'goodIdx'}} ]; 
    my @ta_2 = @ta[ @{$indv_2{'goodIdx'}} ]; 
  
    my ( %al, %al_1, %al_2 ); 
    my %geno; 
    my $is_bad = 0; 
    for (my $i=0; $i<@ta_1; $i++ ) { 
      if ($ta_1[$i] =~ m/^([ATGCN])$/) { 
        $al_1{$1} += 2; 
        $al{$1}   += 2; 
      } elsif ($ta_1[$i] =~ m/^([ATGC])([ATGC])$/) { 
        $al_1{$1} ++; $al_1{$2} ++; 
        $al{$1}   ++; $al{$2}   ++; 
      } elsif ( $ta_1[$i] =~ m/^[ATGC]{3,}$/ ) { 
        $is_bad = 1; 
        last; 
      } else { 
        my @bb = &SNP_tbl::dna_d2b( $ta_1[$i] ); 
        if ( @bb == 2 ) {
          $al_1{$bb[0]} ++; $al_1{$bb[1]} ++; 
          $al{$bb[0]}   ++; $al{$bb[1]}   ++; 
          $ta_1[$i] = "$bb[0]$bb[1]"; 
        } else {
          if ( !(defined $glob{'bad_geno'}{$ta_1[$i]}) ) {
            $glob{'bad_geno'}{$ta_1[$i]} = 1; 
            &tsmsg("[Wrn] Skip site with bad genotype ta_1 [$ta_1[$i]][@bb]\n"); 
          }
          $is_bad = 1; 
          last; 
        }
      }
    }
    delete $al_1{'N'}; 
    scalar(keys %al_1) > 0 or $is_bad = 1; 
    $is_bad == 1 and next; 
    for (my $i=0; $i<@ta_2; $i++ ) { 
      if ($ta_2[$i] =~ m/^([ATGCN])$/) { 
        $al_2{$1} += 2; 
        $al{$1}   += 2; 
      } elsif ($ta_2[$i] =~ m/^([ATGC])([ATGC])$/) { 
        $al_2{$2} ++; $al_2{$2} ++; 
        $al{$2}   ++; $al{$2}   ++; 
      } elsif ( $ta_2[$i] =~ m/^[ATGC]{3,}$/ ) { 
        $is_bad = 1; 
        last; 
      } else {
        my @bb = &SNP_tbl::dna_d2b( $ta_2[$i] ); 
        if ( @bb == 2 ) {
          $al_2{$bb[0]} ++; $al_2{$bb[1]} ++; 
          $al{$bb[0]}   ++; $al{$bb[1]}   ++; 
          $ta_2[$i] = "$bb[0]$bb[1]"; 
        } else {
          if ( !(defined $glob{'bad_geno'}{$ta_2[$i]}) ) {
            $glob{'bad_geno'}{$ta_2[$i]} = 1; 
            &tsmsg("[Wrn] Skip site with bad genotype ta_2 [$ta_2[$i]][@bb]\n"); 
          }
          $is_bad = 1; 
          last; 
        }
      }
    }
    delete $al_2{'N'}; 
    scalar(keys %al_2) > 0 or $is_bad = 1; 
    $is_bad == 1 and next; 
    delete $al{'N'}; 
    scalar( keys %al ) == 2 or next; 
    my @aa = sort { $al{$b} <=> $al{$a} || $a cmp $b } keys %al; 
    $geno{$aa[0]}         = '1 1'; 
    $geno{"$aa[0]$aa[0]"} = '1 1'; 
    $geno{$aa[1]}         = '0 0'; 
    $geno{"$aa[1]$aa[1]"} = '0 0'; 
    $geno{"$aa[0]$aa[1]"} = '1 0'; 
    $geno{"$aa[1]$aa[0]"} = '1 0'; 
    $geno{"N"}            = '9 9'; 
  
    
    my $chrID = $ta[0]; 
    my $chrN ; 
    if ( defined $glob{'chr_id2num'}{$chrID} ) {
      $chrN = $glob{'chr_id2num'}{$chrID}; 
    } else {
      $chrN = &wm97Sunhh::chrID_to_number( $chrID , 'WM97_Chr'); 
      $chrN =~ m!^\d+$! or &stopErr("[Err] Failed to convert chrID [$chrID] to number [$chrN]\n"); 
      defined $glob{'chr_num2id'}{$chrN} and &stopErr("[Err] Repeat chrN [$chrN] for differnt chrID [$glob{'chr_num2id'}{$chrN} $chrID]\n"); 
      $glob{'chr_id2num'}{$chrID} = $chrN; 
      $glob{'chr_num2id'}{$chrN}  = $chrID; 
    }
  
    $chrN =~ m!^\d+$! or &stopErr("[Err] Bad chrID [$chrN] from [$ta[0]]\n"); 
  
    my $gmP = $ta[$gmCn]; 
    if (defined $prev{'gmID'} and $prev{'gmID'} eq $ta[0]) {
      $prev{'gmP'} < $gmP or next; 
    }
    $prev{'gmID'} = $ta[0]; 
    $prev{'gmP'}  = $gmP; 
  
    if ($useID == 0) {
      unless ( defined $fho{'geno1'} ) {
        $fho{'geno1'} = &openFH("${o_pref}_g1.geno", '>'); 
        $fho{'geno2'} = &openFH("${o_pref}_g2.geno", '>'); 
        $fho{'SNP'}   = &openFH("${o_pref}.snp", '>'); 
      }
      print {$fho{'geno1'}}  join(' ', map { $geno{$_} } @ta_1)."\n"; 
      print {$fho{'geno2'}}  join(' ', map { $geno{$_} } @ta_2)."\n"; 
      my $mrkID = "$ta[0]_$ta[1]";
      print {$fho{'SNP'}}    join( "\t", $mrkID, $chrN, $gmP, $ta[1], $aa[0], $aa[1] )."\n"; 
    } else {
      unless ( defined $fho{$chrID}{'geno1'} ) {
        $fho{$chrID}{'geno1'} = &openFH("${o_pref}.${chrID}_g1.geno", '>'); 
        $fho{$chrID}{'geno2'} = &openFH("${o_pref}.${chrID}_g2.geno", '>'); 
        $fho{$chrID}{'SNP'}   = &openFH("${o_pref}.${chrID}.snp", '>'); 
      }
      print {$fho{$chrID}{'geno1'}}  join(' ', map { $geno{$_} } @ta_1)."\n"; 
      print {$fho{$chrID}{'geno2'}}  join(' ', map { $geno{$_} } @ta_2)."\n"; 
      my $mrkID = "$ta[0]_$ta[1]";
      print {$fho{$chrID}{'SNP'}}    join( "\t", $mrkID, $chrN, $gmP, $ta[1], $aa[0], $aa[1] )."\n"; 
    }

  }
  if ($useID == 0) {
    for my $k1 (qw/geno1 geno2 SNP/) {
      close($fho{$k1}); 
    }
  } else {
    for my $k1 (keys %fho) {
      for my $k2 (qw/geno1 geno2 SNP/) {
        close($fho{$k1}{$k2}); 
      }
    }
  }
  close($fh_snp); 
}# proc_file() 
