#!/usr/bin/perl
# 20211130: Apply QUAST to remove duplicated fragments within and among input fasta files. 
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use fastaSunhh; 
use Parallel::ForkManager; 
use Getopt::Long; 
my $fas_obj = fastaSunhh->new(); 

my %opts; 
GetOptions(\%opts, 
  "in_fa_list:s", # in_fa_list # 20211124_asm/Beauregard_hiCanuUtg.fa.gz 
  "out_prefix:s", 
  "exe_quast:s", # quast.py
  "para_quast:s", # -t 50 --fragmented --min-identity 90.0 --min-alignment 300 --fast 
  "min_kept_len:i", # 300; 
  "joinFinal!", 

  "rmSelfDup!", # With this parameter, I'll try to remove duplications within each file after normal running. 
  "selfSplitN:i",  # 10; 

  "cpuN:i", "nprocF:s", "wait_sec:i",   # 
  "help!", 
); 

my $help_txt = <<HH; 
##########################################################################################
## Apply QUAST to remove duplicated fragments within and among input fasta files. 
perl $0 -in_fa_list in_fa_list -out_prefix out_prefix

  -out_prefix    [opref] The final output will be 'opref.pan/final_xx.fa'. Additional file 'opref.consensus.fa' will be available if '-joinFinal' exists. 
  -in_fa_list    [fasta_file_list] The first column is used to find input fasta files. Could be one line only. 

  -exe_quast     [quast.py] Could be '/data/Sunhh/src/General/anaconda/install/envs/quast/bin/quast.py'
  -para_quast    [ -t 50 --fragmented --min-identity 90.0 --min-alignment 300 --min-contig 300 --fast --eukaryote --no-html --no-icarus ]
  -min_kept_len  [300] The minimal length accepted to keep in unaligned sequences. This should be no less than "--min-alignment" value used in QUAST. 

  -joinFinal     [Boolean] Generate a file having all the sequences in opref.pan/

  -rmSelfDup     [Boolean] Given this parameter, deduplication within each file will be applied after performed between files. Please make sure "-min_kept_len" is no less than QUAST parameter "--min-contig". 
  -selfSplitN    [20] 

  -cpuN          [1] This is only useful when "-rmSelfDup" exists. 
  -nprocF        [Nproc]
  -wait_sec      [0]

  -help          [Boolean]
##########################################################################################
HH

defined $opts{'in_fa_list'} or &LogInforSunhh::usage($help_txt); 

$opts{'out_prefix'} //= 'opref'; 
$opts{'exe_quast'}  //= 'quast.py'; 
$opts{'para_quast'} //= ' -t 50 --fragmented --min-identity 90.0 --min-alignment 300 --min-contig 300 --fast --eukaryote '; 
$opts{'selfSplitN'} //= 20; 
$opts{'min_kept_len'} //= 300; 
$opts{'cpuN'}     //= 1; 
$opts{'nprocF'}   //= 'Nproc'; 
$opts{'wait_sec'} //= 0; 

my $MAX_PROCESSES = $opts{'cpuN'} ; # Sometimes $parm{'cpuN'} - 1 may be better.
my $pm = new Parallel::ForkManager($MAX_PROCESSES); 


my $in_param = ''; 
{
  my %txtParam = map { $_ => 1 } qw(exe_quast para_quast); 
  my %valParam = map { $_ => 1 } qw(selfSplitN min_kept_len); 
  my %booParam = map { $_ => 1 } qw(rmSelfDup); 
  PARAM: 
  for my $k1 (sort keys %opts) {
    defined $txtParam{$k1} and $in_param .= " -$k1 '$opts{$k1}' "; 
    defined $valParam{$k1} and $in_param .= " -$k1 $opts{$k1} "; 
    defined $booParam{$k1} and $in_param .= " -$k1 "; 
  }
}

# Load sequences in each fa file. 
my (@faList, @faSeqs); 
{
  &tsmsg("[Msg] Loading fasta files.\n"); 
  my $fh_fa_list = &openFH($opts{'in_fa_list'}); 
  while (<$fh_fa_list>) {
    chomp;
    my @ta=split(/\t/, $_); 
    push(@faList, $ta[0]); 
    push(@faSeqs, $fas_obj->save_seq_to_hash('faFile' => $ta[0]) ); 
  }
  close($fh_fa_list); 
  for (@faSeqs) {
    for my $k (keys %$_) {
      $_->{$k}{'seq'} =~ s!\s!!g; 
      $_->{$k}{'len'} = length( $_->{$k}{'seq'} ); 
    }
  }
}
my @currList = @faList; 


# Create working directory. 
my $wdir = &fileSunhh::new_tmp_dir( 'create' => 1 ); 
mkdir("$wdir/oriseq/"); 
mkdir("$wdir/pan/"); 

# Write original sequence files with modified file names for QUAST runs. 
{
  &tsmsg("[Msg] Writing original fasta files.\n"); 
  for (my $i=0; $i<@faList; $i++) {
    my $fh_otmp = &openFH("$wdir/oriseq/$i.fa", '>'); 
    for my $k1 (sort { $faSeqs[$i]{$b}{'len'} <=> $faSeqs[$i]{$a}{'len'} } keys %{$faSeqs[$i]}) {
      print {$fh_otmp} ">$k1\n$faSeqs[$i]{$k1}{'seq'}\n"; 
    }
    close($fh_otmp); 
    $currList[$i] = "$wdir/oriseq/$i.fa"; 
  }
}

# Determine the cycling and run QUAST. 
&tsmsg("[Msg] Begin to run QUAST.\n"); 
for (my $i=0; $i<$#faList; $i++) {
  &tsmsg("[Msg] Using [$i] as reference.\n"); 
  # Reference fa: $currList[$i]
  # Query fa: @currList[$i+1 .. $#currList]
  # QUAST command: quast.py asm1.fa asm2.fa asm3.fa -R ref.fa --fragmented --min-identity 90.0 --min-alignment 300 --fast

  &fileSunhh::_copy($currList[$i], "$wdir/pan/final_${i}.fa"); 

  # Run QUAST 
  if (-z $currList[$i]) {
    next; 
  }
  -e "$wdir/quastrun/$i/" and &fileSunhh::_rmtree("$wdir/quastrun/$i/"); 
  my $cmd = "$opts{'exe_quast'} $opts{'para_quast'} -o $wdir/quastrun/$i/ -R $currList[$i] "; 
  my @good_fc; 
  for my $fc1 (@currList[$i+1 .. $#currList]) {
    -z $fc1 and next; 
    push(@good_fc, $fc1); 
    # $cmd .= " $fc1 "; 
  }
  if (scalar(@good_fc) >= 1) {
    $cmd .= " "; 
    $cmd .= join(" ", @good_fc); 
    $cmd .= " "; 
  } else {
    # There is no good sample file, so we should skip this QUAST. 
    # We should exit this for loop, but we don't do this because we want to copy the files to pan/final_; 
    next; 
  }
  # $cmd .= join(' ', '', @currList[$i+1 .. $#currList]); 
  $cmd .= " 1> $wdir/quastrun_std.$i 2> $wdir/quastrun_err.$i"; 
  &runCmd($cmd); 
  # Extract unaligned sequences one by one (modified ID)
  #   file path: quast_results/latest/contigs_reports/contigs_report_y.unaligned.info
  #   It seems that: 
  #     (1) When there is no file 'contigs_reports/contigs_report_y.unaligned.info' created, none of the sequences in the sample file can be aligned to the reference. 
  #     (2) When there is a file 'contigs_repports/contigs_report_y.unaligned.info' created, those sequences that are not aligned at all will be assigned as 'full' in 'Unaligned_type'.
  #         I'd like not to add start-to-end information for these fully unaligned sequences. 
  for (my $j=$i+1; $j<=$#currList; $j++) {
    my $fn_info = "$wdir/quastrun/$i/contigs_reports/contigs_report_${j}.unaligned.info"; 
    if (-f $fn_info) {
      # Extract unaligned sequence to fasta file in currList, and then update @faSeqs; 
      open F1,'<',"$fn_info" or die; 
      open F2,'>',"$currList[$j]" or die; 
      <F1>; 
      my %tmpSeqs; 
      while (<F1>) {
        chomp; 
        my @ta = split(/\t/, $_); 
        if ($ta[1] == $ta[2]) {
          # There are many ways to judge this fully unaligned sequences: 
          #   $ta[3] eq 'full'
          #   $ta[4] eq "1-$ta[1]"
          # I just choose one with minimal number comparisons. 
          print F2 ">$ta[0]\n$faSeqs[$j]{$ta[0]}{'seq'}\n"; $tmpSeqs{$ta[0]}{'seq'} = $faSeqs[$j]{$ta[0]}{'seq'}; $tmpSeqs{$ta[0]}{'len'} = $faSeqs[$j]{$ta[0]}{'len'}; 
          next; 
        }
        my @tb = split(/,/, $ta[4]); 
        for my $tc (@tb) {
          $tc =~ m!^(\d+)\-(\d+)$! or do { &tsmsg("[Wrn] Skip unaligned record: $_\n"); next; }; 
          my ($s,$e) = ($1, $2); 
          $e-$s+1 >= $opts{'min_kept_len'} or next; 
          my $newID = &resolve_id($ta[0], $s, $e); 
          my $newSeq = substr($faSeqs[$j]{$ta[0]}{'seq'}, $s-1, $e-$s+1); 
          print F2 ">$newID\n$newSeq\n"; $tmpSeqs{$newID}{'seq'} = $newSeq; $tmpSeqs{$newID}{'len'} = $e-$s+1; 
        }
      }
      close F2; 
      close F1; 
      ## Update @faSeqs; 
      %{$faSeqs[$j]} = (); 
      $faSeqs[$j] = \%tmpSeqs; 
    } else {
      ; # Nothing to do. 
    }
  }
}
&fileSunhh::_copy($currList[-1], "$wdir/pan/final_$#{currList}.fa"); 

# Remove duplication from each file after removing duplicates between files. 
#   I'll replace the files in '$wdir/pan/final_$i.fa' directory with the same file name after cleaning them. 
if ($opts{'rmSelfDup'}) {
  my $selfDir = "$wdir/selfClean/"; 
  mkdir( $selfDir ); 
  # De-duplicate each file one by one.
  #   Get "$selfDir/dedupped/$i.fa"; 
  for (my $i=0; $i<@faSeqs; $i++) {
    $MAX_PROCESSES = &LogInforSunhh::change_procN($pm, $opts{'nprocF'}, $MAX_PROCESSES); 
    if ($opts{'wait_sec'} > 0) { sleep($opts{'wait_sec'}); }
    my $pid = $pm->start and next; 
    mkdir( "$selfDir/$i/" ); # Add this layer for multi-threads running. 
    mkdir( "$selfDir/$i/dedupped/" ); 
    mkdir( "$selfDir/$i/oriseq/" ); 
    # The resulting file "$selfDir/dedupped/$i.fa" should be already self-dedupped. 

    # step 1: Split file "pan/final_$i.fa" into $opts{'selfSplitN'} files for further use. 
    ### If the input file has <=10 sequences or its has <= 50 sequences summing up to <= 1 Mb, I want to separate each sequence to a file. 
    my $seqN_cur = scalar(keys %{$faSeqs[$i]}); # Number of sequences in this file $i. 
    my $seqL_cur = 0; # Total bps in this file. 
    for (keys %{$faSeqs[$i]}) {
      $seqL_cur += $faSeqs[$i]{$_}{'len'}; 
    }
    my $j_time = -1; 

    my @j_rmFn; 
    &fileSunhh::write2file("$selfDir/$i/in_fa_list", "", '>'); push(@j_rmFn, "$selfDir/$i/in_fa_list"); 
    if ( $seqN_cur <= 1 ) {
      &fileSunhh::_move( "$wdir/pan/final_${i}.fa", "$selfDir/$i/dedupped/$i.fa" ); 
      next; 
    }
    if ($seqN_cur <= 10 or ($seqN_cur <= 50 and $seqL_cur <= 1e6)) {
      # Separate each sequence to a single file for redundancy removal. 
      for my $k1 (sort { $faSeqs[$i]{$b}{'len'} <=> $faSeqs[$i]{$a}{'len'} } keys %{$faSeqs[$i]}) {
        $j_time++; 
        &fileSunhh::write2file("$selfDir/$i/oriseq/${j_time}.fa", ">$k1\n$faSeqs[$i]{$k1}{'seq'}\n", '>'); push(@j_rmFn, "$selfDir/$i/oriseq/${j_time}.fa"); 
        &fileSunhh::write2file("$selfDir/$i/in_fa_list", "$selfDir/$i/oriseq/${j_time}.fa\n", '>>'); 
      }
    } else {
      # Separate input file into selfSplitN files with even sequence number. 
      my $j=-1; 
      my $j_step = int($seqN_cur/$opts{'selfSplitN'}+0.5); 
      my $j_ofh; 
      for my $k1 (sort { $faSeqs[$i]{$b}{'len'} <=> $faSeqs[$i]{$a}{'len'} } keys %{$faSeqs[$i]}) {
        $j++; 
        if ($j % $j_step == 0) {
          $j_time = int($j / $j_step); 
          defined $j_ofh and close($j_ofh); 
          open($j_ofh, '>', "$selfDir/$i/oriseq/${j_time}.fa") or &stopErr("[Err] Failed to write to file [$selfDir/$i/oriseq/${j_time}.fa]\n"); push(@j_rmFn, "$selfDir/$i/oriseq/${j_time}.fa"); 
          &fileSunhh::write2file("$selfDir/$i/in_fa_list", "$selfDir/$i/oriseq/${j_time}.fa\n", '>>'); 
        }
        print {$j_ofh} ">$k1\n$faSeqs[$i]{$k1}{'seq'}\n"; 
      }
      defined $j_ofh and close($j_ofh); 
    }# if ($seqN_cur <= 10 or ($seqN_cur <= 50 
    # step 2: Run deduplication with QUAST. 
    &runCmd("perl $0 -joinFinal $in_param -in_fa_list $selfDir/$i/in_fa_list -out_prefix $selfDir/$i/dedupped/$i"); push(@j_rmFn, "$selfDir/$i/dedupped/${i}.pan/"); 
    # step 3: Save deduplicated sequences to one file. 
    &fileSunhh::_move("$selfDir/$i/dedupped/$i.consensus.fa", "$selfDir/$i/dedupped/$i.fa"); 
    # step 4: Clean temporary files. 
    for my $f1 (@j_rmFn) { &fileSunhh::_rmtree($f1); }

    $pm->finish; 
    $MAX_PROCESSES = &LogInforSunhh::change_procN($pm, $opts{'nprocF'}, $MAX_PROCESSES); 
  }# for (my $i=0; $i<@faSeqs; $i++)
  $pm->wait_all_children; 

  # Replace the input fasta file by moving; 
  ### I don't integrate into the former for loop so I am able to check how many files have been deduplicated by checking files in the folder $selfDir/$i/dedupped/ . 
  for (my $i=0; $i<@faSeqs; $i++) { 
    &fileSunhh::_move("$selfDir/$i/dedupped/$i.fa", "$wdir/pan/final_$i.fa"); 
    &tsmsg("[Msg] Self-dedup has been done for in file [$wdir/pan/final_$i.fa] for [$faList[$i]]\n"); 
  } 
}# if rmSelfDup


# Move final sequences files out. I want to keep file information for further use. 
&tsmsg("[Msg] Moving pan directory out.\n"); 
&fileSunhh::_move("$wdir/pan/", "$opts{'out_prefix'}.pan/"); 

# Create final joined fasta file if required. 
if ($opts{'joinFinal'}) {
  &tsmsg("[Msg] Join sequences.\n"); 
  open(OJ1,'>',"$opts{'out_prefix'}.consensus.fa") or &stopErr("[Err] Failed to write file [$opts{'out_prefix'}.consensus.fa]\n"); 
  for (my $i=0; $i<@faList; $i++) {
    open(F1,'<',"$opts{'out_prefix'}.pan/final_$i.fa") or &stopErr("[Err] Failed to read file [$opts{'out_prefix'}.pan/final_$i.fa]\n"); 
    while (<F1>) { print OJ1 $_; }
    close(F1); 
  }
  close(OJ1); 
}

# Clean temporary directory. 
&fileSunhh::_rmtree($wdir); 



##############################################
############## Subfunctions ##################
#
# Write unaligned sequence files, and go back to the next cycle. 
sub resolve_id {
  my ($oldID, $s, $e) = @_; 
  $oldID =~ m!^(\S+)__(\d+)_(\d+)$! or return("${oldID}__${s}_${e}");
  my ($basic, $os, $oe) = ($1, $2, $3); 
  my $news = $s + $os - 1; 
  my $newe = $e + $os - 1; 
  return("${basic}__${news}_${newe}"); 
}# resolve_id() 

sub runCmd {
  &exeCmd_1cmd($_[0]) and &stopErr("[Err] $_[0]\n"); 
}# runCmd() 


