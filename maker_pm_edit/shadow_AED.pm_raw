#------------------------------------------------------------------------
#----                          shadow_AED.pm                         ---- 
#------------------------------------------------------------------------
package shadow_AED;
use strict;
use Carp;
use Bit::Vector;

sub get_abAED{
    my $hits = shift;
    my $tran = shift;

    my $sum = 0;
    $tran->{_touch}++; #keeps self from being counted in score
    foreach my $h (@$hits){
	next if($tran->{_touch}++);
	$sum += get_AED([$h], $tran);
    }
    
    foreach my $h (@$hits, $tran) {
	delete($h->{_touch});
    }

    return 1 if(! @$hits);
    return $sum/@$hits;
}

sub get_eAED {
   my $hits = shift;
   my $tran = shift;
   my $seq = shift;

   return 1 if(! @{$hits} || ! $tran);

   #seperate out hit types
   my @prots;
   my @ests;
   my @others;

   foreach my $h (@$hits){
       if($h->algorithm =~ /protein2genome|^protein_gff|^blastx/i){
	   push(@prots, $h);
       }
       elsif($h->algorithm =~ /est2genome|^est_gff|^altest_gff|^blastn|^tblastx/i){
   	   push(@ests, $h);
       }
       else{
	   push(@others, $h);
       }
   }

   #get coordinates to build array in memory
   my ($start, $end) = ($tran->start('query'), $tran->end('query'));
   foreach my $hit (@{$hits}, $tran){
       my ($hs, $he) = ($hit->start('query'), $hit->end('query'));

       $start = $hs if($hs < $start);
       $end = $he if($he > $end);
   }

   #build array in memory
   my $length = $end - $start + 1;
   my $offset = $start; # do not - 1 so as to give array space coors
   my @b_seq = map {0} (1..$length); 

   #map out hit space
   foreach my $hit (@{$hits}){
       my @hsps = $hit->hsps();

       foreach my $hsp (@hsps){
	   my $s = $hsp->start('query') - $offset;
	   my $e = $hsp->end('query') - $offset;
	   
	   #array space coors
	   confess "ERROR: Start value not permited!!\n" if($s >= $length || $s < 0);
	   confess "ERROR: End value not permited!!\n" if($e < 0 || $e >= $length);
	   
	   @b_seq[$s..$e] = map {1} ($s..$e);
       }
   }

   #fill in space between splice site crossing reads
   if($seq){
       foreach my $set (@{PhatHit_utils::splice_infer_exon_coors(\@ests, $seq)}){
	   #array space coors
	   my $s = $set->[0] - $offset;
	   my $e = $set->[1] - $offset;
	   
	   confess "ERROR: Start value not permited!!\n" if($s >= $length || $s < 0);
	   confess "ERROR: End value not permited!!\n" if($e < 0 || $e >= $length);
	   
	   @b_seq[$s..$e] = map {1} ($s..$e);
       }
   }

   #==calculate bp in evidence
   my %index = (0 => 0, #empty
		1 => 0, #all evidence
		2 => 0, #all trans
		3 => 0, #overlap
	       );

   foreach my $i (@b_seq){
       $index{$i}++ if($i == 1);
   }

   #trans hsps
   my @hsps = sort {$a->start <=> $b->start} $tran->hsps();
   
   #==do EST filtering
   my @keepers;

   #build splice site index for transcript
   my %splices;
   for( my $i = 0; $i < @hsps; $i++){
       my $hsp = $hsps[$i];
       $splices{start}{$hsp->start} = 1 unless($i == 0);
       $splices{end}{$hsp->end} = 1 unless($i == @hsps - 1);
   }      

   #filter out EST hsps that don't have a matching splice site
   foreach my $e (@ests, @others){
       my @ehsps = sort {$a->start <=> $b->start} $e->hsps;
       
       EHSP: for (my $i = 0; $i < @ehsps; $i++){
	   my $ehsp = $ehsps[$i];
	   my $aB = $ehsp->start;
	   my $aE = $ehsp->end;
	   my $aL = abs($aE - $aB) + 1;
	   
	   #splice site matches so keep EST HSP
	   if($i != 0 && $splices{start}{$aB}){ #a first splice site in EST also first splice site in gene
	       $splices{start}{$aB} = 2 if($e->algorithm =~ /est2genome|^est_gff|^altest_gff|^blastn|^tblastx/i); #flag for verified
	       push(@keepers, $ehsp);
	       next EHSP;
	   }
	   elsif($i != @hsps - 1 && $splices{end}{$aE}){ #a last splice site in EST also last splice site in gene
	       $splices{end}{$aE} = 2 if($e->algorithm =~ /est2genome|^est_gff|^altest_gff|^blastn|^tblastx/i); #flag for verified
	       push(@keepers, $ehsp);
	       next EHSP;
	   }
	   
	   #no splice site match so check overlap of EST HSP with exon
	   for(my $j = 0; $j < @hsps; $j++){
	       my $hsp = $hsps[$j]; 
	       my $bB = $hsp->start;
	       my $bE = $hsp->end;
	       my $bL = abs($bE - $bB) + 1;
	       
	       my $class = compare::compare($aB, $aE, $bB, $bE);
	       next if($class eq '0');
	       
	       my $oB = ($aB > $bB) ? $aB : $bB; #overlap begin
	       my $oE = ($aE < $bE) ? $aE : $bE; #overlap end
	       my $oL = abs($oE - $oB) + 1;
	       
	       if(@hsps == 1 && @ehsps == 1 && $oL/$bL >= .9){ #at least 90% of single exon gene overlapped by single HSP EST
		   push(@keepers, $ehsp);
		   next EHSP;
	       }
	       
	       if($oL/$bL >= .9 && $oL/$aL >= .9){ #at least 90% of exon and 90% of single HSP EST overlap
		   push(@keepers, $ehsp);
		   next EHSP;
	       }
	   }
       }
   }

   #confirm CDS by filling in middle between confirmed splice sites
   #this fixes a bug where long exons are punished because only
   #the extremes overlap evidence (splice site crossing reads).
   #so the center of the exon is punished for specificity.
   foreach my $hsp (@hsps){
       #both splice sites are confirmed, so consider everything in the middle as verified
       if($splices{start}{$hsp->start('query')} == 2 && $splices{end}{$hsp->end('query')} == 2){
	   my $s = $hsp->start('query') - $offset;
	   my $e = $hsp->end('query') - $offset;
	   
	   #array space coors
	   confess "ERROR: Start value not permited!!\n" if($s >= $length || $s < 0);
	   confess "ERROR: End value not permited!!\n" if($e < 0 || $e >= $length);

	   foreach my $i ($s..$e){
	       $index{1}++ if($b_seq[$i] == 0); #fix evidence count
	       $b_seq[$i] = 3; #make as verified overlap
	   }
       }
   }

   #map out transcript space
   #done after CDS bug fix so as to properly acount for missing evidence
   foreach my $hsp (@hsps){
      my $s = $hsp->start('query') - $offset;
      my $e = $hsp->end('query') - $offset;
      
      #array space coors
      confess "ERROR: Start value not permited!!\n" if($s >= $length || $s < 0);
      confess "ERROR: End value not permited!!\n" if($e < 0 || $e >= $length);

      foreach my $i ($s..$e){
	  $b_seq[$i] = 2 if($b_seq[$i] != 3); #make transcript (replaces evidence)
      }
   }

   #==calculate bp in hit
   foreach my $i (@b_seq){
       $index{2}++ if($i == 2 || $i == 3); #must account for 3 because of CDS overlap coorrection
   }

   #map keeper EST HSPs
   foreach my $hsp (@keepers){
      my $s = $hsp->start('query') - $offset;
      my $e = $hsp->end('query') - $offset;
      
      #array space coors
      confess "ERROR: Start value not permited!!\n" if($s >= $length || $s < 0);
      confess "ERROR: End value not permited!!\n" if($e < 0 || $e >= $length);

      #only add to transcript overlap
      foreach my $i ($s..$e){
	  $b_seq[$i] += 1 if($b_seq[$i] == 2);
      }
   }

   #==map protein hits by phase
   if(@prots){
       my @ok_frames;
       foreach my $p (@prots){
	   foreach my $phsp ($p->hsps){
	       my $start = $phsp->start('query') - $offset;
	       my $end = $phsp->end('query') - $offset;
	       
	       my $pos = $start; #array position
	       my $cigar = $phsp->cigar_string();
	       if(! $cigar){ #if no gap attribute than we assume translation begins at first bp
		   my $length = abs($end - $start) + 1;
		   $cigar .= 'M'.int($length/3);
	       }
	       
	       my @gap = $cigar =~ /([A-Z]\d+)/g;
	       foreach my $g (@gap){
		   $g =~ /([A-Z])(\d+)/;
		   if($1 eq 'F'){
		       $pos += $2;
		   }
		   elsif($1 eq 'R'){
		       $pos -= $2;
		   }
		   elsif($1 eq 'D'){
		       $pos += ($2 * 3);
		   }
		   elsif($1 eq 'M'){
		       my $go = $2;
		       while($go--){
			   if($p->strand('query') == 1){
			       $ok_frames[$pos+0]->{0}++;
			       $ok_frames[$pos+1]->{1}++;
			       $ok_frames[$pos+2]->{2}++;
			   }
			   else{
			       $ok_frames[$pos+2]->{0}++;
			       $ok_frames[$pos+1]->{1}++;
			       $ok_frames[$pos+0]->{2}++;
			   }
			   $pos += 3;
		       }
		   }
	       }
	   }
       }
       
       #compare phase to transcript
       my $aB = $tran->{_TSTART}{query} - $offset;
       my $aE = $tran->{_TEND}{query} - $offset;
       my $phase = 0;
       ($aB, $aE) = ($aE, $aB) if($aB > $aE);
       foreach my $hsp (@{PhatHit_utils::sort_hits($tran, 'query')}){
	   my $bB = $hsp->start('query') - $offset;
	   my $bE = $hsp->end('query') - $offset;
	   my $class = compare::compare($aB, $aE, $bB, $bE);
	   next if ($class eq '0');
	   
	   $bB = $aB if($aB > $bB);
	   $bE = $aE if($aE < $bE);
	   
	   
	   my @select =  ($tran->strand('query') == 1) ? ($bB..$bE) : reverse($bB..$bE);
	   foreach my $i (@select){
	       $b_seq[$i] += 1 if($b_seq[$i] == 2 && $ok_frames[$i]->{$phase});
	       $phase = ($phase + 1) % 3;
	   }
       }
   }

   #==calculate overlap
   foreach my $i (@b_seq){
       $index{3}++ if($i == 3);
   }

   my @error = grep {$_ > 3} keys %index; #should not have keys greater than 3

   #catch error caused by bad GFF3 input (i.e. hits with no HSPs)
   die "ERROR: The feature being compared appears to be missing\n".
       "some of it's structure.  This can happen when you use\n".
       "a malformed GFF3 file as input to one of MAKER's evidence\n".
       "passthrough options. Failed on ". $tran->name."\n".
       "(from shadow_AED::get_eAED)\n"
       if($index{1} == 0 || $index{2} == 0 || @error);

   my $spec = $index{3}/$index{2}; #specificity
   my $sens = $index{3}/$index{1}; #sensitivity
   my $eAED = 1 - ($spec + $sens)/2;

   return $eAED;
}

sub get_AED {
   my ($sn, $sp) = get_SN_SP(@_);
   my $AED = 1 - ($sn + $sp)/2;

   return $AED;
}

sub get_SN_SP {
   my $hits = shift;
   my $tran = shift;

   return (0, 0) if(! @{$hits} || ! $tran);

   my ($start, $end) = ($tran->start('query'), $tran->end('query'));

   foreach my $hit (@{$hits}){
      my ($hs, $he) = ($hit->start('query'), $hit->end('query'));

      $start = $hs if($hs < $start);
      $end = $he if($he > $end);
   }

   my $length = $end - $start + 1;
   my $offset = $start; # do not - 1 so as to give array space coors
   my $t_vec = Bit::Vector->new($length);
   my $h_vec = Bit::Vector->new($length);
   
   #map out hit space
   foreach my $hit (@{$hits}){
      my @hsps = $hit->hsps() if defined $hit->hsps();      
      
      foreach my $hsp (@hsps){
	 my $s = $hsp->start('query') - $offset;
	 my $e = $hsp->end('query') - $offset;

	 #array space coors
	 confess "ERROR: Start value not permited!!\n" if($s >= $length || $s < 0);
	 confess "ERROR: End value not permited!!\n" if($e < 0 || $e >= $length);

	 $h_vec->Interval_Fill($s, $e);
      }
   }

   #map out transcript space
   my @hsps = $tran->hsps() if defined $tran->hsps();      
   
   foreach my $hsp (@hsps){
      my $s = $hsp->start('query') - $offset;
      my $e = $hsp->end('query') - $offset;
      
      #array space coors
      confess "ERROR: Start value not permited!!\n" if($s >= $length || $s < 0);
      confess "ERROR: End value not permited!!\n" if($e < 0 || $e >= $length);
      
      $t_vec->Interval_Fill($s, $e);
   }
   
   #calculate AED
   my $inter = Bit::Vector->new($length);
   $inter->Intersection($t_vec, $h_vec);
   my $sens = $inter->Norm()/$h_vec->Norm(); #sensitivity
   my $spec = $inter->Norm()/$t_vec->Norm(); #specificity

   return ($sens, $spec);
}

1;
