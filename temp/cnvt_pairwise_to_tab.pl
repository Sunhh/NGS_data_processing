#!/usr/bin/perl
use strict; 
use warnings; 
use Bio::SearchIO; 

!@ARGV and die "perl $0 in.blast > out.tab\n"; 

my $inBlastFile = shift; 

my $in = new Bio::SearchIO(
  -format => 'blast', 
  -file   => $inBlastFile
); 

while( my $result = $in->next_result ) {
  ## $result is a Bio::Search::Result::ResultI compliant object
  while( my $hit = $result->next_hit ) {
    ## $hit is a Bio::Search::Hit::HitI compliant object
    while( my $hsp = $hit->next_hsp ) {
      ## $hsp is a Bio::Search::HSP::HSPI compliant object
      print STDOUT join("\t", $result->query_name(), $hit->name(), $hit->description(), $hsp->score(), $hsp->evalue())."\n"; 
#      if( $hsp->length('total') > 50 ) {
#        if ( $hsp->percent_identity >= 75 ) {
#          print "Query=",   $result->query_name,
#            " Hit=",        $hit->name,
#            " Length=",     $hsp->length('total'),
#            " Percent_id=", $hsp->percent_identity, "\n";
#        }
#      }
    }  
  }
}

