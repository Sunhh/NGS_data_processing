#!/usr/bin/perl
# [5/27/2022] Use parameter to turn off some alignments.
use strict;
use warnings;
use fileSunhh;
use LogInforSunhh;
use ConfigSunhh;
my $cs_obj = ConfigSunhh->new();
use Getopt::Long;
my %opts;
GetOptions(\%opts,
  "onlyAHRD!",
  "onlyB2G!",
  "onlyNt!",
  "notRun:s", #
  "noCds!",
  "help!",
);

my $htxt = <<HHH;
####################################################################################################
# perl $0 annot.cfg
#
#   -onlyAHRD / -onlyB2G  / -onlyNt / -onlyNr
#
#   -notRun    [] Could be: ahrd,bl2go,x2Nr,c2Nt,p2NrBp6,c2Nr,p2NrB2G
#   -noCds     [Boolean] Turn off cds query alignments.
####################################################################################################
HHH

my %todo = qw(
ahrd   1
bl2go  1
toNr   1
toNt   1
);

my %notRun;
{
  my $cnt = 0;
  for (qw/onlyAHRD onlyB2G onlyNt/) {
    defined $opts{$_} and $cnt ++;
  }
  $cnt > 1 and &stopErr("[Err] Only less than one is allowed: onlyAHRD onlyB2G onlyNt\n");
  defined $opts{'onlyAHRD'} and do { %todo = qw(ahrd 1); };
  defined $opts{'onlyB2G'}  and do { %todo = qw(bl2go 1); };
  defined $opts{'onlyNr'}   and do { %todo = qw(toNr 1); };
  defined $opts{'onlyNt'}   and do { %todo = qw(toNt 1); };
  $opts{'notRun'} //= '';
  for (split(/,/, $opts{'notRun'})) {
    s!\s!!g;
    $_ eq '' and next;
    $notRun{$_} = 1;
  }
}

!@ARGV and &LogInforSunhh::usage($htxt);

my $fn_cfg = shift;

my %par;
$cs_obj->getConfig( 'cfg_file' => $fn_cfg, 'replace' => 1, 'hash_r'=> \%par );

# Section: Three blast for AHRD annotation.
if ($todo{'ahrd'} and !(defined $notRun{'ahrd'})) {
  # blast2sprot
  -e $par{'dir_bl2sprot'} or mkdir($par{'dir_bl2sprot'});
  &runCmd("$par{'exe_bp2sprot'} $par{'para_bp2sprot'} -d $par{'db_bl2sprot'} -q $par{'in_pepFaFn'} -o $par{'dir_bl2sprot'}/$par{'out_pref'}.p2sprot.bp6");
  $opts{'noCds'} or &runCmd("$par{'exe_bx2sprot'} $par{'para_bx2sprot'} -d $par{'db_bl2sprot'} -q $par{'in_cdsFaFn'} -o $par{'dir_bl2sprot'}/$par{'out_pref'}.c2sprot.bp6");
  
  # blast2trembl
  -e $par{'dir_bl2trembl'} or mkdir($par{'dir_bl2trembl'});
  &runCmd("$par{'exe_bp2trembl'} $par{'para_bp2trembl'} -d $par{'db_bl2trembl'} -q $par{'in_pepFaFn'} -o $par{'dir_bl2trembl'}/$par{'out_pref'}.p2trembl.bp6");
  $opts{'noCds'} or &runCmd("$par{'exe_bx2trembl'} $par{'para_bx2trembl'} -d $par{'db_bl2trembl'} -q $par{'in_cdsFaFn'} -o $par{'dir_bl2trembl'}/$par{'out_pref'}.c2trembl.bp6");
  
  # blast2tair10
  -e $par{'dir_bl2tair10'} or mkdir($par{'dir_bl2tair10'});
  &runCmd("$par{'exe_bp2tair10'} $par{'para_bp2tair10'} -d $par{'db_bl2tair10'} -q $par{'in_pepFaFn'} -o $par{'dir_bl2tair10'}/$par{'out_pref'}.p2tair10.bp6");
  $opts{'noCds'} or &runCmd("$par{'exe_bx2tair10'} $par{'para_bx2tair10'} -d $par{'db_bl2tair10'} -q $par{'in_cdsFaFn'} -o $par{'dir_bl2tair10'}/$par{'out_pref'}.c2tair10.bp6");
}

# Section: Blast2Nr for blast2go
# blast2Nr
if ($todo{'toNr'} and !(defined $notRun{'x2Nr'})) {
  -e $par{'dir_bl2Nr'} or mkdir($par{'dir_bl2Nr'});
  defined $notRun{'p2NrBp6'} or &runCmd("$par{'exe_bp2Nr'} $par{'para_bp2Nr'} -d $par{'db_bl2Nr'} -q $par{'in_pepFaFn'} -o $par{'dir_bl2Nr'}/$par{'out_pref'}.p2Nr.bp6");
  $opts{'noCds'} or defined $notRun{'c2Nr'} or &runCmd("$par{'exe_bx2Nr'} $par{'para_bx2Nr'} -d $par{'db_bl2Nr'} -q $par{'in_cdsFaFn'} -o $par{'dir_bl2Nr'}/$par{'out_pref'}.c2Nr.bp6");
}
if ($todo{'bl2go'} and !(defined $notRun{'bl2go'})) {
  ### Run blast for blast2go
  -e $par{'dir_bl2Nr'} or mkdir($par{'dir_bl2Nr'});
  defined $notRun{'p2NrB2G'} or &runCmd("$par{'exe_bp2Nr'} $par{'para_bp2NrB2G'} -d $par{'db_bl2Nr'} -q $par{'in_pepFaFn'} -o $par{'dir_bl2Nr'}/$par{'out_pref'}.p2Nr.xml");
}

# Section: CDS blast for view.
if ($todo{'toNt'} and !(defined $notRun{'c2Nt'})) {
  # blast2Nt: This is too slow and not an urgent need, so I decide to do it at last.
  -e $par{'dir_bl2Nt'} or mkdir($par{'dir_bl2Nt'});
  $opts{'noCds'} or &runCmd("$par{'exe_bl2Nt'} $par{'para_bn2Nt'} -db $par{'db_bn2Nt'} -query $par{'in_cdsFaFn'} -out $par{'dir_bl2Nt'}/$par{'out_pref'}.c2Nt.bn6");
}

