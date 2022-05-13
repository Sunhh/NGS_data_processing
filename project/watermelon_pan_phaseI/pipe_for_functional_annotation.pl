#!/usr/bin/perl
use strict;
use warnings;
use fileSunhh;
use LogInforSunhh;
use ConfigSunhh;
my $cs_obj = ConfigSunhh->new();


!@ARGV and die "perl $0 annot.cfg\n";
my $fn_cfg = shift;

my %par;
$cs_obj->getConfig( 'cfg_file' => $fn_cfg, 'replace' => 1, 'hash_r'=> \%par );

# blast2Nt
-e $par{'dir_bl2Nt'} or mkdir($par{'dir_bl2Nt'});
&runCmd("$par{'exe_bl2Nt'} $par{'para_bn2Nt'} -db $par{'db_bn2Nt'} -query $par{'in_cdsFaFn'} -out $par{'dir_bl2Nt'}/$par{'out_pref'}.c2Nt.bn6");

# blast2Nr
-e $par{'dir_bl2Nr'} or mkdir($par{'dir_bl2Nr'});
&runCmd("$par{'exe_bp2Nr'} $par{'para_bp2Nr'} -d $par{'db_bl2Nr'} -q $par{'in_pepFaFn'} -o $par{'dir_bl2Nr'}/$par{'out_pref'}.p2Nr.bp6");
&runCmd("$par{'exe_bx2Nr'} $par{'para_bx2Nr'} -d $par{'db_bl2Nr'} -q $par{'in_cdsFaFn'} -o $par{'dir_bl2Nr'}/$par{'out_pref'}.c2Nr.bp6");
### Run blast for blast2go
&runCmd("$par{'exe_bp2Nr'} $par{'para_bp2NrB2G'} -d $par{'db_bl2Nr'} -q $par{'in_pepFaFn'} -o $par{'dir_bl2Nr'}/$par{'out_pref'}.p2Nr.xml");

# blast2sprot
-e $par{'dir_bl2sprot'} or mkdir($par{'dir_bl2sprot'});
&runCmd("$par{'exe_bp2sprot'} $par{'para_bp2sprot'} -d $par{'db_bl2sprot'} -q $par{'in_pepFaFn'} -o $par{'dir_bl2sprot'}/$par{'out_pref'}.p2sprot.bp6");
&runCmd("$par{'exe_bx2sprot'} $par{'para_bx2sprot'} -d $par{'db_bl2sprot'} -q $par{'in_cdsFaFn'} -o $par{'dir_bl2sprot'}/$par{'out_pref'}.c2sprot.bp6");

# blast2trembl
-e $par{'dir_bl2trembl'} or mkdir($par{'dir_bl2trembl'});
&runCmd("$par{'exe_bp2trembl'} $par{'para_bp2trembl'} -d $par{'db_bl2trembl'} -q $par{'in_pepFaFn'} -o $par{'dir_bl2trembl'}/$par{'out_pref'}.p2trembl.bp6");
&runCmd("$par{'exe_bx2trembl'} $par{'para_bx2trembl'} -d $par{'db_bl2trembl'} -q $par{'in_cdsFaFn'} -o $par{'dir_bl2trembl'}/$par{'out_pref'}.c2trembl.bp6");

# blast2tair10
-e $par{'dir_bl2tair10'} or mkdir($par{'dir_bl2tair10'});
&runCmd("$par{'exe_bp2tair10'} $par{'para_bp2tair10'} -d $par{'db_bl2tair10'} -q $par{'in_pepFaFn'} -o $par{'dir_bl2tair10'}/$par{'out_pref'}.p2tair10.bp6");
&runCmd("$par{'exe_bx2tair10'} $par{'para_bx2tair10'} -d $par{'db_bl2tair10'} -q $par{'in_cdsFaFn'} -o $par{'dir_bl2tair10'}/$par{'out_pref'}.c2tair10.bp6");


