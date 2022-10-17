#!/usr/bin/perl
# [3/17/2022] Need to add the cumulative coverage length of the query gene.
# [3/22/2022] The output gff3 of liftoff needs to be filtered because there are still gene models with gene coverage lower than required.
#             Strand-specific overlapping has been fixed.
# [3/23/2022] I decide not to filter coverage in this program, but I'll output the coverage bases in the output intersection file.
# [10/11/2022] Duplicate input files of Liftoff to avoid using of a single file by multiple Liftoff processes.
use strict;
use warnings;
use LogInforSunhh;
use fileSunhh;
use Getopt::Long;
my %opts;
GetOptions(\%opts,
  "from_tag:s",   # 'CLpan'
  "from_gff3:s",  # "in.from.trim2CDS.gff3"
  "from_fas:s",   # "in.from.genome.fa"
  "to_tag:s",     # 'CMpan'
  "to_gff3:s",    # "in.to.trim2CDS.gff3"
  "to_fas:s",     # "in.to.genome.fa"
  "out_dir:s",       # "output"
  "para_minCov:f",   # 0.5
  "para_minIdent:f", # 0.75
  "para_liftoff:s",  # " -p 20 -s $opts{'para_minIdent'} -a $opts{'para_minCov'} -copies "
  "pl_cntR2Q:s",     # "perl /home/Sunhh/tools/github/NGS_data_processing/annot_tools/liftoff_tools/cnt_R2Q_liftoff_info.pl "
  "exe_liftoff:s",   # "liftoff"
  "help!"
);
my $htxt = <<HH;
####################################################################################################
# perl $0  \\
#   -from_tag CLpan  -from_gff3 CLpan.trim2CDS.gff3  -from_fas CLpan.genome.fa \\
#   -to_tag   CMpan  -to_gff3   CMpan.trim2CDS.gff3  -to_fas   CMpan.genome.fa \\
#   -out_dir  output
# 
# Optional parameters:
# -para_minCov    [0.5]
# -para_minIdent  [0.75]
# -para_liftoff   [ -p 20 -s {para_minIdent} -a {para_minCov} -copies ]
# -pl_cntR2Q      [perl /home/Sunhh/tools/github/NGS_data_processing/annot_tools/liftoff_tools/cnt_R2Q_liftoff_info.pl]
# -exe_liftoff    [liftoff]
####################################################################################################
HH

for my $k1 (qw/from_tag from_gff3 from_fas to_tag to_gff3 to_fas/) {
  defined $opts{$k1} or &LogInforSunhh::usage($htxt);
}
$opts{'out_dir'}       //= 'output';
$opts{'para_minCov'}   //= 0.5;
$opts{'para_minIdent'} //= 0.75;
$opts{'para_liftoff'}  //= " -p 20 -s $opts{'para_minIdent'} -a $opts{'para_minCov'} -copies ";
$opts{'pl_cntR2Q'}     //= "perl /home/Sunhh/tools/github/NGS_data_processing/annot_tools/liftoff_tools/cnt_R2Q_liftoff_info.pl";
$opts{'exe_liftoff'}   //= 'liftoff';

my %dta;
$dta{'tmp_dir'}  = &fileSunhh::new_tmp_dir('create' => 1);
for (qw/para_minCov para_minIdent para_liftoff pl_cntR2Q exe_liftoff/) {
  $dta{$_} = $opts{$_};
}
for (qw/out_dir from_tag from_gff3 from_fas to_tag to_gff3 to_fas/) {
  $dta{$_} = $opts{$_};
}
-e $dta{'out_dir'} or mkdir($dta{'out_dir'});
$dta{'out_gff3'}  = "$dta{'out_dir'}/mapCDS.$dta{'from_tag'}.to.$dta{'to_tag'}.gff3";
$dta{'out_tbl'}   = "$dta{'out_dir'}/mapCDS.$dta{'from_tag'}.to.$dta{'to_tag'}.tbl";
$dta{'out_unmap'} = "$dta{'out_dir'}/mapCDS.$dta{'from_tag'}.to.$dta{'to_tag'}.unmapped";
$dta{'toRM'} = [];

# Copy input files of Liftoff to 'tmp_dir' to avoid conflicts between different processes.
&fileSunhh::_copy($dta{'from_gff3'}, "$dta{'tmp_dir'}/from.gff3");
&fileSunhh::_copy($dta{'from_fas'}, "$dta{'tmp_dir'}/from.fas");
&fileSunhh::_copy($dta{'to_gff3'}, "$dta{'tmp_dir'}/to.gff3");
&fileSunhh::_copy($dta{'to_fas'}, "$dta{'tmp_dir'}/to.fas");

# &runCmd("$dta{'exe_liftoff'} $dta{'para_liftoff'} -dir $dta{'tmp_dir'} -g $dta{'from_gff3'} -o $dta{'out_gff3'} -u $dta{'out_unmap'} $dta{'to_fas'} $dta{'from_fas'}");
&runCmd("$dta{'exe_liftoff'} $dta{'para_liftoff'} -dir $dta{'tmp_dir'} -g $dta{'tmp_dir'}/from.gff3 -o $dta{'out_gff3'} -u $dta{'out_unmap'} $dta{'tmp_dir'}/to.fas $dta{'tmp_dir'}/from.fas");
# &runCmd("$dta{'pl_cntR2Q'}   $dta{'from_gff3'}  $dta{'to_gff3'}  $dta{'out_gff3'} > $dta{'out_tbl'}");
&runCmd("$dta{'pl_cntR2Q'}   $dta{'tmp_dir'}/from.gff3  $dta{'tmp_dir'}/to.gff3  $dta{'out_gff3'} > $dta{'out_tbl'}");
push(@{$dta{'toRM'}}, $dta{'tmp_dir'});

for (@{$dta{'toRM'}}) {
  &fileSunhh::_rmtree($_);
}

