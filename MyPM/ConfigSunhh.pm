package ConfigSunhh; 
#BEGIN {
#	push(@INC,'/usr/local/share/perl5/');
#}

use strict; 
use warnings; 
use LogInforSunhh; 

sub new {
	my $class = shift; 
	my $self = {}; 
	bless $self, $class; 
	
	$self->_initialize(@_);

	return $self; 
}# new() 

sub _initialize {
	my $self = shift; 
	my %parm = @_; 
	return; 
}

=head2 getConfig( 'cfg_file'=>'path.conf', 'replace'=> 1, 'hash_r'=> {} )

Function : 
  This is a method in object. 
  Read 'cfg_file' file and store key,val pairs in hash_reference 'hash_r'. 
  If 'replace' == TRUE, repeated value will be replaced by the newly input one. 

Return   : $self->{'hash_r'} 

 Sample file : 
# [Sunhh@wwz repAnno_tools]$ more path.conf
# dir2                    /home/Sunhh/tools/github/NGS_data_processing
# pl_deal_fasta           __dir2__/deal_fasta.pl
# pl_deal_table           __dir2__/deal_table.pl
#
# dir1                       /home/Sunhh/tools/github/NGS_data_processing/repAnno_tools
# pl_ch_gff_to_tab           __dir1__/ch_gff_to_tab.pl
# pl_ch_seqID                __dir1__/ch_seqID.pl
# pl_filter_tab_byPBSPPT     __dir1__/filter_tab_byPBSPPT.pl
# pl_name_from_tab           __dir1__/name_from_tab.pl
# pl_filter_flank            __dir1__/filter_flank.pl
# pl_filter_RepMsk_out       __dir1__/filter_RepMsk_out.pl
# pl_build_Examplar_byFa     __dir1__/build_Examplar_byFa.pl
# pl_lis_masked_RepMsk_out   __dir1__/lis_masked_RepMsk_out.pl
# pl_get_LTR_wi_Termi        __dir1__/get_LTR_wi_Termi.pl
#
# exe_RepeatMasker          /data/Sunhh/src/Annot/repeatmasker/RepeatMasker/RepeatMasker
# exe_gt                    /data/Sunhh/src/Annot/genometools/gt-1.5.3-complete/bin/gt
# exe_makeblastdb           /usr/local/bin/makeblastdb
# exe_blastn                /usr/local/bin/blastn
# exe_mv                    /usr/bin/mv
### Input  : ('cfg_file'=>'path.conf', 'replace'=>1, 'hash_r'=>\%tool)
### Return : \%tool 
=cut
sub getConfig {
	my $self = shift; 
	my %parm = @_; 
	$parm{'cfg_file'} = $parm{'cfg_file'} // $parm{'cfg'} // 'path.conf'; 
	$parm{'replace'}  = $parm{'replace'} // 1; 
	$parm{'hash_r'}   = $parm{'hash_r'} // {}; 

	open CF,'<',"$parm{'cfg_file'}" or &stopErr("[Err] file [$parm{'cfg_file'}] $!\n"); 
	while (<CF>) {
		chomp; 
		m/^\s*(#|$)/ and next; 
		s/[^\S\t]+$//; 
		$_ =~ s!^(\S+)\s+!!; 
		my $tk = $1; 
		my $tv = $_; 
		while ($tv =~ m/__([^\s_]+)__/) {
			my $pk = $1; 
			defined $parm{'hash_r'}{$pk} or &stopErr("[Err] Unknown key [$pk]\n"); 
			my $pv = $parm{'hash_r'}{$pk}; 
			$tv =~ s/__${pk}__/$pv/; 
		}
		$parm{'hash_r'}{$tk} //= $tv; 
		$parm{'replace'} and $parm{'hash_r'}{$tk} = $tv; 
		&tsmsg("[Msg] Setting $tk=$tv\n"); 
	}
	close CF; 

	return $parm{'hash_r'}; 
}# getConfig() 

=head2 writeConfig( 'cfg_file'=>'path.conf', 'hash_r'=> {} )

Function : 
  This is a method in object. 
  Write 'cfg_file' file and store key,val pairs in hash_reference 'hash_r'. 

Return   : ()

=cut
sub writeConfig {
	my $self = shift; 
	my %parm = @_; 
	open O,'>', "$parm{'cfg_file'}" or &stopErr("[Err] file [$parm{'cfg_file'}] $!\n");
	for my $k (sort keys %{$parm{'hash_r'}}) {
		print O join("\t", $k, $parm{'hash_r'}{$k})."\n"; 
	}
	close(O); 
}# writeConfig() 



1; 

