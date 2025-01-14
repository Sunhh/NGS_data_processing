#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 

-t and !@ARGV and die "gzip -cd gene_ontology_edit.obo.2018-05-01.gz | perl $0 > gene_ontology_edit.obo.2018-05-01.tab\n"; 

my %gg; 
$gg{'nmspace2s'}{'biological_process'} = 'BP'; 
$gg{'nmspace2s'}{'molecular_function'} = 'MF'; 
$gg{'nmspace2s'}{'cellular_component'} = 'CC'; 
$gg{'nmspace2s'}{'external'}           = 'ET'; 

my %curr; 
print STDOUT join("\t", qw/GO_ID ALT_ID Type GO_name GO_def/)."\n"; 
while (<>) {
	chomp; 
	if (m!^\[Term\]! or m!^$!) {
		&out_curr(); 
		%curr = (); 
		next; 
	} elsif (m!^id:\s+(\S+)\s*$!) {
		&setV('id', $1); 
	} elsif (m!^name:\s+(.+?)$!) {
		&setV('name', $1); 
	} elsif (m!^namespace:\s+(\S+)$!) {
		&setV('namespace', $1); 
	} elsif (m!^alt_id:\s+(GO:\d+)$!) {
		&setV('alt_id', $1); 
	} elsif (m!^def:\s+(.+?)$!) {
		&setV('def', $1); 
	}
}
scalar(keys %curr) > 0 and &out_curr(); 
%curr=(); 

sub out_curr {
	scalar(keys %curr) == 0 and return; 
	for my $k1 (qw/alt_id/) {
		$curr{$k1} //= ['NA']; 
	}
	$curr{'alt_id'} = join(';', @{$curr{'alt_id'}}); 
	for my $k2 (qw/id name namespace/) {
		defined $curr{$k2} or &stopErr("[Err] Lack $k2 for $curr{'id'}\n"); 
	}
	for my $k3 (qw/def/) {
		$curr{$k3} //= 'NA'; 
	}
	$curr{'namespace'} = $gg{'nmspace2s'}{ $curr{'namespace'} }; 
	print STDOUT join("\t", @curr{qw/id alt_id namespace name def/})."\n"; 
	%curr = (); 
	return; 
}# out_curr() 
sub setV {
	my ($k, $v) = @_; 
	if ($k eq 'alt_id') {
		push(@{$curr{$k}}, $v); 
	} else {
		defined $curr{$k} and &stopErr("[Err] Repeat curr_$k : $k:$curr{$k}:$v\n"); 
		$curr{$k} = $v; 
	}
	return; 
}
