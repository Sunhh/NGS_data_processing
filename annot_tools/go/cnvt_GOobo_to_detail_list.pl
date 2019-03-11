#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 

-t and !@ARGV and die "gzip -cd gene_ontology_edit.obo.2018-05-01.gz | perl $0 > gene_ontology_edit.obo.2018-05-01.tab\n"; 

my %gg; 
$gg{'nmspace2s'}{'biological_process'} = 'BP'; 
$gg{'nmspace2s'}{'molecular_function'} = 'MF'; 
$gg{'nmspace2s'}{'cellular_component'} = 'CC'; 
$gg{'nmspace2s'}{'external'}           = 'ET'; 

my %id2info; 
my %curr; 
while (<>) {
	chomp; 
	if (m!^\[(Term|Typedef)\]! or m!^$!) {
		&save_curr(); 
		%curr = (); 
		next; 
	} elsif (m!^id:\s+(\S+)\s*$!) {
		&setV('id', $1); 
	} elsif (m!^alt_id:\s+(GO:\d+)$!) {
		&setV('alt_id', $1); 
	} elsif (m!^is_a:\s+(\S+)\s!) {
		&setV('is_a', $1); 
	} elsif (m!^replaced_by:\s+(\S+)\s*$!) {
		&setV('replaced_by', $1); 
	} elsif (m!^relationship:\s+(\S+)\s+(\S+)\s!) {
		# &setV("RS_", $2); 
		&setV("RS_$1", $2); 
	} elsif (m!^namespace:\s+(\S+)$!) {
		&setV('namespace', $1); 
	} elsif (m!^name:\s+(.+?)$!) {
		&setV('name', $1); 
	} elsif (m!^def:\s+(.+?)$!) {
		&setV('def', $1); 
	} elsif (m!^is_obsolete: (\S+)$!) {
		&setV('is_obsolete', $1); 
	} elsif (m!^(comment|synonym|xref|subset|subsetdef|transitive_over|remark|synonymtypedef|is_transitive|ontology|is_metadata_tag|is_class_level|holds_over_chain|format\-version|default\-namespace|data\-version|consider|):\s!) {
		next; 
	} else {
		&stopErr("[Err] Unknown line: $_\n"); 
	}
}
scalar(keys %curr) > 0 and &save_curr(); 
%curr=(); 
&fill_id2info(); 
&out_id2info(); 

sub out_id2info {
	# print STDOUT join("\t", qw/GO_ID Obsolete ALT_ID namespace GO_name GO_def MinR2N AncestorIDs AncestorChains/)."\n"; 
	# GO_ID   ALT_ID  Type    GO_name GO_def
	print STDOUT join("\t", qw/GO_ID ALT_ID namespace GO_name GO_def Obsolete MinR2N AncestorIDs AncestorChains/)."\n"; 
	my %hasOutID; 
	for my $id (sort keys %id2info) {
		my $repID = $id2info{$id}{'repres'}; 
		$id2info{$repID}{'namespace'} eq 'external' and next; 
		$repID eq '' and die "I1:$id\n"; 
		defined $hasOutID{$repID} and next; 
		my $is_obsolete = (defined $id2info{$repID}{'is_obsolete'}) ? "Y" : "N" ; 
		$id2info{$repID}{'otherID'} //= {}; 
		$id2info{$repID}{'ancestorIDs'} //= {}; 
		$id2info{$repID}{'ancestorChain'} //= []; 
		my $otxt_otherID = join(";", sort keys %{$id2info{$repID}{'otherID'}}); 
		my $otxt_nmspace = $gg{'nmspace2s'}{ $id2info{$repID}{'namespace'} }; 
		my $otxt_ancestorIDs = join(";", sort keys %{$id2info{$repID}{'ancestorIDs'}}); 
		my $otxt_ancestorChain = (scalar(@{$id2info{$repID}{'ancestorChain'}}) > 0) ? join(";;", map { join(";", @$_) } @{$id2info{$repID}{'ancestorChain'}}) : ""; 

		$otxt_otherID eq '' and $otxt_otherID = "NA"; 
		$otxt_nmspace eq '' and &stopErr("[Err] Bad namespace [$id2info{$repID}{'namespace'}]\n"); 
		$id2info{$repID}{'minChainLen'} //= -1; 
		$otxt_ancestorIDs eq '' and $otxt_ancestorIDs = "NA"; 
		$otxt_ancestorChain eq '' and $otxt_ancestorChain = "NA"; 
		print STDOUT join("\t", 
			$repID, 
			$otxt_otherID, 
			$otxt_nmspace, 
			$id2info{$repID}{'name'}, 
			$id2info{$repID}{'def'}, 
			$is_obsolete, 
			$id2info{$repID}{'minChainLen'}, 
			$otxt_ancestorIDs,
			$otxt_ancestorChain
		)."\n"; 
		$hasOutID{$repID} = 1; 
	}
	return; 
}# out_id2info() 

sub fill_id2info {
	# Check if the representative is wrong. 
	for my $id (keys %id2info) {
		defined $id2info{$id}{'repres'} or die "I2:$id\n"; 
		my $repID = $id2info{$id}{'repres'}; 
		$id ne $repID and $id2info{$repID}{'otherID'}{$id}++; 
		defined $id2info{$id}{'is_removed'} and next; 

		if (defined $id2info{$id}{'tag_repres'}) {
			# This is main Term listed in OBO. 
			if ($id2info{$id}{'tag_repres'} eq 'self') {
				$id2info{$id}{'repres'} eq $id or die "self:$id\n"; 
			} elsif ($id2info{$id}{'tag_repres'} eq 'replaced') {
				my $repID = $id2info{$id}{'repres'}; 
				$id2info{$repID}{'repres'} eq $repID or die "self:$id\n"; 
			} elsif ($id2info{$id}{'tag_repres'} eq 'removed') {
				; 
			}
		} else {
			# This is alt_id or some other id only shown in other ID's term. 
			my $repID = $id2info{$id}{'repres'}; 
			$id2info{$repID}{'repres'} eq $repID or die "new_self:$id:$repID\n"; 
			$id2info{$repID}{'tag_repres'} eq 'self' or die "sss:$id:$repID\n"; 
		}
	}
	# Get parent-child information and set ALT_ID from represIDs; 
	for my $id (keys %id2info) {
		defined $id2info{$id}{'is_removed'} and next; 
		my $repID = $id2info{$id}{'repres'}; 
		if (defined $id2info{$id}{'is_a'}) {
			for my $pID (@{$id2info{$id}{'is_a'}}) {
				my $p_repID = $id2info{$pID}{'repres'}; 
				$p_repID eq '' and die "bad pID:$pID\n"; 
				$id2info{$repID}{'good_parent'}{$p_repID} = 1; 
			}
		}
		if (defined $id2info{$id}{'RS_part_of'}) {
			for my $pID (@{$id2info{$id}{'RS_part_of'}}) {
				my $p_repID = $id2info{$pID}{'repres'}; 
				$p_repID eq '' and die "bad pID:$pID\n"; 
				$id2info{$repID}{'good_parent'}{$p_repID} = 1; 
			}
		}
	}
	# Find all ancestor terms for each repID_term; 
	for my $id (keys %id2info) {
		defined $id2info{$id}{'is_removed'} and next; 
		$id2info{$id}{'repres'} eq $id or next; 
		$id2info{$id}{'ancestorIDs'} = {}; 
		$id2info{$id}{'ancestorChain'} = [ &chain_root2node([$id]) ]; 
		my $minChainLen; 
		for my $a1 (@{$id2info{$id}{'ancestorChain'}}) {
			my $cl = scalar(@$a1); 
			$id2info{$id}{'minChainLen'} //= $cl; 
			$id2info{$id}{'minChainLen'} > $cl and $id2info{$id}{'minChainLen'} = $cl; 
			for my $a2 (@$a1) {
				$a2 eq $id and next; 
				$id2info{$id}{'ancestorIDs'}{$a2}++; 
			}
		}
		$id2info{$id}{'minChainLen'} //= -1; 
	}

	return; 
}# fill_id2info() 
sub chain_root2node {
	my $has_new = 0; 
	my @back; 
	for my $a1 (@_) {
		if ( defined $id2info{ $a1->[0] }{'is_removed'} ) {
			push(@back, [@$a1]); 
			next; 
		} elsif ( ! defined $id2info{ $a1->[0] }{'good_parent'} ) {
			push(@back, [@$a1]); 
			next; 
		}
		for my $pID (keys %{ $id2info{ $a1->[0] }{'good_parent'} }) {
			push(@back, [$pID, @$a1]); 
		}
		$has_new = 1; 
	}
	$has_new == 1 and @back = &chain_root2node(@back); 
	return(@back); 
}# chain_root2node() 

sub save_curr {
	scalar(keys %curr) == 0 and return; 
	# GO_ID  : $curr{'id'} 
	# ALT_ID : @{$curr{'alt_id'}} + inferred
	# Good_ID: $curr{'replaced_by'}
	# Parent_ID : @{$curr{'is_a'}} + @{$curr{'RS_*'}}
	# IsGood : $curr{'is_obsolete'}
	# Type   : $curr{'namespace'}
	# Name   : $curr{'name'}
	# Descr  : $curr{'def'}
	for my $k1 (keys %curr) {
		$k1 eq 'id' and next; 
		$id2info{$curr{'id'}}{$k1} = $curr{$k1}; 
	}
	if (defined $curr{'replaced_by'}) {
		defined $curr{'is_obsolete'} or &stopErr("[Err]id:$curr{'id'}\n"); 
		$id2info{$curr{'id'}}{'tag_repres'} = 'replaced'; 
		$id2info{$curr{'id'}}{'repres'} //= $curr{'replaced_by'}; 
		if (defined $curr{'alt_id'}) {
			for my $id (@{$curr{'alt_id'}}) {
				$id2info{$id}{'repres'} //= $curr{'replaced_by'}; 
			}
		}
	} elsif (defined $curr{'is_obsolete'}) {
		# This GO term is completely removed. 
		# Such as : GO:0001511; 
		# &stopErr("[Err] Bad is_obsolete: $curr{'id'}\n"); 
		$id2info{$curr{'id'}}{'tag_repres'} = 'removed'; 
		$id2info{$curr{'id'}}{'is_removed'} = 1; 
		$id2info{$curr{'id'}}{'repres'}     = $curr{'id'}; 
		if (defined $curr{'alt_id'}) {
			for my $id (@{$curr{'alt_id'}}) {
				$id2info{$id}{'is_removed'} = 1; 
				$id2info{$id}{'repres'}     //= $curr{'id'}; 
			}
		}
	} elsif (defined $curr{'alt_id'}) {
		$id2info{$curr{'id'}}{'repres'} //= $curr{'id'}; 
		$id2info{$curr{'id'}}{'tag_repres'} = 'self'; 
		for my $id (@{$curr{'alt_id'}}) {
			$id2info{$id}{'repres'} //= $curr{'id'}; 
		}
	} else {
		# There is only one 'id'; 
		$id2info{$curr{'id'}}{'repres'} //= $curr{'id'}; 
		$id2info{$curr{'id'}}{'tag_repres'} = 'self'; 
	}

	return; 
}# save_curr() 
sub setV {
	my ($k, $v) = @_; 
	if ($k =~ m!^(alt_id|is_a|RS_(part_of|regulates|positively_regulates|negatively_regulates|))$!) {
		push(@{$curr{$k}}, $v); 
	} else {
		defined $curr{$k} and &stopErr("[Err] Repeat curr_$k : $k:$curr{$k}:$v\n"); 
		$curr{$k} = $v; 
	}
	return; 
}

