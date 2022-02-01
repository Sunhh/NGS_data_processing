#!/usr/bin/perl
use strict;
use warnings;

-t and !@ARGV and die "perl $0 ipr.SearchResults-transposases.json > ipr.SearchResults-transposases.json.tbl\n";

print STDOUT join("\t", qw/id source source_database name description/)."\n";
while (<>) {
  m!^\s*(#|$)! and next;
  chomp;
  my @lines = &parse_iprJson($_);
  for my $a1 (@lines) {
    for my $k1 (keys %$a1) {
      $a1->{$k1} //= "NA";
    }
    print STDOUT join("\t", @{$a1}{qw/id source source_database name description/})."\n";
  }
}

sub parse_iprJson {
  my ($txt1) = @_;
  $txt1 =~ s!^\s*\[(.*)\]\s*$!$1! or die "Err 1:$txt1\n";
  my @back;
  while ($txt1 =~ s!^\s*\{ "id":"(\S+)", "source":"([^"]+)", "fields":\{ "description":\[(?:"(.*?)")?\], "name":\["([^"]+)"\], "source_database":\["([^"]+)"\] \} \s* \}\s*,*!!x) {
    my %h;
    @h{qw/id source description name source_database/} = ($1, $2, $3, $4, $5);
    push(@back, \%h);
  }
  $txt1 =~ m!^\s*$! or die "Err 2: |$txt1|\n";
  return(@back);
}# parse_iprJson()

# "id":"IPR005021",
#"source":"interpro7_family",
#"fields":{
#  "description":["&lt;p&gt;Terminase large subunit (TerL) from bacteriophages and evolutionarily related viruses, is an important component of the DNA packing machinery and comprises an ATPase domain, which powers DNA translocation and a nuclease domain that cuts concatemeric DNA [[cite:PUB00082566], [cite:PUB00098421], [cite:PUB00098422]]. TerL forms pentamers in which the ATPase domains form a ring distal to the capsid.  The ATPase domain contains a C-terminal subdomain that sits above the ATPase active site, called the \"Lid subdomain\" with reference to analogous lid subdomains found in other ATPases [[cite:PUB00098422]]. It contains a hydrophobic patch (Trp and Tyr residues) that mediates critical interactions in the interface between adjacent ATPase subunits and assists the positioning of the arginine finger residue that catalyses ATP hydrolysis [[cite:PUB00098422]]. The endonuclease cuts concatemeric DNA first in the initiation phase in a sequence specific site and later in the completion stage of the DNA packaging process when the capsid is full [[cite:PUB00098421], [cite:PUB00098422]]. Cryo-EM studies indicate that TerL forms a pentamer that binds to a dodecameric assembly called portal and attaches to the capsid. It has been proposed that nuclease domains form a radially arranged ring that is proximal to portal, playing a key role in pentamer assembly [[cite:PUB00098422]]. The nuclease domain has a RNAse H-like fold and it has been proposed to utilise a two-metal catalysis mechanism like in other RNAse H-like endonucleases such as RNase H, transposases, retroviral integrases and RuvC Holliday junction resolvases [[cite:PUB00098421]]. This entry also includes uncharacterised bacterial sequences.&lt;/p&gt;"
#  ],
#  "name":["Terminase large subunit, Lambdalikevirus-type"],
#  "source_database":["InterPro"]
#}

