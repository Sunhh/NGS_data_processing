#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts,
	"help!", 
	"colN:i", # 4 
	"bg_geneList:s", # all genes needed; 
); 
$opts{'colN'} //= 4; 

!@ARGV and die "perl $0 KEGG_PWY.txt -colN 4 -bg_geneList Cma_geneID_list > Cma_geneID_list.kgPWY_bg\n"; 

my %bg_genLis = %{ &load_bg_geneList( $opts{'bg_geneList'} ) }; 

my $fh_txt = &openFH( $ARGV[0], '<' ); 
my %pwy_annot; 
my @hh = &splitL("\t", &wantLineC($fh_txt)); 
while (&wantLineC($fh_txt)) {
	my @ta = &splitL("\t", $_); 
	$ta[$opts{'colN'}] eq '' and next; 
	my @gen = &splitL(",", $ta[$opts{'colN'}]); 
	my $annot_txt = "map$ta[0]!!$ta[3]!!$ta[1]"; 
	for my $gid (@gen) {
		$gid =~ s!^\s+|\s+$!!g; 
		defined $pwy_annot{$gid}{'annot_txt'}{$annot_txt} and next; 
		$pwy_annot{$gid}{'annot_cnt'}++; 
		$pwy_annot{$gid}{'annot_txt'}{$annot_txt} = $pwy_annot{$gid}{'annot_cnt'}; 
	}
}
close($fh_txt); 

for my $gid (@{$bg_genLis{'id_a'}}) {
	my $annot = ''; 
	defined $pwy_annot{$gid} and $annot = join(';;', sort { $pwy_annot{$gid}{'annot_txt'}{$a} <=> $pwy_annot{$gid}{'annot_txt'}{$b} } keys %{$pwy_annot{$gid}{'annot_txt'}}); 
	print STDOUT join("\t", $gid, $annot)."\n"; 
}


sub load_bg_geneList {
	my $fn = shift; 
	my $fh = &openFH( $fn, '<' ); 
	my %back; 
	while (&wantLineC($fh)) {
		my @ta=&splitL("\t", $_); 
		push(@{$back{'id_a'}}, $ta[0]); 
		$back{'id_h'}{$ta[0]} //= $.; 
	}
	close($fh); 
	return (\%back); 
}#



