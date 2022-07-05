#!/usr/bin/perl
use strict; 
use warnings; 
use fileSunhh; 
use LogInforSunhh; 

my $help_txt = <<HH; 
######################################
perl $0 outGATK_filtSNPs_PASS.vcf > outGATK_filtSNPs_PASS.vcf.depBySampleBySite

In outGATK_filtSNPs_PASS.vcf: The FORMAT should be 'GT:AD:DP...' for DP; 

HH

-t and !@ARGV and &LogInforSunhh::usage( $help_txt ); 


our @InFp = () ;
if ( !@ARGV )
{
	-t or @InFp = (\*STDIN);
}
else
{
	for (@ARGV) {
		push( @InFp, &openFH($_,'<') );
	}
}

my $has_header = 0; 
my @ha; 
for my $ifh (@InFp) {
	while (<$ifh>) {
		if (m!^##!) {
			# chomp; 
			# print STDOUT "$_\n"; 
			next; 
		}
		if (m!^#CHROM\t!) {
			chomp; 
			my @ta = split(/\t/, $_); 
			if ($has_header == 1) {
				my @tb = @ta[0,1,9..$#ta]; 
				scalar(@ha) == scalar(@tb) or &stopErr("[Err] Different header in VCF files.\n"); 
				for (my $i=0; $i<@ha; $i++) {
					$ha[$i] eq $tb[$i] or &stopErr("[Err] Different header in VCF files.\nha:@ha\ntb:@tb\n");
				}
				next; 
			}
			print STDOUT join("\t", @ta[0,1,9..$#ta])."\n"; 
			@ha = @ta[0,1,9..$#ta]; 
			$has_header = 1; 
			next; 
		}
		chomp; 
		my @ta = split(/\t/, $_); 
		my @oa = @ta[0,1]; 
		for (my $i=9; $i<@ta; $i++) {
			my @tb = split(/:/, $ta[$i]); 
			push(@oa, $tb[2]); 
		}
		print STDOUT join("\t", @oa)."\n"; 
	}
}




