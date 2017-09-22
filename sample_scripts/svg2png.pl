#!/usr/bin/perl -w 
use strict; 
use Cwd; 

!@ARGV and die "perl $0 fdsaf\n"; 

my $origin_dir = getcwd(); 

for my $dn (`ls -d *_data`) {
	chomp($dn); 
	print STDOUT "DIR=[$dn]\n"; 
	-d $dn or next; 
	chdir($dn); 
	my $rilN = $dn; 
	$rilN =~ s/_data$//; 
	
	mkdir("PNG", 0755); 
	mkdir("SVG", 0755); 
	my @png_files; 
	my @svg_files; 
	for my $fn (`ls *.svg`) {
		chomp($fn); 
		$fn =~ s/\.svg$//; 
		system "convert $fn.svg $fn.png"; 
		print "convert $fn.svg $fn.png\n"; 
		push(@png_files, "$fn.png"); 
		push(@svg_files, "$fn.svg"); 
	}
	my $nn = scalar(@png_files); 
	my $tt   = $" ; 
	local $" = " "; 
	my $merge_cmd = "montage @png_files -tile 1x$nn -geometry -0-0 ${rilN}_Chroms.png"; 
	my $mv_svg_cmd = "mv @svg_files SVG/"; 
	my $mv_png_cmd = "mv @png_files PNG/"; 
	$" = $tt; 
	system("$merge_cmd"); 
	print STDOUT "[Cmd]$merge_cmd\n"; 
	system("$mv_svg_cmd"); 
	print STDOUT "[Cmd]$mv_svg_cmd\n"; 
	system("$mv_png_cmd"); 
	print STDOUT "[Cmd]$mv_png_cmd\n"; 
	chdir($origin_dir); 
	print STDOUT "[Msg]$dn processing done.\n"; 
}

