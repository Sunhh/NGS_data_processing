my $id = shift;
my $file = shift;

open(FH, $file) || die $!;
while(<FH>)
{
	chomp;
	
	$s = <FH>;
	$i = <FH>;
	$q = <FH>;

	$_ =~ s/\s+.*//;
	
	if ($id eq $_)
	{
		print $_."\n".$s.$i.$q;
		exit;
	}

}
close(FH);
