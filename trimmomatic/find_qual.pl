
while(<>)
{
	$s = <>;
	$i = <>;
	$q = <>;
	
	for($i=0; $i<length($s); $i++)
	{
		$b = substr($s, $i, 1);
		$d = substr($q, $i, 1);

		print $i,"\t",$b,"\t$d\t",ord($d)-33 ,"\n";
	}
}
