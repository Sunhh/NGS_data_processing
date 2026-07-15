#!/usr/bin/perl -w
use mathSunhh;

!@ARGV and die "perl $0 background_IPR subset_list\n"; 

my $back=shift;
my $list=shift;
my $out_1= $list."_pvalue";
my $out_2= $list."_adjp";

open(B,$back);
open (L,$list);

open (O,">$out_1");
my %sub;      ### sublist
my %back_list; ### all_list
#my %des;      ### IPR description
my %ipr;      ### each ipr included bachground gene 
my %sub_ipr;  ### each ipr included sub list gene

while (<L>)
{
	chomp; 
	m/^\s*(#|$)/ and next; 
	my @a = split (/\t/,$_);
	$sub{$a[0]}=1;
}
	my $sub_num = scalar keys %sub;
	
while (<B>)
{
	chomp;
	# m/^\s*(#|$)/ and next; 
	my @b = split (/\t/,$_);
	my $gene = $b[0];
	$back_list{$gene}=1;
	defined $b[1] or next; 
	next if ($b[1] eq 'NA' or $b[1] eq '');
	#$b[1]=~s/;/\t/g; 
	$b[1] =~ s/,\s(IPR\d+)/\t$1/g; 
	my @tmp = split (/\t/,$b[1]);
	foreach $t (@tmp)
	{
		push @{$ipr{$t}},$gene;  ##########  1 
		if (exists $sub{$gene})	
		{
			push @{$sub_ipr{$t}},$gene;   #### 2
		}
	}
}
	my $back_num = scalar keys %back_list;
	
	my @pvals;   # p-values in out_1 row order, for BH adjustment below
	foreach $ip ( sort keys %sub_ipr)
	{
		my $sub_gene;
		my $pvalue;
		$back_list_num = scalar @{$ipr{$ip}};     ## 1
#		if (exists $sub_ipr{$ip})
#		{
				$sub_list_num = scalar @{$sub_ipr{$ip}};  ## 2
				$sub_gene= join (";",@{$sub_ipr{$ip}});
				$pvalue = hypergeometric($back_num,$back_list_num,$sub_num,$sub_list_num);
				print O "$ip\t$back_num\t$back_list_num\t$sub_num\t$sub_list_num\t$sub_gene\t$pvalue\n";
				push @pvals, $pvalue;
#		}
	}

	close (O);
	# Benjamini-Hochberg FDR adjustment of the p-values (col 7 of $out_1), in the
	# same row order. Pure-Perl replacement for R p.adjust(...,method="BH").
	my $adjp = &mathSunhh::p_adjust_BH(\@pvals);
	open (A,">$out_2") or die "[Err] Failed to write $out_2: $!\n";
	print A "$_\n" for @$adjp;
	close (A);



sub hypergeometric {
    my $n = $_[0];  #Total number of genes in all the pathways
    my $np = $_[1]; #Total number of genes in a particular pathway
    my $k = $_[2];  #Total number of genes in the input list from all the pathways
    my $r = $_[3];  #total number of genes in the input list from the particular pathway
    my $nq;
    my $top;

    $nq = $n - $np;

    $log_n_choose_k = &lNchooseK( $n, $k );

    $top = $k;
    if ( $np < $k ) {
        $top = $np;
    }

    $lfoo = &lNchooseK($np, $top) + &lNchooseK($nq, $k-$top);
    $sum = 0;

    for ($i = $top; $i >= $r; $i-- ) {
        $sum = $sum + exp($lfoo - $log_n_choose_k);

        if ( $i > $r) {
            $lfoo = $lfoo + log($i / ($np-$i+1)) +  log( ($nq - $k + $i) / ($k-$i+1)  )  ;
        }
    }
    return $sum;
}

sub lNchooseK {
    my $n = $_[0];
    my $k = $_[1];
    my $answer = 0;

    if( $k > ($n-$k) ){
        $k = ($n-$k);
    }

    for( $i=$n; $i>($n-$k); $i-- ) {
        $answer = $answer + log($i);
    }

    $answer = $answer - &lFactorial($k);
    return $answer;
}

sub lFactorial {
    $returnValue = 0;
    my $number = $_[0];
    for(my $i = 2; $i <= $number; $i++) {
        $returnValue = $returnValue + log($i);
    }
    return $returnValue;
}
