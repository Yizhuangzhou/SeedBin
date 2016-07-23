#!usr/bin/perl -w
use strict;
die "perl $0 [genome algin][output]" unless @ARGV==2;
open(IF,$ARGV[0])||die;
open(OUT,">$ARGV[1]")||die;
my (%total,%hash,%direct);
foreach(1 .. 4){
	<IF>;
}
while(<IF>){
	chomp;
	my @a=split /\t/;
	if(! $total{$a[-1]}{$a[-2]}){
		if($a[7]<$a[8]){
			$total{$a[-1]}{$a[-2]}=$a[7];
		}
		else{
			$total{$a[-1]}{$a[-2]}=$a[8];
		}
	}
	my $min=0;
	my $max=0;
	if($a[7]>$a[8]){
		if($a[2] > $a[3]){
			$min=$a[3];
			$max=$a[2];
		}
		else{
			$min=$a[2];
			$max=$a[3];
		}
	}
	else{
		if($a[0]>$a[1]){
			$min=$a[1];
			$max=$a[0];
		}
		else{
			$min=$a[0];
			$max=$a[1];
		}
	}
	foreach($min .. $max){
		$hash{$a[-1]}{$a[-2]}{$_}++;
	}
}
close IF;
my %flag;
foreach my $k1(keys %hash){
	my $maxp=0;
	my $id="";
	foreach my $k2 (keys %{$hash{$k1}}){
		my $n=scalar keys %{$hash{$k1}{$k2}};
		my $total=$total{$k1}{$k2};
		my $p=$n*100/$total;
		if($maxp ==0){
			$maxp=$p;
			$id=$k2;
		}
		elsif($p >$maxp){
			$id=$k2;
			$maxp=$p;
		}
	}
	$flag{$k1}{$id}=1;
}
open(IF,$ARGV[0])||die;
foreach(1 ..3){
	<IF>;
}
my $f=<IF>;
print OUT "$f";
while(<IF>){
	chomp;
	my @a=split /\t/;
	if($flag{$a[-1]}{$a[-2]}){
		print OUT "$_\n";
	}
}
close IF;
close OUT;

