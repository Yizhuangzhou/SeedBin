#!usr/bin/perl -w
use strict;
die "perl $0 [CP005957_AMPI01.coords][output][output|remain]" unless @ARGV==3;
open(IN,$ARGV[0])||die;
open(OUT,">$ARGV[1]")||die;
open(OF,">$ARGV[2]")||die;
foreach(1 .. 4){
	<IN>;
}
my %flag;
my %len;
while(<IN>){
	chomp;
	my @a=split /\t/;
	$len{$a[-1]}=$a[8];
	my $min=$a[2];
	my $max=$a[3];
	if($a[3] <$a[2]){
		$min=$a[3];
		$max=$a[2];
	}
	foreach($min .. $max){
		$flag{$a[-1]}{$_}=1;
	}
}
close IN;
my %tag;
foreach my $scaf(keys %len){
	my $num=scalar keys %{$flag{$scaf}};
	if($num/$len{$scaf} >=0.9){
		$tag{$scaf}=1;
	}
}
open(IN,$ARGV[0])||die;
foreach(1 .. 4){
	<IN>;
}
while(<IN>){
	chomp;
	my @a=split /\t/;
	if($tag{$a[-1]}){
		print OUT "$_\n";
	}
	else{
		print OF "$_\n";
	}
}
close OUT;
close OF;
close IN;
