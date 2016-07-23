#!usr/bin/perl -w
use strict;
die "perl $0 [classified_allscaffold.xls][CP005957_AMPI01.coords.filter][table][stat]" unless @ARGV==4;
open(IF,$ARGV[0])||die;
open(IN,$ARGV[1])||die;
open(OUT,">$ARGV[2]")||die;
open(OF,">$ARGV[3]")||die;
my %hash;
my (%len,%depth);
my $binnum;
my %bin;
my $total=0;
my $unassign=0;
while(<IF>){
	chomp;
	my @a=split /\t/;
	$hash{$a[0]}=$a[1];
	$len{$a[0]}=$a[2];
	$depth{$a[0]}=$a[3];
	$total+=$len{$a[0]};
	if($a[1] ne "unassigned"){
		$bin{$a[1]}++;
	}
	else{
		$unassign+=$len{$a[0]};
	}
}
close IF;
$binnum=scalar keys %bin;
my %flag;
while(<IN>){
	chomp;
	my @a=split /\t/;
	$flag{$a[-1]}=1;
}
close IN;
my %sp;
my $sum=0;
my $individual_unassigned=0;
foreach my $k(keys %flag){
	next if(!$len{$k});
	print OUT "$k\t$hash{$k}\t$depth{$k}\t$len{$k}\n";
	$sum+=$len{$k};
	if($hash{$k} ne "unassigned"){
		$sp{$hash{$k}}+=$len{$k};
	}
	else{
		$individual_unassigned+=$len{$k};
	}
}
close OUT;
my @sp=sort {$sp{$b} <=> $sp{$a}} keys %sp;
my $sp=$sp[0];
my $sptotal=0;
foreach(keys %hash){
	if($hash{$_} eq $sp){
		$sptotal+=$len{$_};
	}
}
my $precision=sprintf("%.2f",$sp{$sp}*100/$sptotal);
my $sensitivity=sprintf("%.2f",$sp{$sp}*100/$sum);
my $assign=sprintf("%.2f",($total-$unassign)*100/$total);
my $iassign=sprintf("%.2f",($sum-$individual_unassigned)*100/$sum);
print OF "Precision(%)\tSencitivity(%)\tAssign(%)\tIndividual_Assigned\tUnassigned\tBinNum\n";
print OF "$precision\t$sensitivity\t$assign\t$iassign\t$unassign\t$binnum\n";
close OF;
