#!usr/bin/perl -w
use strict;
die "perl $0 [F1RT_allscaffold_source.xls][classified_allscaffold.xls][minlen][output][output|total][log]" unless @ARGV==6;
open(IF,$ARGV[0])||die;
open(IN,$ARGV[1])||die;
my $minlen=$ARGV[2];
open(OUT,">$ARGV[3]")||die;
open(OF,">$ARGV[4]")||die;
open(LOG,">$ARGV[5]")||die;
my %source;
my %plasmid;
<IF>;
while(<IF>){
	chomp;
	my @a=split /\t/;
	$source{$a[0]}=$a[1];
}
close IF;
my %len;
my %hash;
my %flag;
my %depth;
my %unassigned=();
my $sum;
my $totalassigned=0;
my $totalunassigned=0;
while(<IN>){
	chomp;
	my @a=split /\t/;
	$len{$a[0]}=$a[2];
	if($a[1] ne "unassigned"){
		if($source{$a[0]}){
			$totalassigned+=$a[2];
			$hash{$a[1]}{$source{$a[0]}}+=$a[2];
			push @{$flag{$a[1]}},$a[0];
		}
	}
	else{
		if($source{$a[0]}){
			$unassigned{$source{$a[0]}}+=$a[2];
			$totalunassigned+=$a[2];
		}
	}
	if($source{$a[0]}){
		print LOG "$a[0]\t$a[1]\t$source{$a[0]}\t$a[2]\n";
	}
}
close IN;
close LOG;
my %sum1;
foreach my $k(keys %source){
	if($len{$k} && $len{$k} >=$minlen){
		$sum1{$source{$k}}+=$len{$k};
		$sum+=$len{$k};
	}
}
my %tag;

my $total_4precision=0;
foreach my $bin(sort {$a <=> $b} keys %hash){
	my @species=sort {$hash{$bin}{$b} <=> $hash{$bin}{$a}} keys %{$hash{$bin}};
	my $species=$species[0];
	$total_4precision+=$hash{$bin}{$species};
	$tag{$species}{$bin}=$hash{$bin}{$species};
}

#print OUT "Bin\tSpecies\tCompleteness\tPurity\tUnassigned\n";
my $num=0;
my $total_4sensitivity=0;
foreach my $sp(keys %tag){
	my @bin=sort {$tag{$sp}{$b} <=> $tag{$sp}{$a}} keys %{$tag{$sp}};
	my $bin=$bin[0];
	my $n=$hash{$bin}{$sp};
	$total_4sensitivity+=$n;
	my @scaffold=@{$flag{$bin}};
	my $bsum1=0;
	foreach(@scaffold){
		$bsum1+=$len{$_};
	}
	my $individual_precision=$n/$bsum1;
	my $individual_sensitivity=sprintf("%.2f",$n*100/$sum1{$sp});
	$num +=$n*$individual_precision;
#	$num +=$n;
#	my $tmpnum=int($n*$PI);
	$individual_precision=sprintf("%.2f",$individual_precision*100);
	if(!$unassigned{$sp}){
		$unassigned{$sp}=0;
	}
	my $individual_assigned=sprintf("%.2f",($sum1{$sp}-$unassigned{$sp})*100/$sum1{$sp});
	print OUT "$bin\t$sp\t$individual_precision\t$individual_sensitivity\t$individual_assigned\n";
}
close OUT;
my $performance=sprintf("%.2f",$num*100/$sum);
my $precision=sprintf("%.2f",$total_4precision*100/$totalassigned);
my $sensitivity=sprintf("%.2f",$total_4sensitivity*100/$sum);
my $assignedIndex=sprintf("%.2f",($sum-$totalunassigned)*100/$sum);
print OF "$performance\t$precision\t$sensitivity\t$assignedIndex\n";	
close OF;
