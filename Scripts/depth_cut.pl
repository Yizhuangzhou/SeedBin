#!usr/bin/perl -w
use strict;
die "perl $0 [fa][cut.fa][depth][output]" unless @ARGV==4;
open(IF,$ARGV[0])||die;
open(IN,$ARGV[1])||die;
open(DP,$ARGV[2])||die;
open(OUT,">$ARGV[3]")||die;
my $id;
my %seq;
while(<IF>){
	chomp;
	if(/^>(\S+)/){
		$id=$1;
	}
	else{
		$seq{$id}.=$_;
	}
}
close IF;
my $pid;
my %flag;
my %sub;
my @name=();
my %pid;
while(<IN>){
	chomp;
	if(/^>(([^\.\s]+)\.\d+)/){
		$id=$1;
		$pid=$2;
		$pid{$id}=$pid;
		push @name,$id;
	}
	elsif(/^>([^\.\s]+)/){
		$id=$1;
		push @name,$id;
	}
	else{
		$sub{$id}.=$_;
	}
}
close IN;
my %index;
my %len;
foreach my $k(keys %pid){
	if(!$seq{$pid{$k}}){
		print "$k\n";
	}
	my $index=index($seq{$pid{$k}},$sub{$k});
	my $len=length $sub{$k};
	$index{$k}=$index;
	$len{$k}=$len;
}
my %depth;
$/=">";
while(<DP>){
	chomp;
	s/^/>/;
	if(/>(\S+)/){
		$id=$1;
		my @array=split /\n/;
		shift @array;
		my $str=join(" ",@array);
		my @depth=split /\s+/,$str;
		$depth{$id}=\@depth;
	}
}
close DP;
foreach my $id(@name){
	if($pid{$id}){
		$pid=$pid{$id};
		my $index=$index{$id};
		my $len=$len{$id};
		my $sum=0;
		my @depth=@{$depth{$pid}};
		for($index .. ($index+$len-1)){
			$sum+=$depth[$_];
		}
		my $mean=sprintf("%.2f",$sum/$len);
		print OUT "$id\t$mean\n";
	}
	else{
		my @depth=@{$depth{$id}};
		my $sum=0;
		foreach(@depth){
			$sum+=$_;
		}
		my $mean=sprintf("%.2f",$sum/@depth);
		print OUT "$id\t$mean\n";
	}
}
close OUT;
