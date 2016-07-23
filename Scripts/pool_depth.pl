#!usr/bin/perl -w
use strict;
die "perl $0 [HP-_cut.depth][HP+_cut.depth][output]" unless @ARGV==3;
open(IF,$ARGV[0])||die;
open(IN,$ARGV[1])||die;
open(OUT,">$ARGV[2]")||die;
while(<IF>){
	chomp;
	my @a=split /\t/;
	my $line=<IN>;
	chomp($line);
	my @b=split /\t/,$line;
	print OUT "$a[0]\t$a[1]\t$b[1]\n";
}
close IF;
close IN;
close OUT;
