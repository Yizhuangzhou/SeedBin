#!usr/bin/perl -w
use strict;
die "perl $0 [*.depth list][output][morelen]" unless @ARGV==3;
open(LIST,$ARGV[0])||die;
open(OUT,">$ARGV[1]")||die;
my $morelen=$ARGV[2];
my %hash;
my %flag;
my @name=();
while(<LIST>){
	chomp;
	my @tmp=split /\t/;
	open(IF,$tmp[0])||die;
	my $num=$tmp[1];
	$/=">";
	<IF>;
	while(<IF>){
		chomp;
		if(/^(\S+)/){
			my $id=$1;
			my @a=split /\n/;
			shift @a;
			my @b=split /\s+/,join("",@a);
			next if(@b <$morelen);
			$flag{$id}++;
			if($flag{$id} ==1){
				push @name,$id;
			}
			foreach(1 .. $num){
				shift @b;
				pop @b;
			}
			my @c=();
			my $sum=0;
			foreach(@b){
				$sum+=$_;
			}
			my $str="";
			if($sum ==0){
				$str="0\t0";
			}
			else{
				my $n=scalar @b;
				my $mean=sprintf("%.2f",$sum/$n);
				$str="$mean";
				my $data;
				if($n <500){
					$str.="\t$mean";
				}
				else{
					foreach(my $i=0;$i<=$n-500;$i+=250){
						my $tmp=0;
						foreach(my $j=$i;$j<=$i+499;$j++){
							$tmp+=$b[$j];
						}
						$data=$tmp/500;
						$data=sprintf("%.2f",$data);
						$str.="\t$data";
					}
				}
			}
			if($hash{$id}){
				$hash{$id}.="\#$str";
			}
			else{
				$hash{$id}=$str;
			}
		}
	}
	close IF;
}
close LIST;
foreach my $scaf(@name){
	print OUT "$scaf\t$hash{$scaf}\n";
}
close OUT;

sub sd{
	my ($ref,$mean)=@_;
	my @val=@{$ref};
	my $sum=0;
	foreach(@val){
		$sum +=($_-$mean)**2;
	}
	$sum /=$#val;
	$sum =$sum **0.5;
	return $sum;
}
