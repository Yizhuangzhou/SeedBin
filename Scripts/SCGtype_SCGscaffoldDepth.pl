#!usr/bin/perl -w
use strict;
use FindBin qw($Bin);
use File::Basename;
use POSIX qw(strftime);
die "perl $0 [fa][depth][outdir][minlen]" unless @ARGV==4;
####gene prediction
my @suffix=(".fa",".fna",".fasta");
my $base=basename($ARGV[0],@suffix);
my $outdir=$ARGV[2];
mkdir $outdir unless (-e $outdir);
my $minlen=$ARGV[3];
my $fafile="$outdir/$base\.fa";
open(FA,$ARGV[0])||die;
open(OF,">$fafile")||die;
my $name="";
my %length=();
my $current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
print STDERR "Start .... at $current_time\n";
while(<FA>){
	chomp;
	if(/>(\S+)/){	
		$name=$1;
		$length{$name}=0;
	}
	else{
		s/[^ATGC]//gi;
		$length{$name}+=length;
	}
}
close FA;
open(FA,$ARGV[0])||die;
while(<FA>){
	chomp;
	if(/>(\S+)/){
		$name=$1;
		if($length{$name} >=$minlen){
			print OF "$_\n";
		}
	}
	elsif($length{$name} >=$minlen){
		print OF "$_\n";
	}
}
close FA;
close OUT;
$current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
print STDERR "Finish extracting scaffolds with length >=$minlen at $current_time\n";
my $cmd="perl /processing_data/zhouyzh/software/FragGeneScan1.19/run_FragGeneScan.pl -genome=$fafile -out $outdir/$base -complete=0 -train=complete -thread 1 1 >$outdir/$base.frag.out 2>$outdir/$base.frag.err";
system($cmd);
$current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
print STDERR "Finish gene prediction at $current_time\n";
my $proteinfile="$outdir/$base\.faa";
$cmd="/data/software/pool/hmmer-3.1b1-linux-intel-x86_64/binaries/hmmsearch --cut_tc --tblout $outdir/$base\_Pfam.tblout --domtblout $outdir/$base\_Pfam.domtblout -o $outdir/$base\_Pfam.out $Bin/singleCopygene_fromPfam.HMM $proteinfile";
system($cmd);
$cmd="/data/software/pool/hmmer-3.1b1-linux-intel-x86_64/binaries/hmmsearch --cut_tc --tblout $outdir/$base\_TIGR.tblout --domtblout $outdir/$base\_TIGR.domtblout -o $outdir/$base\_TIGR.out $Bin/singleCopygene_fromTIGR.HMM $proteinfile";
system($cmd);
$current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
print STDERR "Finish performing hmmsearch at $current_time\n";
open(TMP,"$Bin/singleCopygene_samegene.xls")||die;
my %flag;
while(<TMP>){
	chomp;
	my @a=split /\t/;
	foreach(@a){
		$flag{$_}="$a[0]|$a[1]";
	}
}
close TMP;
my %minEvalue=();
my %id=();
&SCG2gene("$outdir/$base\_Pfam.tblout");
&SCG2gene("$outdir/$base\_TIGR.tblout");
my %SCG2scaffold=();
open(IN,"$outdir/$base\.faa")||die;
while(<IN>){
	chomp;
	if(/>((\S+)\_\d+\_\d+\_[-+])$/){
		my $scaffold=$2;
		my $gene=$1;
		if($id{$gene}){
			push @{$SCG2scaffold{$id{$gene}}},$scaffold;
		}
	}
}
close IN;
my %depth;
my %depthstr;
open(DP,$ARGV[1])||die;
while(<DP>){
	chomp;
	my @a=split /\t/,$_,2;
	my $id=$a[0];
	my @array=split /\#/,$a[1];
	my $index=0;
	foreach(@array){
		$index++;
		my @data=split /\t/;
		if($index == 1){
			$depth{$id}=$data[0];
			$depthstr{$id}=$data[0];
		}
		else{
			$depthstr{$id}.=",$data[0]";
		}
	}
}
close DP;
my @SCGid=sort {$a cmp $b} keys %SCG2scaffold;
open(OUT,">$outdir/SCGtype_SCGscaffoldDepth.xls")||die;
foreach my $scg(@SCGid){
	print OUT "$scg";
	my @scaffold=@{$SCG2scaffold{$scg}};
	my @tmpscaffold=();
	foreach(@scaffold){
		if($depth{$_}){
			push @tmpscaffold,$_;
		}
	}
	@scaffold=@tmpscaffold;
	@scaffold=sort {$depth{$b} <=> $depth{$a}} @scaffold;
	foreach(@scaffold){
		print OUT "\t$_($depthstr{$_},$length{$_})";
	}
	print OUT "\n";
}
close OUT;
$current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
print STDERR "Finish all at $current_time\n";

sub SCG2gene{
	my ($file)=@_;
	open(IF,$file)||die;
	while(<IF>){
		chomp;
		next if(/^\#/);
		my @a=split /\s+/;
		if(!$minEvalue{$a[0]}){
			$minEvalue{$a[0]}=$a[4];
			if($flag{$a[3]}){
				$id{$a[0]}=$flag{$a[3]};
			}
			else{
				$id{$a[0]}=$a[3];
			}
		}
		elsif($a[4] <$minEvalue{$a[0]}){
			$minEvalue{$a[0]}=$a[4];
			if($flag{$a[3]}){
				$id{$a[0]}=$flag{$a[3]};
			}
			else{
				$id{$a[0]}=$a[3];
			}
		}
	}
	close IF;
}
