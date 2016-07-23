#!usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use FindBin qw($Bin);
use lib "$Bin";
use Statistics::Distributions qw(uprob);
use Cwd;
my $dir=getcwd();
use Cwd 'abs_path';
use File::Path qw(make_path);
use threads;
use threads::shared;
use POSIX qw(strftime);

my ($fafile,$tablefile,$depthfile,$seedfile,$nb_process,$fold,$minlen,$outdir,$LargeLen,$SmallLen,$BinLen,$scgNum,$scgNum4bin,$help);
GetOptions(
	"fafile:s" => \$fafile,
	"tablefile:s" => \$tablefile,
	"depthfile:s" => \$depthfile,
	"seedfile:s" => \$seedfile,
	"nb_process:i" => \$nb_process,
	"minlen:i" => \$minlen,
	"fold:f" => \$fold,
	"outdir:s" =>\$outdir,
	"LargeLen:i" =>\$LargeLen,
	"SmallLen:i" =>\$SmallLen,
	"BinLen:i" =>\$BinLen,
	"scgNum:i" =>\$scgNum,
	"scgNum4bin:i" =>\$scgNum4bin,
	"help" => \$help,
);

$nb_process=20;
$fold ||= 2.58;
$minlen ||=0;
$LargeLen ||= 200000;
$SmallLen ||= 100000;
$BinLen ||= 500000;
$scgNum	||= 7;
$scgNum4bin ||= 10;
$outdir ||=$dir;

unless ($fafile  && $depthfile && $tablefile){
	&usage;
	exit;
}
if($help){
	&usage;
	exit;
}
$outdir=abs_path($outdir);
make_path($outdir) unless (-e $outdir);

my $base=basename($0,".pl");
open(LOG,">$outdir/$base\.log")||die;
my $current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
print LOG "$current_time: Start to conduct $0......\n";
print STDERR "$current_time: Start to conduct $0 ......\n";


my (%mean,%sd,%depth,%rawsd,%min,%max);
my %foldval;
my %scaf2cl;
my @NonSCGscaffold_more10kb=();
my @NonSCGscaffold_less10kb=();
my %prob=();
my %seed2cl;
my ($len,$refxx,$SSxx,$len1,$len2);
my (%refxx,%SSxx,%len,%seedmin,%seedmax,%seedmean);
my @totalseed=();
my %assign=();
my %used;
my %class=();
my $fseed="";
my %clseq;
my @scaffold_more10kb_ScaffoldComb=();
my $IndexStart=0;
my @running=();
my @Threads=();
my %joined=();
my @possible_seed=();
my %psource=();
my %seq;
my %length;
my %Flagscg2cl=();
my %exclude=();
my %SCGassign=();
my %classified=();
my @tmptotalseed=();
my %seedCmin=();
my %seedCmax=();
my %scg2cl=();
my %seedscg=();
my %seedscgNum=();
my %coreassign=();
my %idassign=();
my $INum;
my %freq=();
my %result=();
my @remain_more10kb=();
my @remain_less10kb=();
my %scaffold;

####################################################################### Intializing for Zvalue ##############################################
my @mono=("A","T","G","C");
our @oligo_kmer=();
foreach(@mono){
	my $word=$_;
	foreach (@mono){
		my $di=$word."$_";
		foreach(@mono){
			my $tri=$di."$_";
			foreach(@mono){
				my $tetra=$tri."$_";
				push @oligo_kmer,$tetra;
			}
		}
	}
}
our (@oligo_k_1mer,@oligo_k_2mer)=((),());
foreach(@mono){
	my $word=$_;
	foreach(@mono){
		my $di=$word."$_";
		foreach(@mono){
			my $tri=$di."$_";
			push @oligo_k_1mer,$tri;
		}
	}
}
foreach(@mono){
	my $word=$_;
	foreach(@mono){
		my $di=$word."$_";
		push @oligo_k_2mer,$di;
	}
}

############################################################### Read cutoff file ####################################################
open(CF,"$Bin/cutoff.xls")||die;
my (%smean,%ssd,%dmean,%dsd,%cutoff,%diff);
while(<CF>){
	chomp;
	my @a=split /\t/;
	my @b=split /\s*vs\s*/,$a[0];
	$smean{"$b[0]\t$b[1]"}=$a[1];
	$smean{"$b[1]\t$b[0]"}=$a[1];
	$ssd{"$b[0]\t$b[1]"}=$a[2];
	$ssd{"$b[1]\t$b[0]"}=$a[2];
	$dmean{"$b[0]\t$b[1]"}=$a[3];
	$dmean{"$b[1]\t$b[0]"}=$a[3];
	$dsd{"$b[0]\t$b[1]"}=$a[4];
	$dsd{"$b[1]\t$b[0]"}=$a[4];
	my @tmp=($a[5],$a[6],$a[7]);
	@tmp=sort{$b<=>$a} @tmp;
	$cutoff{"$b[0]\t$b[1]"}=$tmp[0];
	$cutoff{"$b[1]\t$b[0]"}=$tmp[0];
	$diff{"$b[0]\t$b[1]"}=$a[5];
	$diff{"$b[1]\t$b[0]"}=$a[5];
}
close CF;

############################################################### Read fasta file ####################################################
open(FA,$fafile)||die;
my $id;
while(<FA>){
	chomp;
	if(/>(\S+)/){
		$id=$1;
		$seq{$id}="";
		$length{$id}=0;
	}
	else{
		s/[^ATGC]//gi;
		my $temp=uc($_);
		$seq{$id}.=$temp;
		$length{$id}+=length $temp;
	}
}
close FA;

############################################################### Read Table file ####################################################
my %scg;
my @scgid=();
open(TMP,$tablefile)||die;
while(<TMP>){
	chomp;
	my @a=split /\t/;
	my $scg=shift @a;
	push @scgid,$scg;
	foreach(@a){
		if(/^([^\(]+)/){
			$id=$1;
			if($length{$id} >=$minlen){
				$scg{$id}{$scg}=1;
			}
		}
	}
}
close TMP;

my %duplicate=();
foreach my $scaf(keys %scg){
	if(scalar keys %{$scg{$scaf}} >=2){
		foreach my $scaf2(keys %scg){
			if($scaf eq $scaf2 || scalar keys %{$scg{$scaf2}} ==1){
				next;
			}
			else{
				my $n=0;
				foreach my $scg(keys %{$scg{$scaf}}){
					if($scg{$scaf2}{$scg}){
						$n++;
					}
				}
				if($n >=2){
					if($duplicate{$scaf}){
						$duplicate{$scaf}.=";$scaf2";
					}
					else{
						$duplicate{$scaf}="$scaf2";
					}
				}
			}
		}
	}
}

###l############################################################ Read depth file ####################################################
open(DP,$depthfile)||die;
while(<DP>){
	chomp;
	my @a=split /\t/,$_,2;
	my $id=$a[0];
	next if($length{$id} <$minlen);
	my @array=split /\#/,$a[1];
	my $index=0;
	foreach(@array){
		$index++;
		my @data=split /\t/;
		my $mean=shift @data;
		$mean{$id}{$index}=$mean;
		$depth{$id}{$index}=\@data;
		if(scalar @data <10){
			$sd{$id}{$index}="NA";
		}
		else{
			my ($tmin,$tmax)=&depth_process(\@data);
			my @tmpdata=();
			foreach(@data){
				if($_ >=$tmin && $_ <=$tmax){
					push @tmpdata,$_;
				}
			}
			@tmpdata=sort {$a<=>$b} @tmpdata;
			my $sd=&sd(\@tmpdata,$mean);
			if($mean >1){
				$sd /=$mean;
				if($sd <0.38){
					push @{$rawsd{$index}},$sd;
				}
			}
			else{
				$sd=0;
			}
			$sd{$id}{$index}=$sd;
		}
	}
}
close DP;

my (%sd1,%sd2);
my @CovIndex=sort {$a<=>$b} keys %rawsd;
my $IndexNum=scalar @CovIndex;
foreach my $k(@CovIndex){
	my @data=@{$rawsd{$k}};
	@data=sort {$a<=>$b} @data;
	if(@data % 2 ==1){
		my $index=$#data/2;
		$sd1{$k}=sprintf("%.2f",$data[$index]);
	}
	else{
		my $index=@data/2;
		$sd1{$k}=sprintf("%.2f",($data[$index-1]+$data[$index])/2);
	}
	my $index= int(@data*3/4)-1;
	$sd2{$k}=sprintf("%.2f",$data[$index]);
	$current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
	print STDERR "$current_time: Mean SD ratio of Sample$k\t$sd1{$k}\t$sd2{$k}\n";
	print LOG "$current_time: Mean SD ratio of Sample$k\t$sd1{$k}\t$sd2{$k}\n";
}
%rawsd=();
foreach my $scaf (keys %length) {
	next if($length{$scaf} <$minlen);
	foreach my $index (@CovIndex) {
		my $sd1=$sd1{$index};
		my $sd2=$sd2{$index};
		my $sd;
		if($length{$scaf} >=10000){
			if($sd{$scaf}{$index} eq "NA"){
				$sd=$sd1;
			}
			elsif($sd{$scaf}{$index} <$sd1){
				$sd=$sd1;
			}
			elsif($sd{$scaf}{$index} >=$sd1 && $sd{$scaf}{$index} <=$sd2){
				$sd=$sd{$scaf}{$index};
			}
			elsif($sd{$scaf}{$index} >$sd2){
				$sd=$sd2;
			}
			$seedmin{$scaf}{$index}=int($mean{$scaf}{$index}*(1-$sd*$fold));
			$seedmax{$scaf}{$index}=int($mean{$scaf}{$index}*(1+$sd*$fold))+1;
		}
		else{
			$min{$scaf}{$index}=int($mean{$scaf}{$index}*(1-$sd1));
			$max{$scaf}{$index}=int($mean{$scaf}{$index}*(1+$sd1))+1;
		}
	}
}
%sd=();

###l############################################################ Read seed file ####################################################
if($seedfile){
	open(SD,$seedfile)||die;
	while(<SD>){
		chomp;
		my @scaffold=split /;/;
		my $clstr=$_;
		if($fseed eq ""){
			$fseed=$clstr;
		}
		my %sumdepth=();
		my $sumlen=0;
		my $seq="";
		my %data=();
		foreach my $scaf (@scaffold) {
			$seq.=$seq{$scaf};
			$sumlen+=$length{$scaf};
			$used{$scaf}=1;
			foreach my $index(@CovIndex) {
				$sumdepth{$index}+=$mean{$scaf}{$index}*$length{$scaf};
				push @{$data{$index}},@{$depth{$scaf}{$index}};
			}
		}
		$clseq{$clstr}=$seq;
		push @totalseed,$clstr;
		$class{$clstr}=$clstr;
		foreach my $index(@CovIndex) {
			$seedmean{$clstr}{$index}=int($sumdepth{$index}/$sumlen);
			my @data=@{$data{$index}};
			my $sd;
			if(scalar @data >=3){
				my ($tmin,$tmax)=&depth_process(\@data);
				my @tmpdata=();
				foreach(@data){
					if($_ >=$tmin && $_ <=$tmax){
						push @tmpdata,$_;
					}
				}
				$sd=&sd(\@tmpdata,$seedmean{$clstr}{$index});
				if($seedmean{$clstr}{$index} >=5){
					$sd /=$seedmean{$clstr}{$index};
					if($sd <$sd1{$index}){
						$sd=$sd1{$index};
					}
					elsif($sd >$sd2{$index}){
						$sd=$sd2{$index};
					}
				}
				else{
					$sd=$sd1{$index};
				}
			}
			else{
				$sd=$sd1{$index};
			}
			$seedmin{$clstr}{$index}=int($seedmean{$clstr}{$index}*(1-$sd*$fold));
			$seedmax{$clstr}{$index}=int($seedmean{$clstr}{$index}*(1+$sd*$fold))+1;
		}
		my $len=$sumlen;
		if(!$refxx{$clstr}){
			$len=int($len/10000)*10;
			if($len >210){
				$len=210;
			}
			my $zvalueref=&zvalue($seq);
			my $mean=&mean($zvalueref);
			my @d=();
			foreach(@{$zvalueref}){
				my $tmp=$_-$mean;
				push @d,$tmp;
			}
			$refxx=\@d;
			$SSxx=&SS($refxx,$refxx,__LINE__);
			$refxx{$clstr}=$refxx;
			$SSxx{$clstr}=$SSxx;
			$len{$clstr}=$len;
		}
	}
	close SD;
}

############################################################# split scaffolds into large or small ones ############################################
my $m=0;
my @scaffold_less10kb=();
my @scaffold_more10kb=();
my @totalSCGscaffold_more10kb=();
my @totalSCGscaffold_less10kb=();
foreach my $scaf (keys %length) {
	next if($length{$scaf} <$minlen);
	if(!$used{$scaf}){
		if($length{$scaf} <10000){
			push @scaffold_less10kb,$scaf;
			if($scg{$scaf}){
				push @totalSCGscaffold_less10kb,$scaf;
			}
			else{
				push @NonSCGscaffold_less10kb,$scaf;
			}
		}
		else{
			push @scaffold_more10kb,$scaf;
			if($scg{$scaf}){
				push @totalSCGscaffold_more10kb,$scaf;
			}
			else{
				push @NonSCGscaffold_more10kb,$scaf;
			}
		}
	}
}
$m=scalar @scaffold_more10kb;
$current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
print LOG "$current_time: The number of all scaffolds with length >=10kb: $m\n";
print STDERR "$current_time: The number of all scaffolds with length >=10kb: $m\n";
$m=scalar @scaffold_less10kb;
print LOG "$current_time: The number of all scaffolds with length <10kb: $m\n";
print STDERR "$current_time: The number of all scaffolds with length <10kb: $m\n";
$m=scalar @totalSCGscaffold_more10kb;
$current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
print LOG "$current_time: The number of SCG scaffolds with length >=10kb: $m\n";
print STDERR "$current_time: The number of SCG scaffolds with length >=10kb: $m\n";
$m=scalar @totalSCGscaffold_less10kb;
print LOG "$current_time: The number of SCG scaffolds with length <10kb: $m\n";
print STDERR "$current_time: The number of SCG scaffolds with length <10kb: $m\n";
$m=scalar @NonSCGscaffold_more10kb;
$current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
print LOG "$current_time: The number of Non-SCG scaffolds with length >=10kb: $m\n";
print STDERR "$current_time: The number of Non-SCG scaffolds with length >=10kb: $m\n";
$m=scalar @NonSCGscaffold_less10kb;
print LOG "$current_time: The number of Non-SCG scaffolds with length <10kb: $m\n";
print STDERR "$current_time: The number of Non-SCG scaffolds with length <10kb: $m\n";

############################################################### Preprocessing large and small scaffolds ######################################################
foreach my $scaf (@scaffold_more10kb) {
	if(!$refxx{$scaf}){
		my $seq=$seq{$scaf};
		&Prepare4Zvalue($seq,$scaf);
	}
}
foreach my $scaf (@scaffold_less10kb) {
	&kmerFreq($scaf);
}
my @compositionSeed=@totalseed;

############################################################### obtain seeds of large scaffolds ######################################################
print STDERR "$current_time: Start to obtain seeds of large scaffolds\n";
print LOG "$current_time: Start to obtain seeds of large scaffolds\n";
@scaffold_more10kb=sort{$mean{$b}{1} <=>$mean{$a}{1}} @scaffold_more10kb;
my $diffOrcutoff=1;
# 1 for diff and 0 for cutoff
my ($ref)=&largeScaffold_Seed(\@scaffold_more10kb,\@totalseed,$diffOrcutoff);
@totalseed=@{$ref};

open(OF,">$outdir/largeScaffoldSeed_CovInterval.xls")||die;
foreach my $seed (@totalseed) {
	print OF "$seed";
	foreach my $index (@CovIndex){
		print OF "\t$seedmin{$seed}{$index}\t$seedmax{$seed}{$index}";
	}
	print OF "\n";
}
close OF;
my $n=scalar @totalseed;
$current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
print LOG "$current_time: Number of seeds from large scaffolds: $n\n";
print LOG "$current_time: Produce $outdir/largeScaffoldSeed_CovInterval.xls\n";
print STDERR "$current_time: Number of seeds from large scaffolds: $n\n";
print STDERR "$current_time: Produce $outdir/largeScaffoldSeed_CovInterval.xls\n";

############################################################### obtain seeds  by only composition ###########################################################
my ($copref)=&CompositionSeed(\@scaffold_more10kb,\@compositionSeed);
@compositionSeed=@{$copref};
$n=scalar @compositionSeed;
$current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
print LOG "$current_time: Number of seeds from large scaffolds with only composition: $n\n";
print STDERR "$current_time: Number of seeds from large scaffolds with only composition: $n\n";
my %largeScaffold_compositionSeed=();
%clseq=();
%foldval=();
%assign=();
foreach my $cl (@compositionSeed) {
	my @scaffold=split /;/,$cl;
	foreach(@scaffold){
		$clseq{$cl}.=$seq{$_};
	}
}
&Classify_SmallScaffold_withoutCov(\@compositionSeed,\@scaffold_less10kb);
my %num=();
foreach my $cl (@compositionSeed) {
	$largeScaffold_compositionSeed{$cl}=1;
	if($assign{$cl}){
		my @scaffold=split /\;/,$assign{$cl};
		$num{$cl}=scalar @scaffold;
	}
	else{
		$num{$cl}=0;
	}
}
@compositionSeed=sort {$num{$b} <=> $num{$a}} @compositionSeed;
my %maxfoldval=();
foreach my $cl(@compositionSeed){
	next if(!$assign{$cl});
	my @scaffold=split /;/,$assign{$cl};
	@scaffold=sort {$foldval{$cl}{$b} <=> $foldval{$cl}{$a}} @scaffold;
	my @scafcomb=();
	for(my $i=0;$i<$#scaffold;$i++) {
		my $scaf1=$scaffold[$i];
		my $sumlen=$length{$scaf1};
		my $clstr=$scaf1;
		my $seq=$seq{$scaf1};
		for(my $j=$i+1;$j<=$#scaffold;$j++) {
			my $scaf2=$scaffold[$j];
			my $Num=0;
			my $cha=0;
			foreach my $index (@CovIndex) {
				if($mean{$scaf1}{$index} >=$mean{$scaf2}{$index}){
					if($mean{$scaf1}{$index} <5){
						$cha++;
					}
					elsif($mean{$scaf2}{$index} >=$min{$scaf1}{$index} && $mean{$scaf2}{$index} <=$max{$scaf1}{$index}){
						$Num++;
					}
				}
				else{
					if($mean{$scaf2}{$index} <5){
						$cha++;
					}
					elsif($mean{$scaf1}{$index} >=$min{$scaf2}{$index} && $mean{$scaf1}{$index} <=$max{$scaf2}{$index}){
						$Num++;
					}
				}
			}
			if($Num+$cha == $IndexNum){
				my $dist=int(abs($foldval{$cl}{$scaf1}-$foldval{$cl}{$scaf2}));
				if($dist <=7){
					$sumlen+=$length{$scaf2};
					$seq.=$seq{$scaf2};
					$clstr.=";$scaf2";
					if($sumlen >=10000){
						last;
					}
				}
				else{
					last;
				}
			}
			else{
				#print STDERR "$scaf1\t$scaf2\n";
			}
		}
		if($sumlen >=10000){
			&Prepare4Zvalue($seq,$clstr);
			push @scafcomb,$clstr;
		}
	}
	foreach my $scaf (@scaffold) {
		$maxfoldval{$scaf}{$cl}=$foldval{$cl}{$scaf};
		delete $foldval{$cl}{$scaf};
	}
	delete $foldval{$cl};
	my ($refyy,$SSyy,$SSxy);
	my %seed=();
	foreach my $clstr(@scafcomb) {
		$refxx=$refxx{$clstr};
		$len=$len{$clstr};
		$SSxx=$SSxx{$clstr};
		my $flag=1;
		foreach my $seed(@compositionSeed){
			my $refyy=$refxx{$seed};
			my $len2=$len{$seed};
			my $SSyy=$SSxx{$seed};
			my $SSxy=&SS($refxx,$refyy,__LINE__);
			my $cor=&correlation($SSxx,$SSyy,$SSxy);
			if($cor >=$diff{"$len\t$len2"}){
				$flag=0;
				last;
			}
		}
		if($flag ==1){
			push @compositionSeed,$clstr;
			$seed{$clstr}=1;
		}
	}
	foreach my $cl (@scafcomb) {
		if(!$seed{$cl}){
			delete $refxx{$cl};
			delete $SSxx{$cl};
			delete $len{$cl};
		}
	}
}
$n=scalar @compositionSeed;
my %seed=();
foreach my $seed (@totalseed) {
	$seed{$seed}=1;
}
foreach my $seed (@compositionSeed) {
	if(!$seed{$seed} && $seed=~/;/){
		delete $refxx{$seed};
		delete $SSxx{$seed};
		delete $len{$seed};
	}
}
$current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
print LOG "$current_time: Number of seeds from all scaffolds with only composition: $n\n";
print STDERR "$current_time: Number of seeds from all scaffolds with only composition: $n\n";

############################################################### obtain scaffold combinations ###########################################################
%clseq=();
%foldval=();
%assign=();
foreach my $cl (@compositionSeed) {
	my @scaffold=split /;/,$cl;
	foreach(@scaffold){
		$clseq{$cl}.=$seq{$_};
	}
}
&Classify_SmallScaffold_withoutCov(\@compositionSeed,\@scaffold_less10kb);
%maxfoldval=();
%largeScaffold_compositionSeed=();
my @pesudoscaf=();
foreach my $cl(@compositionSeed){
	next if(!$assign{$cl});
	my @scaffold=split /;/,$assign{$cl};
	@scaffold=sort {$foldval{$cl}{$b} <=> $foldval{$cl}{$a}} @scaffold;
	my @scafcomb=();
	for(my $i=0;$i<$#scaffold;$i++) {
		my $scaf1=$scaffold[$i];
		my $sumlen=$length{$scaf1};
		my $clstr=$scaf1;
		my %sumdepth=();
		my %data=();
		foreach my $index(@CovIndex){
			$sumdepth{$index}=$mean{$scaf1}{$index}*$length{$scaf1};
			push @{$data{$index}},@{$depth{$scaf1}{$index}};
		}
		for(my $j=$i+1;$j<=$#scaffold;$j++) {
			my $scaf2=$scaffold[$j];
			my $Num=0;
			my $cha=0;
			foreach my $index (@CovIndex) {
				if($mean{$scaf1}{$index} >=$mean{$scaf2}{$index}){
					if($mean{$scaf1}{$index} <5){
						$cha++;
					}
					elsif($mean{$scaf2}{$index} >=$min{$scaf1}{$index} && $mean{$scaf2}{$index} <=$max{$scaf1}{$index}){
						$Num++;
					}
				}
				else{
					if($mean{$scaf2}{$index} <5){
						$cha++;
					}
					elsif($mean{$scaf1}{$index} >=$min{$scaf2}{$index} && $mean{$scaf1}{$index} <=$max{$scaf2}{$index}){
						$Num++;
					}
				}
			}
			if($Num+$cha == $IndexNum){
				my $dist=int(abs($foldval{$cl}{$scaf1}-$foldval{$cl}{$scaf2}));
				if($dist <=7){
					$sumlen+=$length{$scaf2};
					$clstr.=";$scaf2";
					foreach my $index (@CovIndex) {
						$sumdepth{$index}+=$mean{$scaf2}{$index}*$length{$scaf2};
						push @{$data{$index}},@{$depth{$scaf2}{$index}};
					}
					if($sumlen >=10000){
						last;
					}
				}
				else{
					last;
				}
			}
			else{
				#print STDERR "$scaf1\t$scaf2\n";
			}
		}
		if($sumlen >=10000){
			foreach my $index (@CovIndex) {
				$seedmean{$clstr}{$index}=int($sumdepth{$index}/$sumlen);
				if($seedmean{$clstr}{$index} >=5){
					my @data=@{$data{$index}};
					my $sd;
					if(scalar @data >=3){
						my ($tmin,$tmax)=&depth_process(\@data);
						my @tmpdata=();
						foreach(@data){
							if($_ >=$tmin && $_<=$tmax){
								push @tmpdata,$_;
							}
						}
						$sd=&sd(\@tmpdata,$seedmean{$clstr}{$index});
						$sd /=$seedmean{$clstr}{$index};
						if($sd <$sd1{$index}){
							$sd=$sd1{$index};
						}
						elsif($sd >$sd2{$index}){
							$sd=$sd2{$index};
						}
					}
					else{
						$sd=$sd1{$index};
					}
					$seedmin{$clstr}{$index}=int($seedmean{$clstr}{$index}*(1-$sd*$fold));
					$seedmax{$clstr}{$index}=int($seedmean{$clstr}{$index}*(1+$sd*$fold))+1;
				}
				else{
					$seedmin{$clstr}{$index}=0;
					$seedmax{$clstr}{$index}=0;
				}
			}
			push @scafcomb,$clstr;
		}
	}
	push @pesudoscaf,@scafcomb;
}
%foldval=();
@compositionSeed=();

############################################################### obtain seeds of small scaffolds ###########################################################
@pesudoscaf=sort {$seedmean{$b}{1} <=> $seedmean{$a}{1}} @pesudoscaf;
my ($refyy,$SSyy,$SSxy);
foreach my $clstr(@pesudoscaf) {
	my $seq="";
	my @scaffold=split /\;/,$clstr;
	foreach(@scaffold){
		$seq.=$seq{$_};
	}
	&Prepare4Zvalue($seq,$clstr);
	$refxx=$refxx{$clstr};
	$len=$len{$clstr};
	$SSxx=$SSxx{$clstr};
	my $flag=1;
	foreach my $seed(@totalseed){
		my $Num=0;
		my $cha=0;
		foreach my $index (@CovIndex) {
			if($seedmean{$seed}{$index} <5){
				$cha++;
			}
			elsif($seedmean{$clstr}{$index} >=$seedmin{$seed}{$index} && $seedmean{$clstr}{$index} <=$seedmax{$seed}{$index}){
				$Num++;
			}
		}
		if($Num+$cha == $IndexNum){
			my $refyy=$refxx{$seed};
			my $len2=$len{$seed};
			my $SSyy=$SSxx{$seed};
			my $SSxy=&SS($refxx,$refyy,__LINE__);
			my $cor=&correlation($SSxx,$SSyy,$SSxy);
			if($cor >=$diff{"$len\t$len2"}){
				$flag=0;
				last;
			}
		}
	}
	if($flag ==1){
		push @totalseed,$clstr;
	}
	else{
		delete $refxx{$clstr};
		delete $SSxx{$clstr};
		delete $len{$clstr};
		delete $seedmean{$clstr};
		delete $seedmin{$clstr};
		delete $seedmax{$clstr};
	}
}
$n=scalar @totalseed;
$current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
print LOG "$current_time: The number of all seeds is $n\n";
print STDERR "$current_time: The number of all seeds is $n\n";

############################################################################ output all seeds and their coverage intervals ####################################
@totalseed=sort {$seedmean{$b}{1} <=> $seedmean{$a}{1}} @totalseed;
&CovTreat;
open(OUT,">$outdir/allscaffoldSeed_CovInterval.xls")||die;
foreach my $cl(@totalseed){
	if($seed2cl{$cl}){
		print OUT "$cl\t$seed2cl{$cl}";
	}
	else{
		print OUT "$cl\t$cl";
	}
	foreach my $index (@CovIndex) {
		print OUT "\t$seedmin{$cl}{$index}\t$seedmax{$cl}{$index}";
	}
	print OUT "\n";
}
close OUT;
$current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
print LOG "$current_time: Produce $outdir/allscaffoldSeed_CovInterval.xls\n";
print STDERR "$current_time: Produce $outdir/allscaffoldSeed_CovInterval.xls\n";

############################################################### classify large scaffolds #########################################################
%assign=();
foreach my $cl (keys %class) {
	if($assign{$cl}){
		$assign{$cl}.=";$class{$cl}";
	}
	else{
		$assign{$cl}=$class{$cl};
	}
}
$current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
print STDERR "$current_time: Start to classify large scaffolds\n";
print LOG "$current_time: Start to classify large scaffolds\n";
&Classify_LargeScaffold(\@totalseed,\@scaffold_more10kb);
my $TAG=1;
%clseq=();
%foldval=();
foreach my $cl (@totalseed) {
	my @scaffold=split /;/,$cl;
	foreach(@scaffold){
		$clseq{$cl}.=$seq{$_};
	}
}
$current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
print STDERR "$current_time: Start to classify small scaffolds\n";
print LOG "$current_time: Start to classify small scaffolds\n";
&Classify_SmallScaffold(\@totalseed,\@scaffold_less10kb,$TAG);

############################################################### output raw classified results ######################################################
my %scaf2sp=();
my @clid=sort {$seedmean{$b}{1} <=> $seedmean{$a}{1}} @totalseed;
open(OF,">$outdir/allscaffold_RawClassified.xls")||die;
foreach my $cl(@clid){
	if($assign{$cl}){
		print OF "$cl\t$assign{$cl}\n";
		my @clscaf=split /;/,$assign{$cl};
		foreach my $scaf(@clscaf){
			$scaf2sp{$scaf}=$cl;
		}
	}
	else{
		print OF "$cl\tNA\n";
	}
}
close OF;
$current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
print LOG "$current_time: Produce $outdir/allscaffold_RawClassified.xls\n";
print STDERR "$current_time: Produce $outdir/allscaffold_RawClassified.xls\n";

open(OUT,">$outdir/SCGscaffold_RawClassified.xls")||die;
print OUT "SCG";
foreach(@clid){
	print OUT "\t$_";
}
print OUT "\n";
open(TMP,$tablefile)||die;
while(<TMP>){
	chomp;
	my @a=split /\t/;
	print OUT "$a[0]";
	shift @a;
	my %temp_hash=();
	foreach(@a){
		if(/^([^\(]+)\(/){
			my $scaf=$1;
			my $str="";
			foreach my $index(@CovIndex){
				$str.=",$mean{$scaf}{$index}";
			}
			$str=~s/^\,//;
			if($scaf2sp{$scaf}){
				if($temp_hash{$scaf2sp{$scaf}}){
					$temp_hash{$scaf2sp{$scaf}}.=";$scaf($str,$length{$scaf})";
				}
				else{
					$temp_hash{$scaf2sp{$scaf}}="$scaf($str,$length{$scaf})";
				}
			}
		}
	}
	foreach(@clid){
		if($temp_hash{$_}){
			print OUT "\t$temp_hash{$_}";
		}
		else{
			print OUT "\tNA";
		}
	}
	print OUT "\n";
}
close TMP;
close OUT;
$current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
print LOG "$current_time: Produce $outdir/SCGscaffold_RawClassified.xls\n";
print STDERR "$current_time: Produce $outdir/SCGscaffold_RawClassified.xls\n";

##############################################################################  split  seeds ###########################################################
&splitseed;
%min=();
%max=();
%foldval=();
open(OF,">$outdir/allscaffoldSeed_ClassifiedLength_AfterSplit.xls")||die;
@tmptotalseed=();
my %classifiedFlag=();
%scaf2sp=();
@totalseed=sort {$seedmean{$b}{1} <=> $seedmean{$a}{1}} @totalseed;
foreach my $cl (@totalseed) {
	if($assign{$cl}){
		my @scaffold=split /;/,$assign{$cl};
		my $sumlen=0;
		my $seq="";
		my %scgnum=();
		foreach my $scaf (@scaffold){
			$sumlen+=$length{$scaf};
			$seq.=$seq{$scaf};
			if($scg{$scaf}){
				$scaf2sp{$scaf}=$cl;
				foreach my $scg (keys %{$scg{$scaf}}) {
					$scgnum{$scg}++;
				}
			}
		}
		my %spnum=();
		foreach my $scg (keys %scgnum) {
			my $tmpnum=$scgnum{$scg};
			$spnum{$tmpnum}++;
		}
		my $numsp=0;
		my @scgnum=sort {$b <=> $a} keys %spnum;
		foreach (@scgnum) {
			if($spnum{$_} >$scgNum){
				$numsp=$_;
				last;
			}
			else{
				next;
			}
		}
		my $flag=0;
		if($numsp >=1){
			$flag=1;
		}
		elsif($cl!~/;/ && $sumlen >=$LargeLen){
			$flag=1;
		}
		elsif($cl =~/;/ && $sumlen >=$SmallLen){
			$flag=1;
		}
		if($flag == 1){
			push @tmptotalseed,$cl;
			my $n=scalar @scaffold;
			print OF "$cl\t$sumlen\t$n\t$assign{$cl}\n";
		}
		else{
			foreach my $scaf (@scaffold) {
				if($length{$scaf} >=10000){
					push @remain_more10kb,$scaf;
				}
				else{
					push @remain_less10kb,$scaf;
				}
			}
			my $n=scalar @scaffold;
			print OF "$cl\t$sumlen\t$n\t$assign{$cl}\n";
			delete $assign{$cl};
		}
	}
	else{
		print OF "$cl\tNA\t0\t0\tNA\n";
	}
}
close OF;
$current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
print LOG "$current_time: Produce $outdir/allscaffoldSeed_ClassifiedLength_AfterSplit.xls\n";
print STDERR "$current_time: Produce $outdir/allscaffoldSeed_ClassifiedLength_AfterSplit.xls\n";

open(OUT,">$outdir/SCGscaffold_AfterSplit.xls")||die;
open(TMP,$tablefile)||die;
print OUT "SCG";
foreach(@totalseed){
	if(!$seedmean{$_}{1}){
		print STDERR "$_\n";
	}
}
foreach(@totalseed){
	print OUT "\t$_";
}
print OUT "\n";
while(<TMP>){
	chomp;
	my @a=split /\t/;
	print OUT "$a[0]";
	shift @a;
	my %temp_hash=();
	foreach(@a){
		if(/^([^\(]+)\(/){
			my $scaf=$1;
			my $str="";
			foreach my $index(@CovIndex){
				$str.=",$mean{$scaf}{$index}";
			}
			$str=~s/^\,//;
			if($scaf2sp{$scaf}){
				if($temp_hash{$scaf2sp{$scaf}}){
					$temp_hash{$scaf2sp{$scaf}}.=";$scaf($str,$length{$scaf})";
				}
				else{
					$temp_hash{$scaf2sp{$scaf}}="$scaf($str,$length{$scaf})";
				}
			}
		}
	}
	foreach(@totalseed){
		if($temp_hash{$_}){
			print OUT "\t$temp_hash{$_}";
		}
		else{
			print OUT "\tNA";	
		}
	}
	print OUT "\n";
}
close TMP;
close OUT;
$current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
print LOG "$current_time: Produce $outdir/SCGscaffold_AfterSplit.xls\n";
print STDERR "$current_time: Produce $outdir/SCGscaffold_AfterSplit.xls\n";
@totalseed=@tmptotalseed;
$n=scalar @totalseed;
print STDERR "$current_time: The number of seeds after filtering small bins: $n\n";
print LOG "$current_time: The number of seeds after filtering small bins: $n\n";

####################################################################  Merge  seeds  #################################################################
&mergeseed;
$n=scalar @totalseed;
$current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
print LOG "$current_time: The number of seeds after merging seeds at the first time: $n\n";
print STDERR "$current_time: The number of seeds after merging seeds at the first time: $n\n";

%scaf2sp=();
foreach my $cl(@totalseed){
	if($assign{$cl}){
		my @clscaf=split /;/,$assign{$cl};
		foreach my $scaf(@clscaf){
			$scaf2sp{$scaf}=$cl;
		}
	}
}
open(OUT,">$outdir/SCGscaffold_AfterMerge_FirstTime.xls")||die;
open(TMP,$tablefile)||die;
print OUT "SCG";
@totalseed=sort {$seedmean{$b}{1} <=> $seedmean{$a}{1}} @totalseed;
foreach(@totalseed){
	print OUT "\t$_";
}
print OUT "\n";
while(<TMP>){
	chomp;
	my @a=split /\t/;
	print OUT "$a[0]";
	shift @a;
	my %temp_hash=();
	foreach(@a){
		if(/^([^\(]+)\(/){
			my $scaf=$1;
			my $str="";
			foreach my $index(@CovIndex){
				$str.=",$mean{$scaf}{$index}";
			}
			$str=~s/^\,//;
			if($scaf2sp{$scaf}){
				if($temp_hash{$scaf2sp{$scaf}}){
					$temp_hash{$scaf2sp{$scaf}}.=";$scaf($str,$length{$scaf})";
				}
				else{
					$temp_hash{$scaf2sp{$scaf}}="$scaf($str,$length{$scaf})";
				}
			}
		}
	}
	foreach(@totalseed){
		if($temp_hash{$_}){
			print OUT "\t$temp_hash{$_}";
		}
		else{
			print OUT "\tNA";
		}
	}
	print OUT "\n";
}
close TMP;
close OUT;
$current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
print LOG "$current_time: Produce $outdir/SCGscaffold_AfterMerge_FirstTime.xls\n";
print STDERR "$current_time: Produce $outdir/SCGscaffold_AfterMerge_FirstTime.xls\n";

###########################################################  Merge  seeds by Non-SCG information ######################################################
&MergeDuplicatedSeed;
&Remergeseed;
$n=scalar @totalseed;
$current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
print LOG "$current_time: The number of seeds after merging seeds at the second time: $n\n";
print STDERR "$current_time: The number of seeds after merging seeds at the second time: $n\n";

%scaf2sp=();
foreach my $cl(keys %assign){
	my @clscaf=split /;/,$assign{$cl};
	foreach my $scaf(@clscaf){
		$scaf2sp{$scaf}=$cl;
	}
}
open(OUT,">$outdir/SCGscaffold_AfterMerge_SecondTime.xls")||die;
open(TMP,$tablefile)||die;
print OUT "SCG";
@totalseed=sort {$seedmean{$b}{1} <=> $seedmean{$a}{1}} @totalseed;
foreach(@totalseed){
	print OUT "\t$_";
}
print OUT "\n";
while(<TMP>){
	chomp;
	my @a=split /\t/;
	print OUT "$a[0]";
	shift @a;
	my %temp_hash=();
	foreach(@a){
		if(/^([^\(]+)\(/){
			my $scaf=$1;
			my $str="";
			foreach my $index(@CovIndex){
				$str.=",$mean{$scaf}{$index}";
			}
			$str=~s/^\,//;
			if($scaf2sp{$scaf}){
				if($temp_hash{$scaf2sp{$scaf}}){
					$temp_hash{$scaf2sp{$scaf}}.=";$scaf($str,$length{$scaf})";
				}
				else{
					$temp_hash{$scaf2sp{$scaf}}="$scaf($str,$length{$scaf})";
				}
			}
		}
	}
	foreach(@totalseed){
		if($temp_hash{$_}){
			print OUT "\t$temp_hash{$_}";
		}
		else{
			print OUT "\tNA";
		}
	}
	print OUT "\n";
}
close TMP;
close OUT;
$current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
print LOG "$current_time: Produce $outdir/SCGscaffold_AfterMerge_SecondTime.xls\n";
print STDERR "$current_time: Produce $outdir/SCGscaffold_AfterMerge_SecondTime.xls\n";

###########################################################  Reclassify scaffolds of small Bins ######################################################
@totalseed=sort {$seedmean{$b}{1} <=> $seedmean{$a}{1}} @totalseed;
my @largeseed=();
my @smallseed=();
open(OF,">$outdir/allscaffoldSeed_ClassifiedLength_Final.xls")||die;
foreach my $cl (@totalseed) {
	if($assign{$cl}){
		my %scgnum=();
		my $totallen=0;
		my @clscaf=split /;/,$assign{$cl};
		foreach my $scaf(@clscaf){
			$totallen+=$length{$scaf};
			if($scg{$scaf}){
				foreach my $scg (keys %{$scg{$scaf}}) {
					$scg2cl{$cl}{$scg}++;
					$scgnum{$scg}++;
				}
			}
		}
		my %spnum=();
		foreach my $scg (keys %scgnum) {
			my $tmpnum=$scgnum{$scg};
			$spnum{$tmpnum}++;
		}
		my $numsp=0;
		my @scgnum=sort {$b <=> $a} keys %spnum;
		foreach (@scgnum) {
			if($spnum{$_} >=$scgNum4bin){
				$numsp=$_;
				last;
			}
			else{
				next;
			}
		}
		if($totallen >=$BinLen || $numsp >=1){
			push @largeseed,$cl;
		}
		else{
			push @smallseed,$cl;
		}
		print OF "$cl\t$totallen\n";
	}
	else{
		print OF "$cl\t0\n";
	}
}
$current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
print STDERR "$current_time: Produce $outdir/allscaffoldSeed_ClassifiedLength_Final.xls\n";
print LOG "$current_time: Produce $outdir/allscaffoldSeed_ClassifiedLength_Final.xls\n";
@totalseed=@largeseed;
$n=scalar @totalseed;
$current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
print STDERR "$current_time: The number of seeds after filtering small bins is $n\n";
print LOG "$current_time: The number of seeds after filtering small bins is $n\n";

for(my $i=0;$i<=$#smallseed;$i++){
	my $cl1=$smallseed[$i];
	my @scaffold=split /;/,$assign{$cl1};
	foreach my $scaf (@scaffold) {
		if($length{$scaf} >=10000){
			push @remain_more10kb,$scaf;
		}
		else{
			push @remain_less10kb,$scaf;
		}
	}
	delete $assign{$cl1};
}

##############################################################################  Reclassify SCG scaffolds ###########################################################
my %tmpassign=();
foreach my $cl (@totalseed) {
	if($assign{$cl}){
		my @scaffold=split /;/,$assign{$cl};
		foreach my $scaf (@scaffold) {
			if($scg{$scaf}){
				if($used{$scaf}){
					if($tmpassign{$cl}){
						$tmpassign{$cl}.=";$scaf";
					}
					else{
						$tmpassign{$cl}=$scaf;
					}
				}
			}
			else{
				if($tmpassign{$cl}){
					$tmpassign{$cl}.=";$scaf";
				}
				else{
					$tmpassign{$cl}=$scaf;
				}
			}
		}
	}
}
%assign=%tmpassign;

%SCGassign=();
&PrepareForCorrelation;
foreach my $cl (@totalseed) {
	&Prepare4Bayesin($clseq{$cl},$cl,__LINE__);
}
&SCGscaffold_AllClassify(\@totalseed,\@totalSCGscaffold_more10kb,\@totalSCGscaffold_less10kb);
foreach my $cl (@totalseed) {
	if($SCGassign{$cl}){
		if($assign{$cl}){
			$assign{$cl}.=";$SCGassign{$cl}";
		}
		else{
			$assign{$cl}=$SCGassign{$cl};
		}
	}
}
%scaf2sp=();
foreach my $cl(@totalseed){
	if($assign{$cl}){
		my @clscaf=split /;/,$assign{$cl};
		foreach my $scaf(@clscaf){
			$scaf2sp{$scaf}=$cl;
		}
	}
}
open(OUT,">$outdir/SCGscaffold_FinalClassified.xls")||die;
print OUT "SCG";
@clid=sort {$seedmean{$b}{1} <=> $seedmean{$a}{1}} @totalseed;
foreach(@clid){
	print OUT "\t$_";
}
print OUT "\n";
open(TMP,$tablefile)||die;
while(<TMP>){
	chomp;
	my @a=split /\t/;
	print OUT "$a[0]";
	shift @a;
	my %temp_hash=();
	foreach(@a){
		if(/^([^\(]+)\(/){
			my $scaf=$1;
			my $str="";
			foreach my $index(@CovIndex){
				$str.=",$mean{$scaf}{$index}";
			}
			$str=~s/^\,//;
			if($scaf2sp{$scaf}){
				if($temp_hash{$scaf2sp{$scaf}}){
					$temp_hash{$scaf2sp{$scaf}}.=";$scaf($str,$length{$scaf})";
				}
				else{
					$temp_hash{$scaf2sp{$scaf}}="$scaf($str,$length{$scaf})";
				}
			}
		}
	}
	foreach(@clid){
		if($temp_hash{$_}){
			print OUT "\t$temp_hash{$_}";
		}
		else{
			print OUT "\tNA";
		}
	}
	print OUT "\n";
}
close TMP;
close OUT;
$current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
print LOG "$current_time: Produce $outdir/SCGscaffold_FinalClassified.xls\n";
print STDERR "$current_time: Produce $outdir/SCGscaffold_FinalClassified.xls\n";

###########################################################  Reclassify Non-SCG scaffolds ######################################################
foreach my $cl(@totalseed){
	my @scaffold=split /\;/,$assign{$cl};
	@scaffold=sort {$mean{$b}{1} <=> $mean{$a}{1}} @scaffold;
	my @SCGscaffold=();
	my $SCGseq="";
	my $totalseq="";
	my %hash;
	foreach my $scaf(@scaffold){
		$totalseq.=$seq{$scaf};
		if($scg{$scaf}){
			push @SCGscaffold,$scaf;
			$SCGseq.=$seq{$scaf};
			$hash{$scaf}=1;
		}
	}
	if(scalar @SCGscaffold >=10){
		&SpeciesCovInterval($cl,\@SCGscaffold);
		if(length $SCGseq >=10000){
			$clseq{$cl}=$SCGseq;
		}
		else{
			$clseq{$cl}=$SCGseq;
			my @scaffold=split /;/,$cl;
			foreach my $scaf (@scaffold) {
				if(!$hash{$scaf}){
					$clseq{$cl}.=$seq{$scaf};
				}
			}
		}
	}
	else{
		&CovInterval_4merge($cl,\@scaffold);
		$clseq{$cl}=$totalseq;
	}
}
&CovTreat;
open(OUT,">$outdir/FinalSeed_CovInterval.xls")||die;
foreach my $cl(@totalseed){
     print OUT "$cl";
     foreach my $index (@CovIndex) {
         print OUT "\t$seedmin{$cl}{$index}\t$seedmax{$cl}{$index}";
     }
     print OUT "\n";
}
close OUT;
$current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
print LOG "$current_time: Produce $outdir/FinalSeed_CovInterval.xls\n";
print STDERR "$current_time: Produce $outdir/FinalSeed_CovInterval.xls\n";

%tmpassign=();
foreach my $cl (@totalseed) {
	if($assign{$cl}){
		my @scaffold=split /;/,$assign{$cl};
		foreach my $scaf (@scaffold) {
			if($scg{$scaf}){
				if($tmpassign{$cl}){
					$tmpassign{$cl}.=";$scaf";
				}
				else{
					$tmpassign{$cl}=$scaf;
				}
			}
			else{
				if($used{$scaf}){
					if($tmpassign{$cl}){
						$tmpassign{$cl}.=";$scaf";
					}
					else{
						$tmpassign{$cl}=$scaf;
					}
				}
			}
		}
	}
}
%assign=%tmpassign;
%tmpassign=();

&PrepareForCorrelation;
$current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
print STDERR "$current_time: Start to reclassify large scaffolds\n";
print LOG "$current_time: Start to reclassify large scaffolds\n";
&Classify_LargeScaffold(\@totalseed,\@NonSCGscaffold_more10kb);
$current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
print STDERR "$current_time: Start to reclassify small scaffolds\n";
print LOG "$current_time: Start to reclassify small scaffolds\n";
$TAG=0;
&Classify_SmallScaffold(\@totalseed,\@NonSCGscaffold_less10kb,$TAG);
%prob=();

######################################################################## output finally classified scaffolds##############################################
open(OF,">$outdir/classified_allscaffold.xls")||die;
my $binid=0;
@clid=sort {$seedmean{$b}{1} <=> $seedmean{$a}{1}} @totalseed;
%classifiedFlag=();
foreach my $k(@clid){
	$binid++;
	my @scaffold=split /\;/,$assign{$k};
	@scaffold=sort {$mean{$b}{1} <=> $mean{$a}{1}} @scaffold;
	foreach my $scaf(@scaffold){
		$classifiedFlag{$scaf}=1;
		print OF "$scaf\t$binid\t$length{$scaf}";
		foreach my $index(@CovIndex){
			print OF "\t$mean{$scaf}{$index}";
		}
		print OF "\t$k\n";
	}
}
foreach my $scaf(keys %length) {
	next if($length{$scaf} <$minlen);
	if($classifiedFlag{$scaf}){
		next;
	}
	else{
		print OF "$scaf\tunassigned\t$length{$scaf}";
		foreach my $index (@CovIndex) {
			print OF "\t$mean{$scaf}{$index}";
		}
		print OF "\tunassigned\n";
	}
}
close OF;
$current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
print LOG "$current_time: Produce $outdir/classified_allscaffold.xls\n";
print STDERR "$current_time: Produce $outdir/classified_allscaffold.xls\n";
print STDERR "$current_time: Finish all\n";
close LOG;

######################################################################## Subroutine for getting large scaffold seeds ####################################################
sub largeScaffold_Seed{
	my ($ref1,$ref2,$tag)=@_;
	my @scaffold=@{$ref1};
	my @seed=@{$ref2};
	foreach my $scaf(@scaffold){
		if(scalar @seed == 0){
			push @seed,$scaf;
			foreach my $index (@CovIndex) {
				if(!$seedmean{$scaf}{$index}){
					$seedmean{$scaf}{$index}=$mean{$scaf}{$index};
				}
			}
			$fseed=$scaf;
		}
		else{
			$refxx=$refxx{$scaf};
			$len1=$len{$scaf};
			$SSxx=$SSxx{$scaf};
			my $flag=1;
			foreach my $cl(@seed){
				my $Num=0;
				my $cha=0;
				foreach my $index (@CovIndex) {
					if($seedmean{$cl}{$index} <5){
						$cha++;
					}
					else{
						if($mean{$scaf}{$index} >=$seedmin{$cl}{$index} && $mean{$scaf}{$index} <=$seedmax{$cl}{$index}){
							$Num++;
						}
					}
				}
				if($Num+$cha == scalar @CovIndex){
					my $refyy=$refxx{$cl};
					my $len2=$len{$cl};
					my $SSyy=$SSxx{$cl};
					my $SSxy=&SS($refxx,$refyy,__LINE__);
					my $cor=&correlation($SSxx,$SSyy,$SSxy);
					if($tag == 1 ){ 
						if($cor >=$diff{"$len1\t$len2"}){
							$flag=0;
							last;
						}
					}
					else{
						if($cor >=$cutoff{"$len1\t$len2"}){
							$flag=0;
							last;
						}
					}
				}
			}
			if($flag == 1){
				push @seed,$scaf;
				foreach my $index (@CovIndex) {
					if(!$seedmean{$scaf}{$index}){
						$seedmean{$scaf}{$index}=$mean{$scaf}{$index};
					}
				}				
			}
		}
	}
	return \@seed;
}

######################################################################## Subroutine for getting large scaffold seeds ####################################################
sub CompositionSeed{
	my ($ref1,$ref2,$tag)=@_;
	my @scaffold=@{$ref1};
	my @seed=@{$ref2};
	foreach my $scaf(@scaffold){
		if(scalar @seed == 0){
			push @seed,$scaf;
		}
		else{
			$refxx=$refxx{$scaf};
			$len1=$len{$scaf};
			$SSxx=$SSxx{$scaf};
			my $flag=1;
			foreach my $cl(@seed){
				my $refyy=$refxx{$cl};
				my $len2=$len{$cl};
				my $SSyy=$SSxx{$cl};
				my $SSxy=&SS($refxx,$refyy,__LINE__);
				my $cor=&correlation($SSxx,$SSyy,$SSxy);
				if($cor >=$cutoff{"$len1\t$len2"}){
					$flag=0;
					last;
				}
			}
			if($flag == 1){
				push @seed,$scaf;
			}
		}
	}
	return \@seed;
}

######################################################################## Subroutine for split seeds ####################################################
sub splitseed{
	$current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
	print LOG "$current_time: Start to split seeds\n";
	print STDERR "$current_time: Start to split seeds\n";
	my @tmptotalseed=();
	foreach my $cl(keys %assign){
		my @scaffold=split /;/,$assign{$cl};
		my $sumlen=0;
		my $seq="";
		my %scgnum=();
		my @NonSCGscaffold_more10kb_inside=();
		my @SCGscaffold_more10kb_inside=();
		my @NonSCGscaffold_less10kb_inside=();
		my @SCGscaffold_less10kb_inside=();
		my @more10kb_inside=();
		my @less10kb_inside=();
		my %scafOfscg=();
		my $str="";
		my @SCGscaffold=();
		foreach my $scaf (@scaffold){
			$sumlen+=$length{$scaf};
			if($scg{$scaf}){
				foreach my $scg (keys %{$scg{$scaf}}) {
					$scgnum{$scg}++;
					push @{$scafOfscg{$scg}},$scaf;
				}
			}
			if($used{$scaf}){
				$seq.=$seq{$scaf};
				if($str eq ""){
					$str=$scaf;
				}
				else{
					$str.=";$scaf";
				}
			}
			else{
				if($length{$scaf} <10000){
					if(!$scg{$scaf}){
						push @NonSCGscaffold_less10kb_inside,$scaf;
						push @less10kb_inside,$scaf;
					}
					else{
						push @SCGscaffold,$scaf;
					}
				}
				elsif($length{$scaf} >=10000){
					if(!$scg{$scaf}){
						push @NonSCGscaffold_more10kb_inside,$scaf;
						push @more10kb_inside,$scaf;
					}
					else{
						push @SCGscaffold,$scaf;
					}
				}
			}
		}
		my %spnum=();
		foreach my $scg (keys %scgnum) {
			my $tmpnum=$scgnum{$scg};
			$spnum{$tmpnum}++;
		}
		my $numsp=0;
		my @scgnum=sort {$b <=> $a} keys %spnum;
		foreach (@scgnum) {
			if($spnum{$_} >$scgNum){
				$numsp=$_;
				last;
			}
			else{
				next;
			}
		}
		if($numsp >1){
			if($str ne ""){
				$assign{$cl}=$str;
			}
			else{
				$assign{$cl}="";
			}
			my @tmparray=split /;/,$cl;
			foreach my $scaf (@tmparray) {
				if(!$used{$scaf}){
					$seq.=$seq{$scaf};
				}
			}
			my %tag=();
			foreach my $scg (keys %scgnum) {
				if($scgnum{$scg} >1){
					my @scaf=@{$scafOfscg{$scg}};
					foreach my $scaf (@scaf){
						$tag{$scaf}=1;
					}
				}
			}
			foreach my $scaf (@SCGscaffold) {
				if(!$tag{$scaf}){
					$seq.=$seq{$scaf};
					if($assign{$cl}){
						$assign{$cl}.=";$scaf";
					}
					else{
						$assign{$cl}=$scaf;
					}
				}
				else{
					if($length{$scaf} >=10000){
						push @SCGscaffold_more10kb_inside,$scaf;
						push @more10kb_inside,$scaf;
					}
					else{
						push @SCGscaffold_less10kb_inside,$scaf;
						push @less10kb_inside,$scaf;
					}
				}
			}
			&Prepare4Bayesin($seq,$cl);
			&Prepare4Zvalue($seq,$cl);
			my @scaffold=sort {$foldval{$cl}{$b} <=> $foldval{$cl}{$a}} @less10kb_inside;
			my @scafcomb=();
			for(my $i=0;$i<$#scaffold;$i++) {
				my $scaf1=$scaffold[$i];
				my $sumlen=$length{$scaf1};
				my $clstr=$scaf1;
				my $seq=$seq{$scaf1};
				my %sumdepth=();
				my %data=();
				foreach my $index(@CovIndex){
					$sumdepth{$index}=$mean{$scaf1}{$index}*$length{$scaf1};
					push @{$data{$index}},@{$depth{$scaf1}{$index}};
				}
				for(my $j=$i+1;$j<=$#scaffold;$j++) {
					my $scaf2=$scaffold[$j];
					my $Num=0;
					my $cha=0;
					foreach my $index (@CovIndex) {
					if($mean{$scaf1}{$index} >=$mean{$scaf2}{$index}){
						if($mean{$scaf1}{$index} <5){
								$cha++;
							}
							elsif($mean{$scaf2}{$index} >=$min{$scaf1}{$index} && $mean{$scaf2}{$index} <=$max{$scaf1}{$index}){
								$Num++;
							}
						}
						else{
							if($mean{$scaf2}{$index} <5){
								$cha++;
							}
							elsif($mean{$scaf1}{$index} >=$min{$scaf2}{$index} && $mean{$scaf1}{$index} <=$max{$scaf2}{$index}){
								$Num++;
							}
						}
					}
					if($Num+$cha == $IndexNum){
						my $dist=int(abs($foldval{$cl}{$scaf1}-$foldval{$cl}{$scaf2}));
						if($dist <=7){
							$sumlen+=$length{$scaf2};
							$seq.=$seq{$scaf2};
							$clstr.=";$scaf2";
							foreach my $index (@CovIndex) {
								$sumdepth{$index}+=$mean{$scaf2}{$index}*$length{$scaf2};
								push @{$data{$index}},@{$depth{$scaf2}{$index}};
							}
							if($sumlen >=10000){
								last;
							}
						}
						else{
							last;
						}
					}
					else{
						#print STDERR "$scaf1\t$scaf2\n";
					}
				}
				if($sumlen >=10000){
					foreach my $index (@CovIndex) {
						$seedmean{$clstr}{$index}=int($sumdepth{$index}/$sumlen);
						$seedmin{$clstr}{$index}=$seedmin{$cl}{$index};
						$seedmax{$clstr}{$index}=$seedmax{$cl}{$index};
					}
					&Prepare4Zvalue($seq,$clstr);
					push @scafcomb,$clstr;
					$seq{$clstr}=$seq;
				}
			}
			##############form seeds #################
			my @scafcomb_largescaf=(@scafcomb,@more10kb_inside);
			foreach my $scaf (@more10kb_inside) {
				foreach my $index (@CovIndex) {
					$seedmin{$scaf}{$index}=$seedmin{$cl}{$index};
					$seedmax{$scaf}{$index}=$seedmax{$cl}{$index};
					$seedmean{$scaf}{$index}=$mean{$scaf}{$index};
				}
			}
			my @possible_seed=($cl);
			my $num=1;
			my %pvalue=();
			while($num <$numsp){
				my %totalpvalue=();
				foreach my $cl1 (@scafcomb_largescaf) {
					my $refxx=$refxx{$cl1};
					my $len=$len{$cl1};
					my $SSxx=$SSxx{$cl1};
					my $totalpvalue=0;
					foreach my $cl2 (@possible_seed) {
						if($pvalue{"$cl1\t$cl2"}){
							$totalpvalue+=$pvalue{"$cl1\t$cl2"};
						}
						elsif($pvalue{"$cl2\t$cl1"}){
							$totalpvalue+=$pvalue{"$cl2\t$cl1"};
						}
						else{
							my $refyy=$refxx{$cl2};
							my $len2=$len{$cl2};
							my $SSyy=$SSxx{$cl2};
							my $SSxy=&SS($refxx,$refyy,__LINE__);
							my $cor=&correlation($SSxx,$SSyy,$SSxy);
							my $same_mean=$smean{"$len\t$len2"};
							my $same_sd=$ssd{"$len\t$len2"};
							my $val=($cor-$same_mean)/$same_sd;
							my $pvalue=uprob($val);
							$pvalue=1-$pvalue;
							$pvalue{"$cl1\t$cl2"}=$pvalue;
							$totalpvalue+=$pvalue;
						}
					}
					$totalpvalue{$cl1}=$totalpvalue;
				}
				@scafcomb_largescaf= sort {$totalpvalue{$a} <=> $totalpvalue{$b}} keys %totalpvalue;
				my $seed=shift @scafcomb_largescaf;
				push @possible_seed,$seed;
				$assign{$seed}="";
				&Prepare4Bayesin($seq{$seed},$seed,__LINE__);
				$num++;
				if(scalar @scafcomb_largescaf ==0){
					last;
				}
			}
			%SCGassign=();
			&SCGscaffold_AllClassify(\@possible_seed,\@SCGscaffold_more10kb_inside,\@SCGscaffold_less10kb_inside);
			%prob=();
			foreach my $seed (@possible_seed) {
				if($SCGassign{$seed}){
					@scaffold=split /;/,$SCGassign{$seed};
					if($assign{$seed}){
						$assign{$seed}.=";$SCGassign{$seed}";
					}
					else{
						$assign{$seed}=$SCGassign{$seed};
					}
				}
				my %hash=();
				my $seq="";
				@scaffold=split /;/,$assign{$seed};
				foreach my $scaf (@scaffold) {
					$seq.=$seq{$scaf};
					$hash{$scaf}=1;
				}
				my $totallen=length $seq;
				if($totallen <10000){
					my @clscaf=split /;/,$seed;
					foreach my $scaf (@clscaf) {
						if(!$hash{$scaf}){
							$seq.=$seq{$scaf};
						}
					}
				}
				&Prepare4Zvalue($seq,$seed);
				$clseq{$seed}=$seq;
			}
			&Classify_LargeScaffold(\@possible_seed,\@NonSCGscaffold_more10kb_inside);
			my $TAG=0;
			&Classify_SmallScaffold(\@possible_seed,\@NonSCGscaffold_less10kb_inside,$TAG);
			push @tmptotalseed,@possible_seed;
		}
		else{
			push @tmptotalseed,$cl;
		}
	}
	@totalseed=@tmptotalseed;
	$n=scalar @totalseed;
	$current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
	print LOG "$current_time: The number of seeds after splitting seeds: $n\n";
	print STDERR "$current_time: The number of seeds after splitting seeds: $n\n";	
}

################################################################################### Subroutine for Cov treatment ####################################################
sub CovTreat{
	foreach my $index (@CovIndex) {
		my @meanid=sort {$mean{$a}{$index} <=>$mean{$b}{$index}} keys %mean;
		my @tmpID=sort {$seedmin{$a}{$index} <=> $seedmin{$b}{$index}} @totalseed;
		if($seedmin{$tmpID[0]}{$index} >$mean{$meanid[0]}{$index}){
			$seedmin{$tmpID[0]}{$index}=int($mean{$meanid[0]}{$index});
		}
	}
	foreach my $index (@CovIndex) {
		my @depth=();
		foreach my $cl (@totalseed) {
			push @depth,$seedmin{$cl}{$index};
			push @depth,$seedmax{$cl}{$index};
		}
		@depth=sort {$a <=> $b} @depth;
		my %num;
		for(my $i=0;$i<$#depth;$i++){
			my $min=$depth[$i];
			my $max=$depth[$i+1];
			foreach my $cl (@totalseed) {
				if($seedmin{$cl}{$index} <=$min && $seedmax{$cl}{$index} >=$max){
					$num{"$min\t$max"}++;
				}
			}
		}
		for(my $i=0;$i<$#depth;$i++){
			my $min=$depth[$i];
			my $max=$depth[$i+1];
			if(!$num{"$min\t$max"}){
				foreach my $cl (@totalseed) {
					if($seedmin{$cl}{$index} == $max){
						$seedmin{$cl}{$index}=$min;
					}
					if($seedmax{$cl}{$index} == $min){
						$seedmax{$cl}{$index}=$max;
					}
				}
			}
		}
	}
}

################################################################################### Subroutine Cov interval computation ####################################################
sub CovInterval{
	my ($cl,$ref)=@_;
	my ($tmax,$tmin,$sd);
	my @scaffold=@{$ref};
	my %data=();
	my %sumdepth=();
	my $sumlen=0;
	foreach my $scaf(@scaffold){
		$sumlen+=$length{$scaf};
		foreach my $index (@CovIndex) {
			push @{$data{$index}},$mean{$scaf}{$index};
			$sumdepth{$index}+=$mean{$scaf}{$index}*$length{$scaf};
		}
	}
	foreach my $index (@CovIndex) {
		my $mean=sprintf("%.2f",$sumdepth{$index}/$sumlen);
		$seedmean{$cl}{$index}=$mean;
		if($seedmean{$cl}{$index} >=5){
			my @data=@{$data{$index}};
			@data=sort {$a <=> $b} @data;
			if(scalar @data >= 3){
				($tmin,$tmax)=&depth_process(\@data);
				my @tmpdata=();
				if(scalar @tmpdata >=2){
					$sd=&sd(\@tmpdata,$mean);
					$sd /=$mean;
					if($sd >$sd2{$index}){
						$sd=$sd2{$index};
					}
					elsif($sd <$sd1{$index}){
						$sd=$sd1{$index};
					}
				}
				else{
					$sd=$sd1{$index};
				}
			}
			else{
				$sd=$sd1{$index};
			}
			$seedmax{$cl}{$index}=int($seedmean{$cl}{$index}*(1+$sd*$fold))+1;
			$seedmin{$cl}{$index}=int($seedmean{$cl}{$index}*(1-$sd*$fold));
		}
		else{
			$seedmax{$cl}{$index}=0;
			$seedmin{$cl}{$index}=0;
		}
	}
}

################################################################################### Subroutine for species Cov interval ####################################################
sub SpeciesCovInterval{
	my ($cl,$ref)=@_;
	my @scaffold=@{$ref};
	my %data=();
	foreach my $scaf(@scaffold){
		foreach my $index (@CovIndex) {
			push @{$data{$index}},$mean{$scaf}{$index};
		}
	}
	foreach my $index (@CovIndex) {
		my @data=@{$data{$index}};
		my ($tmin,$tmax)=&depth_process(\@data);
		my @tmpdata=();
		my $sumdepth=0;
		my $sumlen=0;
		foreach my $scaf (@scaffold) {
			if($mean{$scaf}{$index} >=$tmin && $mean{$scaf}{$index} <=$tmax){
				$sumdepth+=$mean{$scaf}{$index}*$length{$scaf};
				$sumlen+=$length{$scaf};
				#push @tmpdata,$mean{$scaf}{$index};
			}
		}
		foreach(@data){
			if($_ >=$tmin && $_ <=$tmax){
				push @tmpdata,$_;
			}
		}
		my $mean=sprintf("%.2f",$sumdepth/$sumlen);
		$seedmean{$cl}{$index}=$mean;
		if($seedmean{$cl}{$index} >=5){
			@tmpdata=sort {$a<=>$b} @tmpdata;
			my $sd=&sd(\@tmpdata,$mean);
			$sd /=$mean;
			if(scalar @scaffold <20){
				my $tmpsd=$sd;
				if($tmpsd >$sd2{$index}){
					$tmpsd=$sd2{$index};
				}
				elsif($tmpsd <$sd1{$index}){
					$tmpsd=$sd1{$index};
					$sd=$sd1{$index};
				}
				my $max=int($seedmean{$cl}{$index}*(1+$tmpsd*$fold))+1;
				if($max >$tmpdata[-1]){
					$tmpdata[-1]=$max;
				}
				my $min=int($seedmean{$cl}{$index}*(1-$tmpsd*$fold));
				if($min <$tmpdata[0]){
					$tmpdata[0]=$min;
				}
			}
			$seedmin{$cl}{$index}=int($tmpdata[0]*(1-$sd));
			$seedmax{$cl}{$index}=int($tmpdata[-1]*(1+$sd))+1;
		}
		else{
			$seedmin{$cl}{$index}=0;
			$seedmax{$cl}{$index}=0;
		}
	}
}

################################################################################### Subroutine for Calculate Depth Intervals ####################################################
sub CovInterval_4merge{
	my ($cl,$ref)=@_;
	my @scaffold=@{$ref};
	my $sd;
	my %data=();
	foreach my $scaf(@scaffold){
		foreach my $index (@CovIndex) {
			push @{$data{$index}},$mean{$scaf}{$index};
		}
	}
	foreach my $index (@CovIndex) {
		my @data=@{$data{$index}};
		@data=sort {$a <=> $b} @data;
		my $sumdepth=0;
		my $sumlen=0;
		foreach my $scaf (@scaffold) {
			$sumdepth+=$mean{$scaf}{$index}*$length{$scaf};
			$sumlen+=$length{$scaf};
		}
		my $mean=sprintf("%.2f",($sumdepth/$sumlen));
		if(scalar @data >=2){
			$sd=&sd(\@data,$mean);
		}
		else{
			$sd=$sd1{$index}*$mean;
		}
		$seedmean{$cl}{$index}=$mean;
		$seedmax{$cl}{$index}=int($data[-1])+1;
		$seedmin{$cl}{$index}=int($data[0]);
	}
}

################################################################################### Subroutine for multiple-thread Bayesin####################################################
sub MultipleThreads_AllBayesin{
	my ($cl,$ref,$line)=@_;
	if(!$clseq{$cl}){
		print STDERR "$line\n";
	}
	&Prepare4Bayesin($clseq{$cl},$cl,__LINE__);
	my %tmpresult=();
	my @smallscaffold=@{$ref};
	foreach my $scaf (@smallscaffold) {
		my $len=$length{$scaf};
		my $val=0;
		foreach my $sub (keys %{$freq{$scaf}}) {
			if($prob{$cl}{$sub}){
				$val +=$freq{$scaf}{$sub}*$prob{$cl}{$sub};
			}
			else{
				print STDERR "$line\n";
			}
		}
		$tmpresult{$scaf}{$cl}=$val*500/$length{$scaf};
	}
	delete $prob{$cl};
	return (\%tmpresult);
}

################################################################################### Subroutine for multiple-thread Bayesin####################################################
sub foldval{
	my ($cl,$ref,$line)=@_;
	&Prepare4Bayesin($clseq{$cl},$cl,__LINE__);
	my %tmpfoldval=();
	my @smallscaffold=@{$ref};
	foreach my $scaf (@smallscaffold) {
		my $len=$length{$scaf};
		my $val=0;
		foreach my $sub (keys %{$freq{$scaf}}) {
			if($prob{$cl}{$sub}){
				$val +=$freq{$scaf}{$sub}*$prob{$cl}{$sub};
			}
			else{
				print STDERR "$line\n";
			}
		}
		$tmpfoldval{$cl}{$scaf}=$val*500/$length{$scaf};
	}
	delete $prob{$cl};
	return \%tmpfoldval;
}

######################################################### Subroutine for bayesin analysis #######################################################################
sub bayesin{
	my ($scaf,$reference,$REF,$line)=@_;
	my %prob=%{$REF};
	my $len=$length{$scaf};
	my %val=();
	my @refID=@{$reference};
	foreach my $cl (@refID) {
		my $val=0;
		foreach my $sub (keys %{$freq{$scaf}}) {
			if($prob{$cl}{$sub}){
				$val +=$freq{$scaf}{$sub}*$prob{$cl}{$sub};
			}
			else{
				print STDERR "$line\n";
			}
		}
		$val{$cl}=$val*500/$length{$scaf};
	}
	@refID=sort{$val{$b} <=> $val{$a}} @refID;
	my $clstr=shift @refID;
	my $fclstr=$clstr;
	foreach my $cl (@refID) {
		my $tmp=$val{$fclstr}-$val{$cl};
		if($tmp >=7){
			last;
		}
		else{
			$clstr.="\#cl";
		}
	}
	return($clstr,$val{$fclstr});
}

################################################################################### classify small scaffolds ####################################################
sub Classify_SmallScaffold{
	my ($ref1,$ref2,$TAG)=@_;
	#When $TAG =1, remain foldvals,otherwise, do not  remain foldvals
	my @seed=@{$ref1};
	my @smallscaffold=@{$ref2};
	my %scaffold=();
	my @tmpseed=();
	foreach my $scaf (@smallscaffold) {
		my $tag=0;
		foreach my $cl (@seed) {
			my $num=0;
			my $cha=0;
			foreach my $index (@CovIndex) {
				if($seedmean{$cl}{$index} <5){
					$cha++;
				}
				elsif($mean{$scaf}{$index} >=$seedmin{$cl}{$index} && $mean{$scaf}{$index} <$seedmax{$cl}{$index}){
					$num++;
				}
			}
			if($num+$cha == $IndexNum){
				$tag++;
				push @{$scaffold{$cl}},$scaf;
			}
		}
		if($tag == 0){
			push @{$scaffold{$fseed}},$scaf;
		}
	}
	foreach my $cl (@seed) {
		if($scaffold{$cl}){
			push @tmpseed,$cl;
		}
	}
	@seed=@tmpseed;
	%result=();
	$IndexStart=0;
	@running=();
	@Threads=();
	%joined=();
	while($IndexStart <=$#seed){
		while(scalar @running< $nb_process && $IndexStart <=$#seed){
			my $seed=$seed[$IndexStart];
			$IndexStart ++;
			next if(!$scaffold{$seed});
			my @scaffold=@{$scaffold{$seed}};
			my $thread = threads->create({'context' => 'list'},\&MultipleThreads_AllBayesin,$seed,\@scaffold,__LINE__);
			push @Threads, $thread;
			foreach my $thr (@Threads) {
				if ($thr->is_joinable()) {
					my ($ref) = $thr->join();
					my %tmpfoldval=%{$ref};
					foreach my $scaf (keys %tmpfoldval) {
						foreach my $cl (keys %{$tmpfoldval{$scaf}}) {
							$result{$scaf}{$cl}=$tmpfoldval{$scaf}{$cl};
						}
					}
					$joined{$thr}=1;
				}
				@running = threads->list(threads::running);
			}
		}
		@running = threads->list(threads::running);
	}
	while (scalar @running != 0) {
		foreach my $thr (@Threads) {
			if ($thr->is_joinable()){
				my ($ref) = $thr->join();
				my %tmpfoldval=%{$ref};
				foreach my $scaf (keys %tmpfoldval) {
					foreach my $cl (keys %{$tmpfoldval{$scaf}}) {
						$result{$scaf}{$cl}=$tmpfoldval{$scaf}{$cl};
					}
				}
				$joined{$thr}=1;
			}
		}
		@running = threads->list(threads::running);
	}
	foreach my $thr (@Threads){
		if($joined{$thr}){
			next;
		}
		else{
			if($thr->is_joinable()){
				my ($ref) = $thr->join();
				my %tmpfoldval=%{$ref};
				foreach my $scaf (keys %tmpfoldval) {
					foreach my $cl (keys %{$tmpfoldval{$scaf}}) {
						$result{$scaf}{$cl}=$tmpfoldval{$scaf}{$cl};
					}
				}
				$joined{$thr}=1;
			}
			else{
				my $tid = $thr->tid;
				print STDERR "$tid is not joined\n";
			}
		}
	}
	%scaffold=();
	foreach my $scaf (@smallscaffold) {
		if($result{$scaf}){
			my @source=sort {$result{$scaf}{$b} <=> $result{$scaf}{$a}} keys %{$result{$scaf}};
			my $cl=$source[0];
			if($assign{$cl}){
				$assign{$cl}.=";$scaf";
			}
			else{
				$assign{$cl}="$scaf";
			}
			if($TAG == 1){
				$foldval{$cl}{$scaf}=$result{$scaf}{$cl};
			}
		}
	}
	%result=();
}

################################################################################### classify small scaffolds ####################################################
sub Classify_SmallScaffold_withoutCov{
	my ($ref1,$ref2)=@_;
	my @seed=@{$ref1};
	my @smallscaffold=@{$ref2};
	%result=();
	$IndexStart=0;
	@running=();
	@Threads=();
	%joined=();
	while($IndexStart <=$#seed){
		while(scalar @running< $nb_process && $IndexStart <=$#seed){
			my $seed=$seed[$IndexStart];
			$IndexStart ++;
			next if($largeScaffold_compositionSeed{$seed});
			my $thread = threads->create({'context' => 'list'},\&MultipleThreads_AllBayesin,$seed,\@smallscaffold,__LINE__);
			push @Threads, $thread;
			foreach my $thr (@Threads) {
				if ($thr->is_joinable()) {
					my ($ref) = $thr->join();
					my %tmpfoldval=%{$ref};
					foreach my $scaf (keys %tmpfoldval) {
						foreach my $cl (keys %{$tmpfoldval{$scaf}}) {
							$result{$scaf}{$cl}=$tmpfoldval{$scaf}{$cl};
						}
					}
					$joined{$thr}=1;
				}
				@running = threads->list(threads::running);
			}
		}
		@running = threads->list(threads::running);
	}
	while (scalar @running != 0) {
		foreach my $thr (@Threads) {
			if ($thr->is_joinable()){
				my ($ref) = $thr->join();
				my %tmpfoldval=%{$ref};
				foreach my $scaf (keys %tmpfoldval) {
					foreach my $cl (keys %{$tmpfoldval{$scaf}}) {
						$result{$scaf}{$cl}=$tmpfoldval{$scaf}{$cl};
					}
				}
				$joined{$thr}=1;
			}
		}
		@running = threads->list(threads::running);
	}
	foreach my $thr (@Threads){
		if($joined{$thr}){
			next;
		}
		else{
			if($thr->is_joinable()){
				my ($ref) = $thr->join();
				my %tmpfoldval=%{$ref};
				foreach my $scaf (keys %tmpfoldval) {
					foreach my $cl (keys %{$tmpfoldval{$scaf}}) {
						$result{$scaf}{$cl}=$tmpfoldval{$scaf}{$cl};
					}
				}
				$joined{$thr}=1;
			}
			else{
				my $tid = $thr->tid;
				print STDERR "$tid is not joined\n";
			}
		}
	}
	foreach my $scaf (@smallscaffold) {
		if($result{$scaf}){
			my @source=sort {$result{$scaf}{$b} <=> $result{$scaf}{$a}} keys %{$result{$scaf}};
			if($maxfoldval{$scaf}){
				my @tmpsource=keys %{$maxfoldval{$scaf}};
				my $cl=$tmpsource[0];
				if($result{$scaf}{$source[0]} <$maxfoldval{$scaf}{$cl}){
					if($assign{$cl}){
						$assign{$cl}.=";$scaf";
					}
					else{
						$assign{$cl}="$scaf";
					}
					$foldval{$cl}{$scaf}=$maxfoldval{$scaf}{$cl};
				}
				else{
					my $cl=$source[0];
					if($assign{$cl}){
						$assign{$cl}.=";$scaf";
					}
					else{
						$assign{$cl}="$scaf";
					}
					$foldval{$cl}{$scaf}=$result{$scaf}{$cl};
				}
			}
			else{
				my $cl=$source[0];
				if($assign{$cl}){
					$assign{$cl}.=";$scaf";
				}
				else{
					$assign{$cl}="$scaf";
				}
				$foldval{$cl}{$scaf}=$result{$scaf}{$cl};
			}
		}
	}
	%result=();
}

################################################################################### Merge classification ####################################################
sub mergeseed{
	my %scg2cl=();
	my @duplicatedSeed=();
	my @singleSeed=();
	%clseq=();
	foreach my $cl(@totalseed){
		my @scaffold=split /;/,$cl;
		my %hash;
		foreach my $scaf (@scaffold) {
			$clseq{$cl}.=$seq{$scaf};
			$hash{$scaf}=1;
		}
		if($assign{$cl}){
			my @clscaf=split /;/,$assign{$cl};
			my %scgnum_inside=();
			foreach my $scaf(@clscaf){
				if($scg{$scaf}){
					foreach my $scg (keys %{$scg{$scaf}}) {
						$scg2cl{$cl}{$scg}++;
						$scgnum_inside{$scg}++;
					}
				}
			}
			my %spnum_inside=();
			foreach my $scg (keys %scgnum_inside) {
				my $tmpnum=$scgnum_inside{$scg};
				$spnum_inside{$tmpnum}++;
			}
			my $numsp_inside=0;
			my @scgnum_inside=sort {$b <=> $a} keys %spnum_inside;
			foreach (@scgnum_inside) {
				if($spnum_inside{$_} >$scgNum){
					$numsp_inside=$_;
					last;
				}
				else{
					next;
				}
			}
			if($numsp_inside >1){
				push @duplicatedSeed,$cl;
			}
			else{
				push @singleSeed,$cl;
				foreach my $scaf (@clscaf) {
					if(!$hash{$scaf}){
						$clseq{$cl}.=$seq{$scaf};
					}
				}
			}
		}
	}	
	foreach my $cl (@singleSeed){
		my @scaffold=split /;/,$assign{$cl};
		&CovInterval_4merge($cl,\@scaffold);
	}
	&PrepareForCorrelation;
	my %merge=();
	while(){
		foreach my $cl (@singleSeed) {
			if(!$seedmean{$cl}{1}){
				print STDERR "$cl have no mean cov\n";
			}
		}
		my @clid=sort {$seedmean{$b}{1} <=> $seedmean{$a}{1}} @singleSeed;
		my $seednum=scalar @singleSeed;
		my @tmptotalseed=();
		for(my $i=0;$i<=$#clid;$i++){
			my $cl1=$clid[$i];
			next if($merge{$cl1});
			my $refxx=$refxx{$cl1};
			my $SSxx=$SSxx{$cl1};
			my $len1=$len{$cl1};
			while(){
				my %pvalue=();
				my @source=();
				for(my $j=$i+1;$j<=$#clid;$j++){
					my $cl2=$clid[$j];
					next if($merge{$cl2});
					next if($cl1 eq $cl2);
					my $refyy=$refxx{$cl2};
					my $SSyy=$SSxx{$cl2};
					my $len2=$len{$cl2};
					my $SSxy=&SS($refxx,$refyy,__LINE__);
					my $cor=&correlation($SSxx,$SSyy,$SSxy);
					my $cutoff=$cutoff{"$len1\t$len2"};
					if($cor >=$cutoff){
						my %allscg=();
						if($scg2cl{$cl1}){
							foreach my $scg(keys %{$scg2cl{$cl1}}) {
								$allscg{$scg}+=$scg2cl{$cl1}{$scg};
							}
						}
						if($scg2cl{$cl2}){
							foreach my $scg(keys %{$scg2cl{$cl2}}) {
								$allscg{$scg}+=$scg2cl{$cl2}{$scg};
							}
						}
						my $dup=0;
						foreach my $scg (keys %allscg) {
							if($allscg{$scg} >1){
								$dup++;
							}
						}
						if($dup <=$scgNum){
							push @source,$cl2;
							my $same_mean=$smean{"$len1\t$len2"};
							my $same_sd=$ssd{"$len1\t$len2"};
							my $val=($cor-$same_mean)/$same_sd;
							my $pvalue=uprob($val);
							$pvalue=1-$pvalue;
							$pvalue{$cl2}=$pvalue;
						}
						else{
							next;
						}
					}
					else{
						next;
					}
				}
				if(scalar @source >0){
					my @tmpsource=();
					my @tmpseed=();
					foreach my $cl2 (@source) {
						my $Num=0;
						my $cha=0;
						my $tmpnum=0;
						foreach my $index (@CovIndex) {
							if($seedmean{$cl1}{$index} <5){
								$cha++;
							}
							else{
								if($seedmean{$cl2}{$index} <=$seedmax{$cl1}{$index} && $seedmean{$cl2}{$index} >=$seedmin{$cl1}{$index}){
									$Num++;
								}
								if($seedmin{$cl2}{$index} <$seedmax{$cl1}{$index} && $seedmax{$cl2}{$index} >$seedmin{$cl1}{$index}){
									$tmpnum++;
								}
							}
						}
						if($Num+$cha == $IndexNum){
							push @tmpsource,$cl2;
						}
						if($tmpnum+$cha == $IndexNum){
							push @tmpseed,$cl2;
						}
					}
					my $seed="";
					if(scalar @tmpsource >0){
						@tmpsource=sort {$pvalue{$b} <=> $pvalue{$a}} @tmpsource;
						$seed=$tmpsource[0];
					}
					elsif(scalar @tmpseed >0){
						my %dist=();
						foreach my $cl2 (@tmpseed) {
							my $sum=0;
							foreach my $index (@CovIndex) {
								$sum+=($seedmean{$cl1}{$index}-$seedmean{$cl2}{$index})**2;
							}
							my $dist=$sum**0.5;
							$dist{$cl2}=$dist;
						}
						@tmpseed=sort {$dist{$b} <=> $dist{$a}} @tmpseed;
						$seed=$tmpseed[0];
					}
					if($seed ne ""){
						$merge{$seed}=1;
						$clseq{$cl1}.=$clseq{$seed};
						my @tmpscaffold=();
						if($assign{$cl1}){
							push @tmpscaffold,split /;/,$assign{$cl1};
						}
						if($assign{$seed}){
							push @tmpscaffold,split /;/,$assign{$seed};
						}
						$assign{$cl1}=join(";",@tmpscaffold);
						&CovInterval_4merge($cl1,\@tmpscaffold);
						&Prepare4Zvalue($clseq{$cl1},$cl1);
						if($scg2cl{$seed}){
							foreach my $scg(keys %{$scg2cl{$seed}}) {
								if($scg2cl{$cl1}{$scg}){
									$scg2cl{$cl1}{$scg}+=$scg2cl{$seed}{$scg};
								}
								else{
									$scg2cl{$cl1}{$scg}=$scg2cl{$seed}{$scg};
								}
							}
						}
						if($scg2cl{$seed}){
							delete $scg2cl{$seed};
						}
						delete $assign{$seed};
						delete $clseq{$seed};
					}
					else{
						last;
					}
				}
				else{
					last;
				}
			}
			push @tmptotalseed,$cl1;	
		}
		@singleSeed=@tmptotalseed;
		if($seednum == scalar @singleSeed){
			last;
		}
	}
	@totalseed=(@singleSeed,@duplicatedSeed);
}

################################################################################### Merge classification ####################################################
sub MergeDuplicatedSeed{
	my %scg2cl=();
	%clseq=();
	foreach my $cl(@totalseed){
		my %hash=();
		my $SCGseq="";
		my $totalseq="";
		my @SCGscaffold=();
		if($assign{$cl}){
			my @clscaf=split /;/,$assign{$cl};
			foreach my $scaf(@clscaf){
				$totalseq.=$seq{$scaf};
				if($scg{$scaf}){
					$SCGseq.=$seq{$scaf};
					$hash{$scaf}=1;
					push @SCGscaffold,$scaf;
					foreach my $scg (keys %{$scg{$scaf}}) {
						$scg2cl{$cl}{$scg}++;
					}
				}
			}
			if(scalar @SCGscaffold >=10){
				&SpeciesCovInterval($cl,\@SCGscaffold);
				if(length $SCGseq >=10000){
					$clseq{$cl}=$SCGseq;
				}
				else{
					$clseq{$cl}=$SCGseq;
					my @scaffold=split /;/,$cl;
					foreach my $scaf (@scaffold) {
						if(!$hash{$scaf}){
							$clseq{$cl}.=$seq{$scaf};
						}
					}
				}
			}
			else{
				&CovInterval_4merge($cl,\@clscaf);
				$clseq{$cl}=$totalseq;
			}
		}
	}
	&PrepareForCorrelation;
	while(){
		my @clid=sort {$seedmean{$b}{1} <=> $seedmean{$a}{1}} @totalseed;
		my $seednum=scalar @totalseed;
		my @tmptotalseed=();
		my %merge=();
		for(my $i=0;$i<=$#clid;$i++){
			my $cl1=$clid[$i];
			next if($merge{$cl1});
			my $refxx=$refxx{$cl1};
			my $SSxx=$SSxx{$cl1};
			my $len1=$len{$cl1};
			while(){
				my %pvalue=();
				my @source=();
				for(my $j=$i+1;$j<=$#clid;$j++){
					my $cl2=$clid[$j];
					next if($merge{$cl2});
					next if($cl1 eq $cl2);
					my $refyy=$refxx{$cl2};
					my $SSyy=$SSxx{$cl2};
					my $len2=$len{$cl2};
					my $SSxy=&SS($refxx,$refyy,__LINE__);
					my $cor=&correlation($SSxx,$SSyy,$SSxy);
					my $cutoff=$cutoff{"$len1\t$len2"};
					if($cor >=$cutoff){
						my %allscg=();
						if($scg2cl{$cl1}){
							foreach my $scg(keys %{$scg2cl{$cl1}}) {
								$allscg{$scg}+=$scg2cl{$cl1}{$scg};
							}
						}
						if($scg2cl{$cl2}){
							foreach my $scg(keys %{$scg2cl{$cl2}}) {
								$allscg{$scg}+=$scg2cl{$cl2}{$scg};
							}
						}
						my $dup=0;
						my $totaldup=0;
						foreach my $scg (keys %allscg) {
							if($allscg{$scg} >1){
								$totaldup++;
							}
							if($scg2cl{$cl1} && $scg2cl{$cl1}{$scg} && $scg2cl{$cl2} && $scg2cl{$cl2}{$scg}){
								$dup++;
							}
						}
						if($dup <10){
							push @source,$cl2;
							my $same_mean=$smean{"$len1\t$len2"};
							my $same_sd=$ssd{"$len1\t$len2"};
							my $val=($cor-$same_mean)/$same_sd;
							my $pvalue=uprob($val);
							$pvalue=1-$pvalue;
							$pvalue{$cl2}=$pvalue;
						}
						else{
							next;
						}
					}
					else{
						next;
					}
				}
				if(scalar @source >0){
					my @tmpsource=();
					my @tmpseed=();
					foreach my $cl2 (@source) {
						my $Num=0;
						my $cha=0;
						my $tmpnum=0;
						foreach my $index (@CovIndex) {
							if($seedmean{$cl1}{$index} <5){
								$cha++;
							}
							else{
								if($seedmean{$cl2}{$index} <=$seedmax{$cl1}{$index} && $seedmean{$cl2}{$index} >=$seedmin{$cl1}{$index}){
									$Num++;
								}
								if($seedmin{$cl2}{$index} <$seedmax{$cl1}{$index} && $seedmax{$cl2}{$index} >$seedmin{$cl1}{$index}){
									$tmpnum++;
								}
							}
						}
						if($Num+$cha == $IndexNum){
							push @tmpsource,$cl2;
						}
						if($tmpnum+$cha == $IndexNum){
							push @tmpseed,$cl2;
						}
					}
					my $seed="";
					if(scalar @tmpsource >0){
						@tmpsource=sort {$pvalue{$b} <=> $pvalue{$a}} @tmpsource;
						$seed=$tmpsource[0];
					}
					elsif(scalar @tmpseed >0){
						my %dist=();
						foreach my $cl2 (@tmpseed) {
							my $sum=0;
							foreach my $index (@CovIndex) {
								$sum+=($seedmean{$cl1}{$index}-$seedmean{$cl2}{$index})**2;
							}
							my $dist=$sum**0.5;
							$dist{$cl2}=$dist;
						}
						@tmpseed=sort {$dist{$b} <=> $dist{$a}} @tmpseed;
						$seed=$tmpseed[0];
					}
					if($seed ne ""){
						$merge{$seed}=1;
						$clseq{$cl1}.=$clseq{$seed};
						my @tmpscaffold=();
						if($assign{$cl1}){
							push @tmpscaffold,split /;/,$assign{$cl1};
						}
						if($assign{$seed}){
							push @tmpscaffold,split /;/,$assign{$seed};
						}
						$assign{$cl1}=join(";",@tmpscaffold);
						&CovInterval_4merge($cl1,\@tmpscaffold);
						&Prepare4Zvalue($clseq{$cl1},$cl1);
						if($scg2cl{$seed}){
							foreach my $scg(keys %{$scg2cl{$seed}}) {
								if($scg2cl{$cl1}{$scg}){
									$scg2cl{$cl1}{$scg}+=$scg2cl{$seed}{$scg};
								}
								else{
									$scg2cl{$cl1}{$scg}=$scg2cl{$seed}{$scg};
								}
							}
						}
						if($scg2cl{$seed}){
							delete $scg2cl{$seed};
						}
						delete $assign{$seed};
						delete $clseq{$seed};
					}
					else{
						last;
					}
				}
				else{
					last;
				}
			}
			push @tmptotalseed,$cl1;	
		}
		@totalseed=@tmptotalseed;
		if($seednum == scalar @totalseed){
			last;
		}
	}
}

###################################################################################  ReMerge classification ####################################################
sub Remergeseed{
	my %scg2cl=();
	%clseq=();
	my %smallseed=();
	my %totallen=();
	foreach my $cl(@totalseed){
		my %hash=();
		my $SCGseq="";
		my $totalseq="";
		my @SCGscaffold=();
		if($assign{$cl}){
			my @clscaf=split /;/,$assign{$cl};
			foreach my $scaf(@clscaf){
				$totalseq.=$seq{$scaf};
				if($scg{$scaf}){
					$SCGseq.=$seq{$scaf};
					$hash{$scaf}=1;
					push @SCGscaffold,$scaf;
					foreach my $scg (keys %{$scg{$scaf}}) {
						$scg2cl{$cl}{$scg}++;
					}
				}
			}
			if(scalar @SCGscaffold >=10){
				&SpeciesCovInterval($cl,\@SCGscaffold);
				if(length $SCGseq >=10000){
					$clseq{$cl}=$SCGseq;
				}
				else{
					$clseq{$cl}=$SCGseq;
					my @scaffold=split /;/,$cl;
					foreach my $scaf (@scaffold) {
						if(!$hash{$scaf}){
							$clseq{$cl}.=$seq{$scaf};
						}
					}
				}
			}
			else{
				&CovInterval_4merge($cl,\@clscaf);
				$clseq{$cl}=$totalseq;
			}
			my $totallen=length $totalseq;
			$totallen{$cl}=$totallen;
			if($totallen <$BinLen){
				$smallseed{$cl}=1;
			}
		}
	}
	&PrepareForCorrelation;
	while(){
		my @clid=sort {$seedmean{$b}{1} <=> $seedmean{$a}{1}} @totalseed;
		my $seednum=scalar @totalseed;
		my @tmptotalseed=();
		my %merge=();
		for(my $i=0;$i<=$#clid;$i++){
			my $cl1=$clid[$i];
			next if($merge{$cl1});
			my $refxx=$refxx{$cl1};
			my $SSxx=$SSxx{$cl1};
			my $len1=$len{$cl1};
			while(){
				my %pvalue=();
				my @source=();
				my %FLAG=();
				for(my $j=$i+1;$j<=$#clid;$j++){
					my $cl2=$clid[$j];
					next if($merge{$cl2});
					next if($cl1 eq $cl2);
					my %allscg=();
					if($scg2cl{$cl1}){
						foreach my $scg(keys %{$scg2cl{$cl1}}) {
							$allscg{$scg}+=$scg2cl{$cl1}{$scg};
						}
					}
					if($scg2cl{$cl2}){
						foreach my $scg(keys %{$scg2cl{$cl2}}) {
							$allscg{$scg}+=$scg2cl{$cl2}{$scg};
						}
					}
					my $dup=0;
					foreach my $scg (keys %allscg) {
						#if($scg2cl{$cl1} && $scg2cl{$cl1}{$scg} && $scg2cl{$cl2} && $scg2cl{$cl2}{$scg}){
						if($allscg{$scg} >=2){
							$dup++;
						}
					}
					if($dup <=$scgNum){
						my $refyy=$refxx{$cl2};
						my $SSyy=$SSxx{$cl2};
						my $len2=$len{$cl2};
						my $SSxy=&SS($refxx,$refyy,__LINE__);
						my $cor=&correlation($SSxx,$SSyy,$SSxy);
						my $cutoff=$cutoff{"$len1\t$len2"};
						if($cor >=$cutoff){
							push @source,$cl2;
							my $same_mean=$smean{"$len1\t$len2"};
							my $same_sd=$ssd{"$len1\t$len2"};
							my $val=($cor-$same_mean)/$same_sd;
							my $pvalue=uprob($val);
							$pvalue=1-$pvalue;
							$pvalue{$cl2}=$pvalue;
						}
						else{
							if($smallseed{$cl1} || $smallseed{$cl2}){
								push @source,$cl2;
								my $same_mean=$smean{"$len1\t$len2"};
								my $same_sd=$ssd{"$len1\t$len2"};
								my $val=($cor-$same_mean)/$same_sd;
								my $pvalue=uprob($val);
								$pvalue=1-$pvalue;
								$pvalue{$cl2}=$pvalue;
								if($smallseed{$cl1}){
									$FLAG{$cl1}=1;
								}
								if($smallseed{$cl2}){
									$FLAG{$cl2}=1;
								}
							}
							else{
								next;
							}
						}
					}
					else{
						next;
					}
				}
				if(scalar @source >0){
					my @tmpsource=();
					foreach my $cl2 (@source) {
						my $Num=0;
						my $cha=0;
						foreach my $index (@CovIndex) {
							if($seedmean{$cl1}{$index} <5){
								$cha++;
							}
							else{
								if($seedmean{$cl2}{$index} <=$seedmax{$cl1}{$index} && $seedmean{$cl2}{$index} >=$seedmin{$cl1}{$index}){
									$Num++;
								}
							}
						}
						if($Num+$cha == $IndexNum){
							push @tmpsource,$cl2;
						}
					}
					my $seed="";
					if(scalar @tmpsource >0){
						@tmpsource=sort {$pvalue{$b} <=> $pvalue{$a}} @tmpsource;
						$seed=$tmpsource[0];
					}
					if($seed ne ""){
						$merge{$seed}=1;
						my @tmpscaffold=();
						my $seq="";
						if($assign{$cl1}){
							if($FLAG{$cl1}){
								my @totalscaffold=split /;/,$assign{$cl1};
								foreach my $scaf (@totalscaffold) {
									if($scg{$scaf}){
										$seq.=$seq{$scaf};
									}
									push @tmpscaffold,$scaf;
								}
							}
							else{
								push @tmpscaffold,split /;/,$assign{$cl1};
								$seq.=$clseq{$cl1};
							}
						}
						if($assign{$seed}){
							if($FLAG{$seed}){
								my @totalscaffold=split /;/,$assign{$seed};
								foreach my $scaf (@totalscaffold) {
									if($scg{$scaf}){
										$seq.=$seq{$scaf};
										foreach my $scg (keys %{$scg{$scaf}}) {
											if($scg2cl{$cl1}{$scg}){
												$scg2cl{$cl1}{$scg}+=$scg2cl{$seed}{$scg};
											}
											else{
												$scg2cl{$cl1}{$scg}=$scg2cl{$seed}{$scg};
											}
										}
									}
									push @tmpscaffold,$scaf;
								}
							}
							else{
								push @tmpscaffold,split /;/,$assign{$seed};
								$seq.=$clseq{$seed};
								my @totalscaffold=split /;/,$assign{$seed};
								foreach my $scaf (@totalscaffold) {
									if($scg{$scaf}){
										foreach my $scg (keys %{$scg{$scaf}}) {
											if($scg2cl{$cl1}{$scg}){
												$scg2cl{$cl1}{$scg}+=$scg2cl{$seed}{$scg};
											}
											else{
												$scg2cl{$cl1}{$scg}=$scg2cl{$seed}{$scg};
											}
										}
									}
								}
							}
						}
						$clseq{$cl1}=$seq;
						$assign{$cl1}=join(";",@tmpscaffold);
						$totallen{$cl1}+=$totallen{$seed};
						if($totallen{$cl1} >$BinLen){
							if($smallseed{$cl1}){
								#delete $smallseed{$cl1};
							}
						}
						&Prepare4Zvalue($clseq{$cl1},$cl1);
						if($scg2cl{$seed}){
							delete $scg2cl{$seed};
						}
						delete $assign{$seed};
						delete $clseq{$seed};
						if($smallseed{$seed}){
							delete $smallseed{$seed};
						}
						delete $totallen{$seed};
					}
					else{
						last;
					}
				}
				else{
					last;
				}
			}
			push @tmptotalseed,$cl1;	
		}
		@totalseed=@tmptotalseed;
		if($seednum == scalar @totalseed){
			last;
		}
	}
}

################################################################################### classify large scaffolds ####################################################
sub Classify_LargeScaffold{
	my ($ref1,$ref2)=@_;
	my @seed=@{$ref1};
	my @largescaffold=@{$ref2};
			
	##### classify scaffold with length >=10kb#########################
	my ($refxx,$SSxx,$len1,$refyy,$SSyy,$len2);
	foreach my $scaf(@largescaffold){
		my @source=();
		foreach my $cl(@seed){
			my $Num=0;
			my $cha=0;
			foreach my $index(@CovIndex){
				if($seedmean{$cl}{$index} <5){
					$cha++;
				}
				elsif($mean{$scaf}{$index}>=$seedmin{$cl}{$index} && $mean{$scaf}{$index} <=$seedmax{$cl}{$index}){
					$Num++;
				}
			}
			if($Num+$cha == $IndexNum){
				push @source,$cl;
			}
		}
		if(scalar @source ==0){
			next;
		}
		elsif(scalar @source ==1){
			if($assign{$source[0]}){
				$assign{$source[0]}.=";$scaf";
			}
			else{
				$assign{$source[0]}=$scaf;
			}
		}
		else{
			$refxx=$refxx{$scaf};
			$SSxx=$SSxx{$scaf};
			$len1=$len{$scaf};
			my $maxpvalue="NA";
			my $clid="";
			my $maxcor=0;
			foreach my $cl(@source){
				$refyy=$refxx{$cl};
				$SSyy=$SSxx{$cl};
				$len2=$len{$cl};
				my $SSxy=&SS($refxx,$refyy,__LINE__);
				my $cor=&correlation($SSxx,$SSyy,$SSxy);
				my $same_mean=$smean{"$len1\t$len2"};
				my $same_sd=$ssd{"$len1\t$len2"};
				my $val=($cor-$same_mean)/$same_sd;
				my $pvalue=uprob($val);
				$pvalue=1-$pvalue;
				if($maxpvalue eq "NA"){
					$maxpvalue=$pvalue;
					$maxcor=$cor;
					$clid=$cl;
				}
				else{
					if($pvalue >$maxpvalue){
						$maxpvalue=$pvalue;
						$clid=$cl;
						$maxcor=$cor;
					}
					elsif($pvalue == $maxpvalue){
						if($cor>$maxcor){
							$clid=$cl;
							$maxcor=$cor;
						}
					}
				}
			}
			if($assign{$clid}){
				$assign{$clid}.=";$scaf";
			}
			else{
				$assign{$clid}=$scaf;
			}
		}
	}
}

####################################################################### trim abnormal depth ################################
sub depth_process{
	my ($ref)=@_;
	my @data=@{$ref};
	@data=sort {$a <=>$b} @data;
	my @tmpdata=();
	my ($IQR,$Q1,$Q3);
	if((@data+1) % 4 ==0){
		my $index1=(@data+1)/4;
		my $index2=3*(@data+1)/4;
		$Q1=$data[$index1-1];
		$Q3=$data[$index2-1];
		$IQR=$Q3-$Q1;
	}    
	else{
		my $val=(@data+1)/4;
		my $int=int($val);
		my $decimal=$val-$int;
		$Q1=$decimal* $data[$int-1]+(1-$decimal)*$data[$int];
		$val=3*(@data+1)/4;
		$int=int($val);
		$decimal=$val-$int;
		$Q3=$decimal* $data[$int-1]+(1-$decimal)*$data[$int];
		$IQR=$Q3-$Q1;
	}    
	my $tmin=$Q1-1.5*$IQR;
	my $tmax=$Q3+1.5*$IQR;
	return ($tmin,$tmax);
}

####################################################################### SD ################################################################################
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

####################################################################### Subroutine for computing mean################################################################
sub mean{
	my ($ref)=@_;
	my @a=@{$ref};
	my $sum=0;
	foreach(@a){
		$sum+=$_;
	}
	my $mean=$sum/@a;
	return $mean;
}

####################################################################### Subroutine for computing median################################################################
sub median{
	my ($ref)=@_;
	my @val=@{$ref};
	@val=sort {$a <=> $b} @val;
	my $median;
	if(@val % 2 ==0){
		$median=($val[@val/2]+$val[(@val/2)-1])/2;
	}
	else{
		$median=$val[$#val/2];
	}
	return $median;
}

####################################################################### Computing zvalue ###################################################################
sub zvalue{
	my ($seq)=@_;
	my (%kmer,%k_1mer,%k_2mer)=((),(),());
	$seq=~s/[^ATGC]//gi;
	my $rev=reverse $seq;
	$rev=~tr/ATGC/TACG/;
	$seq.=$rev;
	my $len=length $seq;
	foreach(@oligo_kmer){
		$kmer{$_}=1;
	}
	foreach(@oligo_k_1mer){
		$k_1mer{$_}=1;
	}
	foreach(@oligo_k_2mer){
		$k_2mer{$_}=1;
	}
	for(my $i=0;$i<=$len-4;$i++){
		my $sub=substr($seq,$i,4);
		$kmer{$sub}++;
	}
	for(my $i=0;$i<=$len-3;$i++){
		my $sub=substr($seq,$i,3);
		$k_1mer{$sub}++;
	}
	for (my $i=0;$i<=$len-2;$i++){
		my $sub=substr($seq,$i,2);
		$k_2mer{$sub}++;
	}
	my @zvalue=();
	foreach(sort keys %kmer){
		my $N_koligo=$kmer{$_};
		my $N_former=$k_1mer{substr($_,0,3)};
		my $N_latter=$k_1mer{substr($_,1,3)};
		my $N_midder=$k_2mer{substr($_,1,2)};
		my $denominator=($N_midder-$N_former)*($N_midder-$N_latter)*$N_former*$N_latter;
		my $sqrt=$denominator**0.5;
		if($sqrt !=0 and $N_midder !=0){
			my $zvalue=$N_midder**0.5*($N_koligo*$N_midder-$N_former*$N_latter)/$sqrt;
			push @zvalue,$zvalue;
		}
		else{
			push @zvalue,0;
		}
	}
	my $ref=\@zvalue;
	return $ref;
}

####################################################################### Computing kmer freq for small scaffold ###################################################################
sub kmerFreq{
	my ($scaf)=@_;
	my $seq=$seq{$scaf};
	$seq=~s/[^ATGC]//gi;
	my $len=length $seq;
	for(my $i=0;$i<=$len-4;$i++){
		my $sub=substr($seq,$i,4);
		$freq{$scaf}{$sub}++;
	}
}

##################################################################### Subroutine for computing covariance or variance#################################################
sub SS{
	my ($ref1,$ref2,$line)=@_;
	if(!$ref1){
		print STDERR "$line\n";
	}
	my @a1=@{$ref1};
	if(!$ref2){
		print STDERR "$line\n";
	}
	my @a2=@{$ref2};
	my $sum=0;
	for (my $k=0;$k<=$#a1;$k++){
		$sum+=$a1[$k]*$a2[$k];
	}
	return $sum;
}

################################################################# Subroutine for computing Pearson correlation coefficient ############################################
sub correlation {
	my ($ssxx,$ssyy,$ssxy)=@_;
	my $correl=0;
	if($ssxy !=0 && $ssxx !=0 && $ssyy!=0){
		my $sign=$ssxy/abs($ssxy);
		$correl=$sign*sqrt($ssxy*$ssxy/($ssxx*$ssyy));
	}
	else{
		$correl="0";
	}
	return $correl;
}

################################################################# Subroutine for Prepare for correlation #################################################################
sub PrepareForCorrelation{
	foreach my $cl(keys %clseq){
		my $seq=$clseq{$cl};
		my $len=length $seq;
		if($len <10000){
			print LOG "$cl: $len\n";
			print STDERR "$cl: $len\n";
		}
		$len = int($len / 10000) * 10;
		if($len > 210){
			$len = 210;
		}
		$len{$cl} = $len;
		my $zvalueref=&zvalue($seq);
		my $mean=&mean($zvalueref);
		my @c=();
		foreach(@{$zvalueref}){
			my $tmp=$_-$mean;
			push @c,$tmp;
		}
		my $refyy=\@c;
		my $SSyy=&SS($refyy,$refyy,__LINE__);
		$refxx{$cl}=$refyy;
		$SSxx{$cl}=$SSyy;
		$len{$cl}=$len;
	}
}

################################################################################### Subroutine prior processing for Zvalue ####################################################
sub Prepare4Zvalue{
	my ($seq,$cl)=@_;
	my $len=length $seq;
	$len=int($len/10000)*10;
	if($len >210){
		$len=210;
	}
	my $zvalueref=&zvalue($seq);
	my $mean=&mean($zvalueref);
	my @d=();
	foreach(@{$zvalueref}){
		my $tmp=$_-$mean;
		push @d,$tmp;
	}
	$refxx=\@d;
	$SSxx=&SS($refxx,$refxx,__LINE__);
	$refxx{$cl}=$refxx;
	$SSxx{$cl}=$SSxx;
	$len{$cl}=$len;
}

################################################################################### Subroutine prior processing for Beyesin ####################################################
sub Prepare4Bayesin{
	my ($seq,$cl,$line)=@_;
	if($seq eq ""){
		print STDERR "$line\n";
	}
	my $length=length $seq;
	my %hash_kmer=();
	my %hash_sum=();
	for(my $i=0;$i<=$length-4;$i++){
		my $sub=substr($seq,$i,4);
		$hash_kmer{$cl}{$sub}++;
	}
	$hash_sum{$cl}=$length-4+1;
	foreach(@oligo_kmer){
		my $val;
		if($hash_kmer{$cl}{$_}){
			$val=$hash_kmer{$cl}{$_}/$hash_sum{$cl};
		}
		else{
			if($hash_sum{$cl}){
				$val=1/$hash_sum{$cl};
			}
		}
		$prob{$cl}{$_}=log($val)/log(10);
	}
}

########################################################################## Subroutine for classify duplicated SCG scaffolds ####################################################
sub duplicated_SCGscaffold_Allclassify{
	my ($ref1,$ref2,$ref3)=@_;
	my @seed=@{$ref1};
	my @duplicatedSCGscaffold_more10kb=@{$ref2};
	my @duplicatedSCGscaffold_less10kb=@{$ref3};
	my ($refxx,$SSxx,$len1,$refyy,$SSyy,$len2);
	my %possiblesource=();
	my @remained_more10kb=();
	my @remained_less10kb=();
	foreach my $scaf(@duplicatedSCGscaffold_more10kb){
		my @source=();
		foreach my $cl(@seed){
			my $Num=0;
			my $cha=0;
			foreach my $index (@CovIndex) {
				if($seedmean{$cl}{$index} <5){
					$cha++;
				}
				elsif($mean{$scaf}{$index}>=$seedmin{$cl}{$index} && $mean{$scaf}{$index} <=$seedmax{$cl}{$index}){
					$Num++;
				}
			}
			if($Num+$cha == $IndexNum){
				push @source,$cl;
			}
		}
		if(scalar @source ==0){
			next;
		}
		elsif(scalar @source ==1){
			if($SCGassign{$source[0]}){
				$SCGassign{$source[0]}.=";$scaf";
			}
			else{
				$SCGassign{$source[0]}=$scaf;
			}
			if($scg{$scaf}){
				foreach my $scg (keys %{$scg{$scaf}}) {
					$Flagscg2cl{$scg}{$source[0]}=1;
				}
			}
			if($duplicate{$scaf}){
				my @excludedScaffold=split /;/,$duplicate{$scaf};
				foreach(@excludedScaffold){
					$exclude{$_}{$source[0]}=1;
				}
			}
		}
		else{
			$possiblesource{$scaf}=\@source;
			push @remained_more10kb,$scaf;
		}
	}

	my %classified=();
	foreach my $scaf(@duplicatedSCGscaffold_less10kb){
		my @source=();
		foreach my $cl(@seed){
			my $Num=0;
			foreach my $index (@CovIndex) {
				if($mean{$scaf}{$index}>=$seedmin{$cl}{$index} && $mean{$scaf}{$index} <=$seedmax{$cl}{$index}){
					$Num++;
				}
			}
			if($Num == $IndexNum){
				push @source,$cl;
			}
		}
		if(scalar @source ==0){
			next;
		}
		elsif(scalar @source ==1){
			if($SCGassign{$source[0]}){
				$SCGassign{$source[0]}.=";$scaf";
			}
			else{
				$SCGassign{$source[0]}=$scaf;
			}
			if($scg{$scaf}){
				foreach my $scg (keys %{$scg{$scaf}}) {
					$Flagscg2cl{$scg}{$source[0]}=1;
				}
			}
			if($duplicate{$scaf}){
				my @excludedScaffold=split /;/,$duplicate{$scaf};
				foreach(@excludedScaffold){
					$exclude{$_}{$source[0]}=1;
				}
			}
		}
		else{
			$possiblesource{$scaf}=\@source;
			push @remained_less10kb,$scaf;
		}
	}
	
	while(scalar @remained_more10kb >0){
		my %pvalue=();
		foreach my $scaf (@remained_more10kb) {
			my @tmpsource=@{$possiblesource{$scaf}};
			my @source=();
			foreach my $cl (@tmpsource) {
				if($exclude{$scaf}{$cl}){
					next;
				}
				else{
					push @source,$cl;
				}
			}
			$refxx=$refxx{$scaf};
			$SSxx=$SSxx{$scaf};
			$len1=$len{$scaf};
			if(scalar @source == 0){
				@source=@tmpsource;
				my $maxpvalue="NA";
				my $clid="";
				my $maxcor=0;
				foreach my $cl(@source){
					$refyy=$refxx{$cl};
					$SSyy=$SSxx{$cl};
					$len2=$len{$cl};
					my $SSxy=&SS($refxx,$refyy,__LINE__);
					my $cor=&correlation($SSxx,$SSyy,$SSxy);
					my $same_mean=$smean{"$len1\t$len2"};
					my $same_sd=$ssd{"$len1\t$len2"};
					my $val=($cor-$same_mean)/$same_sd;
					my $pvalue=uprob($val);
					$pvalue=1-$pvalue;
					if($maxpvalue eq "NA"){
						$maxpvalue=$pvalue;
						$maxcor=$cor;
						$clid=$cl;
					}
					else{
						if($pvalue >$maxpvalue){
							$maxpvalue=$pvalue;
							$clid=$cl;
							$maxcor=$cor;
						}
						elsif($pvalue == $maxpvalue){
							if($cor>$maxcor){
								$clid=$cl;
								$maxcor=$cor;
							}
						}
					}
				}
				if($SCGassign{$clid}){
					$SCGassign{$clid}.=";$scaf";
				}
				else{
					$SCGassign{$clid}=$scaf;
				}
				if($scg{$scaf}){
					foreach my $scg (keys %{$scg{$scaf}}) {
						$Flagscg2cl{$scg}{$clid}=1;
					}
				}
				if($duplicate{$scaf}){
					my @excludedScaffold=split /;/,$duplicate{$scaf};
					foreach(@excludedScaffold){
						$exclude{$_}{$clid}=1;
					}
				}
			}
			elsif(scalar @source == 1){
				my $clid=$source[0];
				if($SCGassign{$clid}){
					$SCGassign{$clid}.=";$scaf";
				}
				else{
					$SCGassign{$clid}=$scaf;
				}
				if($scg{$scaf}){
					foreach my $scg (keys %{$scg{$scaf}}) {
						$Flagscg2cl{$scg}{$clid}=1;
					}
				}
				if($duplicate{$scaf}){
					my @excludedScaffold=split /;/,$duplicate{$scaf};
					foreach(@excludedScaffold){
						$exclude{$_}{$clid}=1;
					}
				}
			}
			else{
				my $maxpvalue="NA";
				my $clid="";
				my $maxcor=0;
				foreach my $cl(@source){
					$refyy=$refxx{$cl};
					$SSyy=$SSxx{$cl};
					$len2=$len{$cl};
					my $SSxy=&SS($refxx,$refyy,__LINE__);
					my $cor=&correlation($SSxx,$SSyy,$SSxy);
					my $same_mean=$smean{"$len1\t$len2"};
					my $same_sd=$ssd{"$len1\t$len2"};
					my $val=($cor-$same_mean)/$same_sd;
					my $pvalue=uprob($val);
					$pvalue=1-$pvalue;
					if($maxpvalue eq "NA"){
						$maxpvalue=$pvalue;
						$maxcor=$cor;
						$clid=$cl;
					}
					else{
						if($pvalue >$maxpvalue){
							$maxpvalue=$pvalue;
							$clid=$cl;
							$maxcor=$cor;
						}
						elsif($pvalue == $maxpvalue){
							if($cor>$maxcor){
								$clid=$cl;
								$maxcor=$cor;
							}
						}
					}
				}
				$pvalue{$clid}{$scaf}=$maxpvalue;
			}
		}
		
		my @tmpremained_more10kb=();
		foreach my $cl (keys %pvalue) {
			my @scaffold=sort {$pvalue{$cl}{$b} <=> $pvalue{$cl}{$a}} keys %{$pvalue{$cl}};
			foreach my $scaf (@scaffold) {
				if($exclude{$scaf}{$cl}){
					push @tmpremained_more10kb,$scaf;
				}
				else{
					if($SCGassign{$cl}){
						$SCGassign{$cl}.=";$scaf";
					}
					else{
						$SCGassign{$cl}=$scaf;
					}
					if($scg{$scaf}){
						foreach my $scg (keys %{$scg{$scaf}}) {
							$Flagscg2cl{$scg}{$cl}=1;
						}
					}
					if($duplicate{$scaf}){
						my @excludedScaffold=split /;/,$duplicate{$scaf};
						foreach(@excludedScaffold){
							$exclude{$_}{$cl}=1;
						}
					}
				}
			}
		}
		@remained_more10kb =@tmpremained_more10kb;
	}	

	while(@remained_less10kb >0){
		my %pvalue=();
		foreach my $scaf (@remained_less10kb) {
			my @tmpsource=@{$possiblesource{$scaf}};
			my @source=();
			foreach my $cl (@tmpsource) {
				if($exclude{$scaf}{$cl}){
					next;
				}
				else{
					push @source,$cl;
				}
			}
			if(scalar @source == 0){
				@source=@tmpsource;
				my ($clstr,$foldval)=&bayesin($scaf,\@source,\%prob,__LINE__);
				my @clstr=split /\#/,$clstr;
				if($SCGassign{$clstr[0]}){
					$SCGassign{$clstr[0]}.=";$scaf";
				}
				else{
					$SCGassign{$clstr[0]}=$scaf;
				}
				if($scg{$scaf}){
					foreach my $scg (keys %{$scg{$scaf}}) {
						$Flagscg2cl{$scg}{$clstr[0]}=1;
					}
				}
				if($duplicate{$scaf}){
					my @excludedScaffold=split /;/,$duplicate{$scaf};
					foreach(@excludedScaffold){
						$exclude{$_}{$clstr[0]}=1;
					}
				}
			}
			elsif(scalar @source == 1){
				if($SCGassign{$source[0]}){
					$SCGassign{$source[0]}.=";$scaf";
				}
				else{
					$SCGassign{$source[0]}=$scaf;
				}
				if($scg{$scaf}){
					foreach my $scg (keys %{$scg{$scaf}}) {
						$Flagscg2cl{$scg}{$source[0]}=1;
					}
				}
				if($duplicate{$scaf}){
					my @excludedScaffold=split /;/,$duplicate{$scaf};
					foreach(@excludedScaffold){
						$exclude{$_}{$source[0]}=1;
					}
				}
			}
			else{
				my ($clstr,$foldval)=&bayesin($scaf,\@source,\%prob,__LINE__);
				my @clstr=split /\#/,$clstr;
				$pvalue{$clstr[0]}{$scaf}=$foldval;
			}
		}
		
		my @tmpremained_less10kb=();
		foreach my $cl (keys %pvalue) {
			my @scaffold=sort {$pvalue{$cl}{$b} <=> $pvalue{$cl}{$a}} keys %{$pvalue{$cl}};
			foreach my $scaf (@scaffold) {
				if($exclude{$scaf}{$cl}){
					push @tmpremained_less10kb,$scaf;
				}
				else{
					if($SCGassign{$cl}){
						$SCGassign{$cl}.=";$scaf";
					}
					else{
						$SCGassign{$cl}=$scaf;
					}
					if($scg{$scaf}){
						foreach my $scg (keys %{$scg{$scaf}}) {
							$Flagscg2cl{$scg}{$cl}=1;
						}
					}
					if($duplicate{$scaf}){
						my @excludedScaffold=split /;/,$duplicate{$scaf};
						foreach(@excludedScaffold){
							$exclude{$_}{$cl}=1;
						}
					}
				}
			}
		}
		@remained_less10kb = @tmpremained_less10kb;
	}
}

################################################################################### Subroutine for classify SCG scaffolds ####################################################
sub SCGscaffold_AllClassify{
	my ($ref1,$ref2,$ref3)=@_;
	my @seed=@{$ref1};
	my @SCGscaffold_more10kb_4class=@{$ref2};
	my @SCGscaffold_less10kb_4class=@{$ref3};
	%Flagscg2cl=();
	%exclude=();
	%SCGassign=();
	my @duplicatedSCGscaffold_less10kb=();
	my @duplicatedSCGscaffold_more10kb=();
	my @OtherSCGscaffold_less10kb=();
	my @OtherSCGscaffold_more10kb=();
	foreach my $scaf (@SCGscaffold_more10kb_4class) {
		if($duplicate{$scaf}){
			push @duplicatedSCGscaffold_more10kb,$scaf;
		}
		else{
			push @OtherSCGscaffold_more10kb,$scaf;
		}
	}
	foreach my $scaf (@SCGscaffold_less10kb_4class) {
		if($duplicate{$scaf}){
			push @duplicatedSCGscaffold_less10kb,$scaf;
		}
		else{
			push @OtherSCGscaffold_less10kb,$scaf;
		}
	}

	&duplicated_SCGscaffold_Allclassify($ref1,\@duplicatedSCGscaffold_more10kb,\@duplicatedSCGscaffold_less10kb);
	##### classify scaffold with length >=10kb#########################
	my ($refxx,$SSxx,$len1,$refyy,$SSyy,$len2);
	foreach my $scaf(@OtherSCGscaffold_more10kb){
		my @source=();
		foreach my $cl(@seed){
			my $Num=0;
			my $cha=0;
			foreach my $index(@CovIndex){
				if($seedmean{$cl}{$index} <5){
					$cha++;
				}
				elsif($mean{$scaf}{$index}>=$seedmin{$cl}{$index} && $mean{$scaf}{$index} <=$seedmax{$cl}{$index}){
					$Num++;
				}
			}
			if($Num+$cha == $IndexNum){
				push @source,$cl;
			}
		}
		if(scalar @source ==0){
			next;
		}
		elsif(scalar @source ==1){
			if($SCGassign{$source[0]}){
				$SCGassign{$source[0]}.=";$scaf";
			}
			else{
				$SCGassign{$source[0]}=$scaf;
			}
		}
		else{
			my @tmpsource=();
			foreach my $cl(@source){
				if($exclude{$scaf}{$cl}){
					next;
				}
				else{
					my $n=0;
					if($scg{$scaf}){
						foreach(keys %{$scg{$scaf}}){
							if($Flagscg2cl{$_}{$cl}){
								$n++;
							}
						}
					}
					if($n <2){
						push @tmpsource,$cl;
					}
				}
			}
			if(scalar @tmpsource ==0){
				@tmpsource=@source;
			}
			if(scalar @tmpsource ==1){
				@source=@tmpsource;
				if($SCGassign{$source[0]}){
					$SCGassign{$source[0]}.=";$scaf";
				}
				else{
					$SCGassign{$source[0]}=$scaf;
				}
				if($scg{$scaf}){
					foreach my $scg (keys %{$scg{$scaf}}) {
						$Flagscg2cl{$scg}{$source[0]}=1;
					}
				}
			}
			else{
				$refxx=$refxx{$scaf};
				$SSxx=$SSxx{$scaf};
				$len1=$len{$scaf};
				@source=@tmpsource;
				my $maxpvalue="NA";
				my $clid="";
				my $maxcor=0;
				foreach my $cl(@source){
					$refyy=$refxx{$cl};
					$SSyy=$SSxx{$cl};
					$len2=$len{$cl};
					my $SSxy=&SS($refxx,$refyy,__LINE__);
					my $cor=&correlation($SSxx,$SSyy,$SSxy);
					my $same_mean=$smean{"$len1\t$len2"};
					my $same_sd=$ssd{"$len1\t$len2"};
					my $val=($cor-$same_mean)/$same_sd;
					my $pvalue=uprob($val);
					$pvalue=1-$pvalue;
					if($maxpvalue eq "NA"){
						$maxpvalue=$pvalue;
						$maxcor=$cor;
						$clid=$cl;
					}
					else{
						if($pvalue >$maxpvalue){
							$maxpvalue=$pvalue;
							$clid=$cl;
							$maxcor=$cor;
						}
						elsif($pvalue == $maxpvalue){
							if($cor>$maxcor){
								$clid=$cl;
								$maxcor=$cor;
							}
						}
					}
				}
				if($SCGassign{$clid}){
					$SCGassign{$clid}.=";$scaf";
				}
				else{
					$SCGassign{$clid}=$scaf;
				}
			}
		}
	}
	##### classify scaffold with length <10kb#########################
	&MultipleThreads_SCGscaffold_AllBayesin(\@seed,\@OtherSCGscaffold_less10kb);
}

################################################################################### Subroutine for multiple-thread confident Bayesin####################################################
sub MultipleThreads_SCGscaffold_AllBayesin{
	my ($ref1,$ref2)=@_;
	my @seed=@{$ref1};
	my @OtherSCGscaffold_less10kb=@{$ref2};
	my %unclass=();
	foreach my $scaf (@OtherSCGscaffold_less10kb) {
		my @source=();
		foreach my $cl (@seed) {
			my $num=0;
			my $cha=0;
			foreach my $index (@CovIndex) {
				if($seedmean{$cl}{$index} <5){
					$cha++;
				}
				elsif($mean{$scaf}{$index} >=$seedmin{$cl}{$index} && $mean{$scaf}{$index} <=$seedmax{$cl}{$index}){
					$num++;
				}
			}
			if($num+$cha == $IndexNum){
				push @source,$cl;
			}
		}
		if(@source == 0){
			next;
		}
		elsif(scalar @source ==1){
			if($SCGassign{$source[0]}){
				$SCGassign{$source[0]}.=";$scaf";
			}
			else{
				$SCGassign{$source[0]}="$scaf";
			}
			if($scg{$scaf}){
				foreach my $scg (keys %{$scg{$scaf}}) {
					$Flagscg2cl{$scg}{$source[0]}=1;
				}
			}
			if($duplicate{$scaf}){
				my @excludedScaffold=split /;/,$duplicate{$scaf};
				foreach(@excludedScaffold){
					$exclude{$_}{$source[0]}=1;
				}
			}
		}
		else{
			my @tmpsource=();
			foreach my $cl(@source){
				if($exclude{$scaf}{$cl}){
					next;
				}
				else{
					my $n=0;
					if($scg{$scaf}){
						foreach (keys %{$scg{$scaf}}) {
							if($Flagscg2cl{$_}{$cl}){
								$n++;
							}
						}
					}
					if($scg{$scaf} && $n != scalar keys %{$scg{$scaf}} && $n <2){
						push @tmpsource,$cl;
					}
				}
			}
			if(scalar @tmpsource ==0){
				my ($clstr,$foldval)=&bayesin($scaf,\@source,\%prob,__LINE__);
				my @clstr=split /\#/,$clstr;
				if($SCGassign{$clstr[0]}){
					$SCGassign{$clstr[0]}.=";$scaf";
				}
				else{
					$SCGassign{$clstr[0]}="$scaf";
				}
				if($scg{$scaf}){
					foreach my $scg (keys %{$scg{$scaf}}) {
						$Flagscg2cl{$scg}{$clstr[0]}=1;
					}
				}
				if($duplicate{$scaf}){
					my @excludedScaffold=split /;/,$duplicate{$scaf};
					foreach(@excludedScaffold){
						$exclude{$_}{$clstr[0]}=1;
					}
				}
			}
			elsif(scalar @tmpsource ==1){
				@source=@tmpsource;
				if($SCGassign{$source[0]}){
					$SCGassign{$source[0]}.=";$scaf";
				}
				else{
					$SCGassign{$source[0]}="$scaf";
				}
				if($scg{$scaf}){
					foreach my $scg (keys %{$scg{$scaf}}) {
						$Flagscg2cl{$scg}{$source[0]}=1;
					}
				}
				if($duplicate{$scaf}){
					my @excludedScaffold=split /;/,$duplicate{$scaf};
					foreach(@excludedScaffold){
						$exclude{$_}{$source[0]}=1;
					}
				}
			}
			else{
				@source=@tmpsource;
				my ($clstr,$foldval)=&bayesin($scaf,\@source,\%prob,__LINE__);
				if($clstr !~/\#/){
					if($SCGassign{$clstr}){
						$SCGassign{$clstr}.=";$scaf";
					}
					else{
						$SCGassign{$clstr}="$scaf";
					}
					if($scg{$scaf}){
						foreach my $scg (keys %{$scg{$scaf}}) {
							$Flagscg2cl{$scg}{$clstr}=1;
						}
					}
					if($duplicate{$scaf}){
						my @excludedScaffold=split /;/,$duplicate{$scaf};
						foreach(@excludedScaffold){
							$exclude{$_}{$clstr}=1;
						}
					}
				}
				else{
					$unclass{$scaf}=\@source;
				}
			}
		}
	}

	####################classify remaining scaffolds ######################
	while(scalar keys %unclass >0){
		my %foldFlag=();
		foreach my $scaf (keys %unclass) {
			my @source=@{$unclass{$scaf}};
			my @tmpsource=();
			foreach my $cl(@source){
				if($exclude{$scaf}{$cl}){
					next;
				}
				else{
					my $n=0;
					if($scg{$scaf}){
						foreach(keys %{$scg{$scaf}}){
							if($Flagscg2cl{$_}{$cl}){
								$n++;
							}
						}
					}
					if($scg{$scaf} && $n != scalar keys %{$scg{$scaf}} && $n <2){
						push @tmpsource,$cl;
					}
				}
			}
			if(scalar @tmpsource ==0){
				my ($clstr,$foldval)=&bayesin($scaf,\@source,\%prob,__LINE__);
				my @clstr=split /\#/,$clstr;
				if($SCGassign{$clstr[0]}){
					$SCGassign{$clstr[0]}.=";$scaf";
				}
				else{
					$SCGassign{$clstr[0]}="$scaf";
				}
				if($scg{$scaf}){
					foreach my $scg (keys %{$scg{$scaf}}) {
						$Flagscg2cl{$scg}{$clstr[0]}=1;
					}
				}
				if($duplicate{$scaf}){
					my @excludedScaffold=split /;/,$duplicate{$scaf};
					foreach(@excludedScaffold){
						$exclude{$_}{$clstr[0]}=1;
					}
				}
			}
			elsif(scalar @tmpsource ==1){
				@source=@tmpsource;
				if($SCGassign{$source[0]}){
					$SCGassign{$source[0]}.=";$scaf";
				}
				else{
					$SCGassign{$source[0]}="$scaf";
				}
				if($scg{$scaf}){
					foreach my $scg (keys %{$scg{$scaf}}) {
						$Flagscg2cl{$scg}{$source[0]}=1;
					}
				}
				if($duplicate{$scaf}){
					my @excludedScaffold=split /;/,$duplicate{$scaf};
					foreach(@excludedScaffold){
						$exclude{$_}{$source[0]}=1;
					}
				}
			}
			else{
				my ($clstr,$foldval)=&bayesin($scaf,\@tmpsource,\%prob,__LINE__);
				my @clstr=split /\#/,$clstr;
				if($clstr !~/\#/){
					if($scg{$scaf}){
						foreach my $scg (keys %{$scg{$scaf}}) {
							$Flagscg2cl{$scg}{$clstr[0]}=1;
						}
					}
					if($duplicate{$scaf}){
						my @excludedScaffold=split /;/,$duplicate{$scaf};
						foreach(@excludedScaffold){
							$exclude{$_}{$clstr[0]}=1;
						}
					}
					if($SCGassign{$clstr[0]}){
						$SCGassign{$clstr[0]}.=";$scaf";
					}
					else{
						$SCGassign{$clstr[0]}="$scaf";
					}
				}
				else{
					my $clid=$clstr[0];
					$foldFlag{$clid}{$scaf}=$foldval;
				}
			}
		}

		my %tmpUnclass=();
		foreach my $seed (keys %foldFlag) {
			my @scaffold=sort {$foldFlag{$seed}{$b} <=> $foldFlag{$seed}{$a}} keys %{$foldFlag{$seed}};
			foreach my $scaf (@scaffold) {
				if($exclude{$scaf}{$seed}){
					$tmpUnclass{$scaf}=$unclass{$scaf};
				}
				else{
					my $num=0;
					if($scg{$scaf}){
						foreach my $scg (keys %{$scg{$scaf}}) {
							if($Flagscg2cl{$scg}{$seed}){
								$num++;
							}
						}
					}
					if($scg{$scaf} && $num != scalar keys %{$scg{$scaf}} && $num <2){
						if($scg{$scaf}){
							foreach my $scg (keys %{$scg{$scaf}}) {
								$Flagscg2cl{$scg}{$seed}=1;
							}
						}
						if($duplicate{$scaf}){
							my @excludedScaffold=split /;/,$duplicate{$scaf};
							foreach(@excludedScaffold){
								$exclude{$_}{$seed}=1;
							}
						}
						if($SCGassign{$seed}){
							$SCGassign{$seed}.=";$scaf";
						}
						else{
							$SCGassign{$seed}="$scaf";
						}
					}
					else{
						$tmpUnclass{$scaf}=$unclass{$scaf};
					}
				}
			}
		}
		%unclass=%tmpUnclass;
		my $num=scalar %unclass;
	}
}

########################################################################## Usage ##########################################################################
sub usage{
	print STDERR <<USAGE;
Name
	$0
Usage
	The following arguments must be provided:
	-fafile <s>: fasta file containing scaffolds assembled from metagenomic reads
	-depthfile <s>: the file containing depths for scaffolds
	-tablefile <s>: the file containing SCG type and their Scaffolds
	The following arguments are optional:
	-seedfile <s>: the file containing seeds
	-nb_process <i>: the number of threads to be used, the default is 20
	-fold <f>: the fold for computing depth range, the default is 2.58
	-minlen <i>: the minimal length of scaffold, not including Ns,the default is 0
	-outdir <s>: the output directory
	-LargeLen <i>: the minimal length of large scaffold derived seeds, the default is 200,000
	-SmallLen <i>: the minimal length of small scaffold derived seeds, the default is 100,000
	-BinLen <i>: the minimal length of bins, the default is 500,000
	-scgNum <i>: the minimal number of scgs which can be considered as an bin, the default is 7
	-scgNum4bin <i>: the minimal number of scgs which can be considered as an final bin after merging, the default is 10
	-help: show the help message
USAGE
}
