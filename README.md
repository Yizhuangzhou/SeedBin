# SeedBin v1.0  
A program to automatically classify metagenomic sequences by using progressing composition and coverage (single or multiple samples), also similarity and features of single copy genes (SCGs). 

## Please Cite ##
If you use SeedBin in your publication, please cite:

## Key features ##
•	The state-of-the-art binning approaches combined only composition and coverage. It is the first time for us to combine four binning-useful information including composition, coverage, similarity and SCG features

•	Our method was devised to select seeds from within metagenomic data without any extra effect such as primer walking. Additionally, one algorithm was devised to select small-scaffold seeds, which is useful when no large scaffold is obtained for some species. Validation shows that this algorithm is very effective.

•	SCGs, which were previously used to assess binning performance, were found to carry two binning-useful features in our study and thus integrated into our method for binning. SCG features are indispensable to classify the scaffolds with indiscernible composition and coverage.

•	Most previous approaches used static coverage and composition. In contrast,  our method used progressing composition and coverage.
•	Our method is a reference-independent algorithm without using phylogenetic affiliations for binning. Similarity is only applied to identify SCGs.

•	Small scaffolds are binned by Naïve Bayesian classifier

## Perl Packages ##

Getopt::Long

File::Basename

FindBin

lib

Statistics::Distributions

Cwd

File::Path

threads

threads::shared

POSIX

The above packages are needed be installed in order to run SeedBin. They can downloaded from http://www.cpan.org/. Please note that Statistics packages are included in the same directory of SeedBin. 

## Cutoffs for seeding and mergeing ##
The file named cutoff.xls must be localized at the same directory of SeedBin. it contains the cutoffs for seeding and merging.

## Usage ##
The following arguments must be provided:
  -fafile \: fasta file containing scaffolds
  -depthfile\: the file containing depths for scaffolds
  -tablefile \<s\>\: the file containing SCG type and their Scaffolds
  
## Coverage calculation ##
After assembly we map the reads of each sample back to the assembly using soap (http://soap.genomics.org.cn/). Of course, you can use other mapping software such as bowtie2. Then you can use soap.coverage (http://soap.genomics.org.cn/) to calculate base coverage. Finally, use Scripts/depth_mean_rawdepth.pl to calculate coverage file for -depthfile. There is one example in Scripts/example/. The first parameter for Scripts/depth_mean_rawdepth.pl is a file with format that [file path] and [trimmed length] separated by tab 
Please note the trimmed length for both end is very similar to the read length of Illimula read. 

## Generation of SCG table ##


## Support ##
If you are having issues, please email me via zhouyizhuang3@163.com


