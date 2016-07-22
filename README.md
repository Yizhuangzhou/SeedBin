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

The above packages are needed be installed in order to run SeedBin. They can downloaded from http://www.cpan.org/

## Support ##
If you are having issues, please email me via zhouyizhuang3@163.com


