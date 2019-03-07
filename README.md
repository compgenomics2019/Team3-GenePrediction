# Team3-GenePrediction

## command for aragorn
```
wget http://mbio-serv2.mbioekol.lu.se/ARAGORN/Downloads/aragorn1.2.38.tgz
tar -xvzf aragorn1.2.38.tgz
cd aragorn1.2.38
gcc -O3 -ffast-math -finline-functions -o aragorn aragorn1.2.38.c
./aragorn -t ../contigs.fasta 
```
## command for rnammer
### install hmmer

```
wget http://eddylab.org/software/hmmer/hmmer-2.2g.tar.gz
uncompress hmmer-2.2.tar.gz
tar zxf hmmer-2.2g.tar.gz
cd hmmer-2.2g/
./configure --host=host\_name
make
make check
```
### install rnammer
_the link is available for 4 hours_
wget http://www.cbs.dtu.dk/download/4AE7BE96-3ED6-11E9-8B58-A1B3B9CD16B5/rnammer-1.2.src.tar.Z
cat rnammer-1.2.tar.Z | gunzip | tar xvf -
cd rnammer-1.2.src
perl rnammer -S bac -m lsu,ssu,tsu -gff - < example/ecoli.fsa

_modify the PATH is the src_
my $INSTALL\_PATH = "install\_src\_path";
	$HMMSEARCH_BINARY = "hammer\_path/hmmer-2.2g/binaries/hmmsearch";
	$PERL = "perl\_path";
_test_
perl rnammer-1.2.src/rnammer -S bac -m lsu,ssu,tsu -multi -gff rnammer.gff -f rnammer.fasta < contigs.fasta


## Bacterial Genome Prediction

This pipeline is meant to predict the gene prediction using genome. Gene prediction is the process of finding which regions of genomic DNA encodes genes. Gene prediction not only predicts the DNA encoding but also protein coding genes like RNA genes.  

## Installing Miniconda

First install the bash file for installing Miniconda. Miniconda is comfortable to install different tools. 

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh
rm  Miniconda3-latest-Linux-x86_64.sh
export PATH=~/miniconda/bin:$PATH
```

## Miniconda Create environment from yml file

<!-- Strong -->

You don't need to worry about downloading individual tools. We already made a environment file where if you just download the yml file, it will download the tools for you. 

```
#Create environment after downloading yml file(on google drive)
conda-env create -f installs/environment2.yml -n gp
source activate gp
```

## Install Genemark-S2 and command and parameters for using gen mark

<!-- Links -->
[Gene Mark](exon.gatech.edu/)

You need to download Genemark-S2 from the source. The website is from exon.gatech.edu.  You can click on Gene Mark. 
