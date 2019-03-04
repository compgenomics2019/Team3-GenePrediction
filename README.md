# Team3-GenePrediction

## command for aragorn
wget http://mbio-serv2.mbioekol.lu.se/ARAGORN/Downloads/aragorn1.2.38.tgz
tar -xvzf aragorn1.2.38.tgz
cd aragorn1.2.38
gcc -O3 -ffast-math -finline-functions -o aragorn aragorn1.2.38.c
./aragorn -t ../contigs.fasta 

## command for rnammer
### install hmmer
wget http://eddylab.org/software/hmmer/hmmer-2.2g.tar.gz
uncompress hmmer-2.2.tar.gz
tar zxf hmmer-2.2g.tar.gz
cd hmmer-2.2g/
./configure --host=host\_name
make
make check

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
