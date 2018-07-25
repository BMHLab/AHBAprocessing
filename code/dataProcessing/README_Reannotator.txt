########################################
# HOW TO GET STARTED WITH RE-ANNOTATOR #
########################################

A. Things to do before the first time using Re-Annotator

########################################
I) Install External Programs
########################################

Following programs should be installed on your system prior to running the Re-Annotator:
a) PERL
it is in Mac
b) BWA        (http://sourceforge.net/projects/bio-bwa/files/)
to install unzip the file and run
make
After copy executable file in the directory you will be running things from and put this path to config.sh script
c) SAMtools   (https://sourceforge.net/projects/samtools/files/)
to install unzip the file and run
./configure
make
d) Annovar    (http://www.openbioinformatics.org/annovar/)
just download and unzip

########################################
II) Get External Data
########################################

0) make sure you have probe sequences. Agilent allows to download this information from: https://earray.chem.agilent.com/earray/
Need to choose 014850 ID and download FASTA file.
It will have an extension .txt, so rename it to have .fasta at the end
I named mine probes2annotate.fasta and put it in the easy to access folder (probably where the Re-annotate is)
CUST sequences were sent by AHBA after messaging them. The file is re-formated to .txt format from excel and saved as probes2annotateALL.txt
Then in order to properly convert txt to fasta format thet works with re-annotator go to:http://iubio.bio.indiana.edu/cgi-bin/readseq.cgi
and input this file, choose output sequence format as Pearson|FASTA|fa and do the conversion. Name the output file appropriately probes2annotateALL.fasta


a) Reference Genome Sequence (chrom folder)
- Download the reference genome, e.g., hg38
http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz (this will be downloaded in step IIc fro the command line)
 - unzip (using gzip)
  make sure every chromosome is in a single file

b) Gene Database
- Download information on gene locations, e.g., RefSeq (this will also be downloaded from the command line in step IIc)
http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz
- unzip (using gzip)

c) Make sure desired databases for annotation are available in Annovar:
These can be for instance
- RefGene
./annotate_variation.pl -downdb refGene -buildver hg38 humandb/
running this command will ask to run two additional commands to generate fa files:

./annotate_variation.pl --buildver hg38 --downdb seq humandb/hg38_seq
./retrieve_seq_from_fasta.pl humandb/hg38_refGene.txt -seqdir humandb/hg38_seq/chroms/ -format refGene -outfile humandb/hg38_refGeneMrna.fa

then copy this file hg38_refGeneMrna.fa to ReAnnotator/wholeRef/ to run further command (described below)
- snpdatabase, e.g., snp129 (available from annovar) - I do not run those, not using SNP information.
./annotate_variation.pl -downdb snp129 -buildver hg19 -webfrom annovar humandb/

########################################
III) Generate the mRNA Reference Sequence
########################################

a) Execute new_exomeBuilding.pl to build the mRNA reference sequence
example:
$> ./BuildExomeReference.pl -i ~/ReAnnotator/refGene.txt -o ~/ReAnnotator/exomeRef/hg19exome -r ~/ReAnnotator/hg19/

IN case when Hg38 is used download chroms from Download reference genome sequence: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/ (file: hg38.chromFa.tar.gz) - same as downloaded by annovar and saves in humandb/hg38_seq/chroms/

./BuildExomeReference.pl -i refGene.txt -o exomeRef/hg38exome -r ../chroms/

b) Use "BWA index" to generate the BWA index files for
i) the exome reference sequence: ./bwa index ReAnnotator/exomeRef/hg38exome_inclUTR.fasta
ii) the whole genome reference sequence: ./bwa index ReAnnotator/wholeRef/hg38_refGeneMrna.fa

Re-annotator expects file named hg38_refGeneMrna.fasta.gz so, rename file in a text editor to .fasta instead of fa and gzip the file: (have not checked if works with unzipped, but probably will)
gzip hg38_refGeneMrna.fasta

########################################
IV) Complete the config.sh
########################################

* provide exact locations to the external programs
in the file config.sh
* config.sh is also the place to change the default settings
 - # of CPUs
 - # of mismatches
 - genome version, e.g., gh19, hg18, mm9, ...
 - genedb (refGene, ensemble, ...)
 - snpdb  (snp135, snp136, ...)

This is how my config looked:
—————————————————————————————
#!/bin/bash

##put here the local paths to your external executables
#e.g.:
##BWA
bwaex="/Volumes/AA_HD/reannotate/bwa"

##SAMTOOLS
samex="/Volumes/AA_HD/reannotate/samtools"

##ANNOVAR
annoex="/Volumes/AA_HD/reannotate/annovar/annotate_variation.pl"

##ANNOVAR DB
annodb="/Volumes/AA_HD/reannotate/annovar/humandb/"

##set these to meet your requirments
##NOTE: genome version, and ideally annogenedb should match the
##      versions used for aligning against the whole genome the
##      versions used for building the exome, respectively

buildver="hg38"
annogenedb="refGene"

#provide here general settings
##number of CPUs for alignment
ncpu=4

##maximal number of mismatches in alignment
maxmis=5

##query length
qlen=50

##do SNPannotation
doSNPanno=false

export bwaex
export samex
export annoex
export annodb
export annogenedb
export annosnpdb
export ncpu
export maxmis
export qlen
export buildver
—————————————————————————————

NOTE: step III) has to be carried out every time the Gene Database is updated!
NOTE: update config.sh to meet the needs for your Re-Annotation

B. Things to do at every run

I) Check that config.sh is set up correctly for the genome and exome you are about to use
 - regenerate the mRNA reference if you are using an updated database version or new genome release

II) Convert probes file into a fasta file
For Illumina probe files the script parse_fastaFromOriginIlmnAnno.pl can be used.

########################################
execute: ./ReAnnotator.sh
with corresponding parameters
########################################

example:
./ReAnnotator.sh probes2annotateALL.fasta exomeRef/hg38exome_inclUTR.fasta wholeRef/hg38_refGeneMrna.fasta.gz refGene.txt outputs tmp

When running the line, make sure the output folder is empty to begin with, otherwise it will use old files to generate the output (and the same result will be generated each time no matter the input).

Prior to running the script, make sure that all the directories exist.

The script will call following scripts in a row:
 a) run_realignment.sh
 b) run_coordinateConversion.sh
 c) run_genome_snp_annotation.sh
