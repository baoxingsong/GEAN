# Introduction

Here, we provide example data to help go throw the genetic variants synchronization pipeline taking advantage of the gene structure annotation. You need to modify some commands or scripts to run your own data.

## Preparation
````
cd input
cp ../denovoVariantCalling/tair10.fa.gz ./
cp ../denovoVariantCalling/TAIR10_GFF3_genes_no_UM.gff.gz ./
gunzip *gz
````
## Preparation
Prepare a file named as vList. The first column is the accession name, and the second column is the path to the sdi file. The following content works for the example data on my computer.
````
PA9996  /media/song/8t2/geneStructureAnnotation/MSA/input/PA9996.sdi
PA9997  /media/song/8t2/geneStructureAnnotation/MSA/input/PA9997.sdi
PA9998  /media/song/8t2/geneStructureAnnotation/MSA/input/PA9998.sdi
PA9999  /media/song/8t2/geneStructureAnnotation/MSA/input/PA9999.sdi
````
Run the following commands to get the chromosome name list
````
grep ">" tair10.fa | sed '~s/>//g' > chrList
````
# Run the pipeline step by step
## Cut the whole genome sequence into fragments
````
gean premsa -i TAIR10_GFF3_genes_no_UM.gff -r tair10.fa -v vList -t 8 -l 100
````
## Perform MSA on each fragments
````
mkdir MsaPreFolder/mafft
cd MsaPreFolder/mafft
perl ../../script/mafftSubmit_V1.1.pll > command
````
You could run command line by line with any way you like.
For example you could run it with `xjobs` in parallel
````
xjobs -j 10 -s command
````
*For step-by-step style, you could use any multiple-sequence-alignment software for this step, but you should name the MSA result with "original name" + ".mafft", and put them into the folders named by chromosome names.
## Create new sdi files
Create a new sdi file for each accession
````
~/Dropbox/gean/gean msatosdi -c chrList -m ./ -o ./ -r ../../tair10.fa -t 8 -v ../../vList
````
*If you do NOT want to generate sdi files for all the chromosome or all the accessions, you could modify file chrList or ../../vList.\
You will find the output sdi files with synchronized variants under the current folder.
## Reformat sdi files into PLINK files.
please go to [Irasis](https://github.com/baoxingsong/Irisas) for this step.
