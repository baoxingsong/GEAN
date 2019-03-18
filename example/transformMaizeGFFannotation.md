Here I am showing how to transform the GFF file from B73 v3 reference to the M017 CAU assemblies.
In this example the minimap2 splice aware sequence alignment model is used, which is much much faster than Mummer

Download the files
````
# B73 genomic sequence
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-31/fasta/zea_mays/dna/Zea_mays.AGPv3.31.dna.genome.fa.gz

# GFF3 file for B73 V3 genome
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-31/gff3/zea_mays/Zea_mays.AGPv3.31.gff3.gz

# Genomic sequences Mo17
https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-Mo17-REFERENCE-CAU-1.0/Zm-Mo17-REFERENCE-CAU-1.0.fa.gz

# uncompress all the files
gunzip Zea_mays.AGPv3.31.dna.genome.fa.gz
gunzip Zea_mays.AGPv3.31.gff3.gz
gunzip Zm-Mo17-REFERENCE-CAU-1.0.fa.gz
````

get the transcripts sequence of B73
````
gean gff2seq -i Zea_mays.AGPv3.31.gff3 -r Zea_mays.AGPv3.31.dna.genome.fa -p b73_protein.fa -c b73_cds.fa -g b73_gene.fa
````

to save computational time, only the primary transcripts sequence are used
````
perl script/getPrimaryTranscriptSequenceFromGeanOutput.pl Zea_mays.AGPv3.31.gff3 b73_cds.fa >primary_cds.fasta
````

align the primary transcript CDS sequence to the Mo17 genome sequencing using minimap2
````
minimap2 -ax splice -t 90 -a -uf -C5 Zm-Mo17-REFERENCE-CAU-1.0.fa primary_cds.fasta >cds.sam
````
begin to run gean
````
split -l 20 cds.sam minimap2
ls minimap2* | awk '{print "gean spltogff -i Zea_mays.AGPv3.31.gff3 -l 100000 -w 5000 -r Zea_mays.AGPv3.31.dna.genome.fa -a "$1" -s Zm-Mo17-REFERENCE-CAU-1.0.fa -o "$1".gff"}' > spltogff_commnd
xjobs -j 50 -s spltogff_commnd
````

then merge those files together
````
cat minimap*gff > mo17.gff
````
do some duplication filtering
````
gean purifygff -i mo17.gff -s Zm-Mo17-REFERENCE-CAU-1.0.fa -o mo17_purify.gff
````
keep only very confident non-redundant records by syntenic analysis
````
gean sinsyn -i b73.gff -s Mo17.fa -a mo17.gff -o Mo17_sinsyn.gff
````
keep only those transcript records only with ORF conserved
````
gean orf -i Mo17_sinsyn.gff -s Mo17.fa -o Mo17_sinsyn_orf.gff
````
###### this is pipeline to transform the gene structure annotation from reference genome sequence to another de novo assembly genome squence. This is not a complete pipeline to annotate the new assembly.
###### If you are doing genome assembly for a new accession using different genome sequence assembly pipelines. The number of successfully transformed transcripts (marked as ConservedFunction in the output gff file) could be used to measure the genome assembly quality.
