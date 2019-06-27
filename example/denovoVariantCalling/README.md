Here I am showing how to project the TAIR10 reference genome sequence gene structure annotation to a de novo assembly Ler-0 genome sequence.
Then base pair resolution whole genome sequence alignment and variant calling would be performed.
The TAIR10 genome sequence and gene structure annotation were downloaded from the TAIR website (https://www.arabidopsis.org/).
The Ler-0 de novo genome sequence was released by reference publication 1.
## Find highly similar fragments
The purpose of this step is to provide anchors for high resolution alignment. Results from some seed-to-chain algorithm is enough.
I am using mummer-4.0.0beta2.

```
nucmer -t 40 --maxmatch -c 100 -b 500 -l 50 -p t2l tair10.fa ler.fa
show-coords -H -T -r -l -c t2l.delta >t2l.tab
cat t2l.tab | awk '$12 == $13 {print $12"\t"$1"\t"$2"\t"$12"\t"$3"\t"$4}' > t2l.tbl
```
for parallelization purpose, split the t2l.tbl into small files
```
head -10000 t2l.tbl >align7.tbl
head -20000 t2l.tbl | tail -10000 >align8.tbl
head -30000 t2l.tbl | tail -10000 >align9.tbl
tail -3100 t2l.tbl > align10.tbl
```
prepare the align1.tbl, align2.tbl, align3.tbl, align4.tbl, align5.tbl manually.
align1.tbl is based on the result from reference publication 1 and it is only important for annotation purpose

## Liftover the reference gene structure annotation to the de novo assembly genome sequence
run the following commands in parallel
```
gean transgff -i TAIR10_GFF3_genes_no_UM.gff -r tair10.fa -a align1.tbl -s ler.fa -o ler1.gff -w 50000 -sl
gean transgff -i TAIR10_GFF3_genes_no_UM.gff -r tair10.fa -a align2.tbl -s ler.fa -o ler2.gff -w 50000 -sl
gean transgff -i TAIR10_GFF3_genes_no_UM.gff -r tair10.fa -a align3.tbl -s ler.fa -o ler3.gff -w 50000 -sl
gean transgff -i TAIR10_GFF3_genes_no_UM.gff -r tair10.fa -a align4.tbl -s ler.fa -o ler4.gff -w 50000 -sl
gean transgff -i TAIR10_GFF3_genes_no_UM.gff -r tair10.fa -a align5.tbl -s ler.fa -o ler5.gff -w 50000 -sl
gean transgff -i TAIR10_GFF3_genes_no_UM.gff -r tair10.fa -a align6.tbl -s ler.fa -o ler6.gff -w 50000 -sl
gean transgff -i TAIR10_GFF3_genes_no_UM.gff -r tair10.fa -a align7.tbl -s ler.fa -o ler7.gff -w 50000 -sl
gean transgff -i TAIR10_GFF3_genes_no_UM.gff -r tair10.fa -a align8.tbl -s ler.fa -o ler8.gff -w 50000 -sl
gean transgff -i TAIR10_GFF3_genes_no_UM.gff -r tair10.fa -a align9.tbl -s ler.fa -o ler9.gff -w 50000 -sl
gean transgff -i TAIR10_GFF3_genes_no_UM.gff -r tair10.fa -a align10.tbl -s ler.fa -o ler10.gff -w 50000 -sl
```
Merge those annotations files and remove duplication gene annotations
```
cat ler*.gff >ler.gff
gean purifygff -i ler.gff -s ler.fa -o ler_purify.gff
```
## Whole genome sequence alignment and variant calling using the gene structure annotation as anchors
```
gean varcall -i TAIR10_GFF3_genes_no_UM.gff -r tair10.fa -t ler_purify.gff -s ler.fa -w 5000 -o zsdp_ler.sdi
```
## For debugging
You could generate a new genome sequence by replacing the reference allele with alternative allele taking the reference
genome sequence and variant calling files as input using the 'pseudogeno' function. And then check the identical of the de novo assembly and the new generated genome sequence file.
```
gean pseudogeno -r tair10.fa -v zsdp_ler.sdi -o ler_zsdp_de_novo.fa
perl compareTwoFastaFile.pl ler_zsdp_de_novo.fa ler.fa
```
## Reference publications
```
1. http://www.pnas.org/content/113/28/E4052
```
