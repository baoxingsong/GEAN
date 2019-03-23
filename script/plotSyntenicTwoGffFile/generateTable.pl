#!perl -w
use strict;
use warnings FATAL => 'all';

my $refGff = $ARGV[0];
my $newGff = $ARGV[1];

my %gene_chr;
my %gene_position;

open INPUT, "$refGff";
while( my $line=<INPUT> ){
    if( $line=~/^([\dMtP]+)\s+\S+\s+gene\s+(\d+).*ID=(.*?);/ ){
        $gene_chr{$3}=$1;
        $gene_position{$3}=$2;
    }
}
close INPUT;

open INPUT, "$newGff";
while( my $line=<INPUT> ){
    if( $line=~/^([\dMtP]+)\s+\S+\s+gene\s+(\d+).*ID=(.*?);/ ){
        if( exists $gene_chr{$3} ){
            print "$gene_chr{$3}\t$gene_position{$3}\t$1\t$2\n";
        }
    }
}
close INPUT;
