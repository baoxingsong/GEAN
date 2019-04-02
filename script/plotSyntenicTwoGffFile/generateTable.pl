#!perl -w
use strict;
use warnings FATAL => 'all';

my $refGff = $ARGV[0];
my $newGff = $ARGV[1];

my %gene_chr;
my %gene_position;

open INPUT, "$refGff";
while( my $line=<INPUT> ){
    if( $line=~/^(\S+)\s+\S+\s+gene\s+(\d+).*ID=(.*?);/ ){
        $gene_chr{$3}=$1;
        $gene_position{$3}=$2;
    }
}
close INPUT;

my $group="";
open INPUT, "$newGff";
while( my $line=<INPUT> ){
    if( $line=~/^#start chain (\d+)/){
        $group = $1;
    }
    if( $line=~/^(\S+)\s+\S+\s+gene\s+(\d+).*ID=(.*?);/ ){
        my $chr=$1;
        my $position=$2;
        my $gene=$3;
        $gene=~s/_\d+$//g;
        if( exists $gene_chr{$gene} && $position < 2147483646 ){
            print "$gene_chr{$gene}\t$gene_position{$gene}\t$chr\t$position\t$group\n";
        }else{
            print STDERR "$gene\n";
        }
    }
}
close INPUT;
