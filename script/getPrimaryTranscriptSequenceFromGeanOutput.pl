#!perl
use strict;
use warnings FATAL => 'all';

my $seq="";
my $name="";
my %seqs;
open INPUT, "$ARGV[1]";
while( my $line=<INPUT> ){
    if( $line=~/^>(\S+)/ ){
        if( length($seq)>0 && length($name)>0 ){
            $seqs{$name}=$seq;
        }
        $name=$1;
        $seq="";
    }else{
        $seq = $seq . $line;
    }
}
if( length($seq)>0 && length($name)>0 ){
    $seqs{$name}=$seq;
}
close INPUT;


my %geneToPrimaryTrasncript;
open INPUT, "$ARGV[0]";
while( my $line=<INPUT> ){
    if( $line=~/^(\S+)\s+(\S+)\s+(mRNA|transcript)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+.*ID=(\S+?);.*Parent=(\S+?)(;|$)/ ){
        if( exists $geneToPrimaryTrasncript{$10} ){

        }elsif ( exists $seqs{$9} ) {
            $geneToPrimaryTrasncript{$10}=$9;
        }
    }elsif( $line=~/^(\S+)\s+(\S+)\s+(mRNA|transcript)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+.*Parent=(\S+?);.*ID=(\S+?)(;|$)/  ){
        if( exists $geneToPrimaryTrasncript{$9}  ){

        }elsif ( exists $seqs{$10}  ) {
            $geneToPrimaryTrasncript{$9}=$10;
        }
    }
}
close INPUT;



while( my ($key, $value) = each %geneToPrimaryTrasncript ){
    print ">$value\n";
    print "$seqs{$value}";
}
