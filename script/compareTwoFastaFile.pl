#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  compareTowFastaFile.pl
#
#        USAGE:  ./compareTowFastaFile.pl  
#
#  DESCRIPTION:  
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  Baoxing Song (songbx.me), song@mpipz.mpg.de
#      COMPANY:  MPIPZ
#      VERSION:  1.0
#      CREATED:  09/14/2017 01:35:30 PM
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;

my %seqs1;
my $name1="";
my $seq1="";
open INPUT, "$ARGV[0]";
while( my $line=<INPUT>  ){
    if( $line=~/>(\S+)/  ){
        if( length($name1) > 0  ){
            $seq1=~s/\s//g;
            $seq1 = uc ($seq1);
            $seqs1{$name1} = $seq1;
        }
        $name1=$1;
        $seq1="";
    }else{
        $seq1 = $seq1 . $line;
    }
}
close INPUT;
if( length($name1) > 0  ){
    $seq1=~s/\s//g;
    $seq1 = uc ($seq1);
    $seqs1{$name1} = $seq1;
}

my %seqs2;
my $name2="";
my $seq2="";
open INPUT, "$ARGV[1]";
while( my $line=<INPUT>  ){
    if( $line=~/>(\S+)/  ){
        if( length($name2) > 0  ){
            $seq2=~s/\s//g;
            $seq2 = uc ($seq2);
            $seqs2{$name2} = $seq2;
        }
        $name2 = $1;
        $seq2 = "";
    }else{
        $seq2 = $seq2 . $line;
    }
}
close INPUT;

if( length($name2) > 0  ){
    $seq2=~s/\s//g;
    $seq2 = uc($seq2);
    $seqs2{$name2} = $seq2;
}

while( my ($k, $v) = each %seqs1  ){
    if( $seqs2{$k} eq $v  ){

    }else{
        print "$ARGV[0]\t$k\n";
    }
}
while( my ($k, $v) = each %seqs2  ){
    if( $seqs1{$k} eq $v  ){

    }else{
        print "$ARGV[0]\t$k\n";
    }
}

