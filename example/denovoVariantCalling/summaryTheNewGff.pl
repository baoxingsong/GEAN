#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  summaryTheNewGff.pl
#
#        USAGE:  ./summaryTheNewGff.pl  
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
#      CREATED:  09/16/2018 05:31:19 PM
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;

my %allPossiableMrna;
open INPUT, "ler.gff";
while( my $line=<INPUT>  ){
    if( $line=~/CDS.*Parent=(AT\dG\d\d\d\d\d\.\d+)/  ){
        $allPossiableMrna{$1}=1;
    }
}
close INPUT;


my %allTranscripts;
my %allGoodTranscripts;
my %liftTranscripts;
my %realignTranscripts;
my $source="";
my $transcript="";
open INPUT, "ler_purify_v2.gff";
while( my $line=<INPUT>  ){
    if( $line=~/^(\S+)\t+(\S+)\t+.*Parent=(AT\dG\d\d\d\d\d\.\d+)/  ){
        $source=$2;
        $transcript=$3;
    }elsif( $line=~/ConservedFunction/ && (exists $allPossiableMrna{$transcript})  ){
        if( $source eq "LIFTOVER"  ){
            $liftTranscripts{$transcript}=1;
        }elsif( $source eq "REALIGNMENT"  ){
            $realignTranscripts{$transcript}=1;
        }
        $allGoodTranscripts{$transcript}=1;
    }
    if( exists $allPossiableMrna{$transcript}  ){
        $allTranscripts{$transcript}=1;
    }
}
close INPUT;

my $size = keys %allTranscripts;
print "allTranscripts\t$size\n";
$size = keys %allGoodTranscripts;
print "allGoodTranscripts\t$size\n";
$size = keys %realignTranscripts;
print "realignTranscripts\t$size\n";
$size = keys %liftTranscripts;
print "liftTranscripts\t$size\n";

