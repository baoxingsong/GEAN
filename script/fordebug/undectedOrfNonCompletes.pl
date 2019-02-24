#!perl
use strict;
use warnings FATAL => 'all';

my %allTranscripts;
my %allGoodTranscripts;
my %liftTranscripts;
my %realignTranscripts;
my $source="";
my $transcript="";

open INPUT, "/workdir/bx674/melon_debug/bs674/melon_gean_purify.gff";
while( my $line=<INPUT>  ){
    if( $line=~/^(\S+)\t+(\S+)\t+.*Parent=([0123456789MELOC]+\.\d\.\d)/  ){
        $source=$2;
        $transcript=$3;
    }elsif( $line=~/ConservedFunction/  ){
        if( $source eq "LIFTOVER"  ){
            $liftTranscripts{$transcript}=1;
        }elsif( $source eq "REALIGNMENT"  ){
            $realignTranscripts{$transcript}=1;
        }
        $allGoodTranscripts{$transcript}=1;
    }
    $allTranscripts{$transcript}=1;
}
close INPUT;

my %lostOfFunctions;
my %withFunctions;
open INPUT, "/workdir/bx674/melon/gean/CM3.6.1_pseudomol_cds.fa";
while( my $line=<INPUT>  ){
    if( $line=~/^>(\S+)\s/ ){
        $transcript=$1;
        if( $line=~/ConservedFunction/ ){
            $withFunctions{$transcript}=1;
        }else{
            $lostOfFunctions{$transcript}=1;
        }
    }
}
close INPUT;

print "those transcripts have re-annotated as functional\n";

while( my ($key, $value) = each %allGoodTranscripts ){
    if( exists $lostOfFunctions{$key} ){
        print "$key\n";
    }
}

print "those functional transcript have not been aligned correctly\n";

while( my ($key, $value) = each %withFunctions ){
    if( exists $allGoodTranscripts{$key} ){

    }else{
        print "$key\n";
    }
}


my $size = keys %allTranscripts;
print "allTranscripts\t$size\n";
$size = keys %allGoodTranscripts;
print "allGoodTranscripts\t$size\n";
$size = keys %realignTranscripts;
print "realignTranscripts\t$size\n";
$size = keys %liftTranscripts;
print "liftTranscripts\t$size\n";

