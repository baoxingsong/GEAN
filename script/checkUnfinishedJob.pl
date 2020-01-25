#!perl -w
use strict;
use warnings FATAL => 'all';

open INPUT, "$ARGV[0]";
while( my $line=<INPUT> ){
	if( $line=~/\s(\S+gff)$/ ){
		my $outputFile = "$1";
		my $good = 0;
		open INPUT0, "$outputFile";
		while( my $line0=<INPUT0> ){
			if( $line0=~/liftover done/ ){
				$good = 1;
			}
		}
		close INPUT0;
		if( $good == 0 ){
			print "$line";
		}
	}
}
close INPUT;
