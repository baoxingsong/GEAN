#!perl -w
use strict;
use warnings FATAL => 'all';
#===============================================================================
#
#         FILE:  mergeSdi.pl
#
#        USAGE:  ./mafftSubmit.pl
#
#  DESCRIPTION:  This is a part of WMSA pipeline. It is designed to generate mafft commands.
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  Baoxing Song (songbx.me), song@mpipz.mpg.de
#      COMPANY:  MPIPZ
#      VERSION:  1.0
#      CREATED:  09/16/2017 02:12:16 AM
#     REVISION:  ---
#===============================================================================

my @chrs;
open INPUT, "../../cl";
while( my $line=<INPUT> ){
	if( $line=~/(\S+)/ ){
		push(@chrs, $1);
	}
}
close INPUT;

open INPUT, "../../al";
while( my $line=<INPUT> ){
	$line=~s/\s//g;
	if( $line=~/\w/ ){
		open OUTPUT, ">$line.sdi";
		my $lineNumber=2;
		for my $chr ( @chrs ){
			my $thissdifile="../$chr/$line.sdi";
			if( -e $thissdifile ){
				open INPUT2, "$thissdifile";
				while( my $line2=<INPUT2> ){
					$line2=~s/\s$//g;
					$line2=~s/^\s//g;
					if( $line2=~/$chr/ ){
						if( $lineNumber>2 ){
							print OUTPUT "\n$line2";
						}else{
							print OUTPUT "$line2";
						}
						$lineNumber++;
					}
				}
				close INPUT2
			}
		}
	}
}
close INPUT;

