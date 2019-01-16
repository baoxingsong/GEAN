#!perl -w
use strict;

#===============================================================================
#
#         FILE:  mafftSubmit.pl
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
#      VERSION:  1.1
#      CREATED:  09/17/2017 06:15:20 AM
#     REVISION:  07/20/2018
#     REVISION:  compare the sequence before MSA and after MSA, If they are different, then realign them
#===============================================================================

my $count=0;
my $command;
opendir(DIR1, "../") || die "Can't open directory ../";
my @dots1 = readdir(DIR1);
foreach my $file1 (@dots1){
	my $chrName = $file1;
	if( ($chrName=~/\w/) && ( $chrName ne "mafft" ) && ( $chrName ne "command") ){
		my $c = "";
		system("mkdir $chrName");
		my $dir_name = "./" . $chrName . "/";
		my $dir_name2 = "../" . $chrName . "/";
		opendir(DIR, $dir_name2) || die "Can't open directory $dir_name";
		my @dots = readdir(DIR);
		foreach my $file (@dots){
			if($file=~/^\d+_\d+$/){
				my $inputfile=$dir_name2 . $file;
				my $outfile=$dir_name . $file . ".mafft";
				my $iftoberun = 1;

				if(-e $outfile){
					my @args = stat ($outfile);
					my $size = $args[7];
					if( $size > 1 ){
						$iftoberun = 0;
						my %before_msa_seqs;
						my %after_msa_seqs;
						my $name="";
						my $seq="";
						open INPUT, "$inputfile";
                        while( my $line=<INPUT> ){
                            if( $line=~/>(\S+)/  ){
                                if( length($name)>0 && length($seq)>0  ){
                                    $seq=~s/\s//g;
                                    $seq=~s/-//g;
                                    $seq = uc($seq);
                                    $before_msa_seqs{$name} = $seq;
                                }
                                $name=$1;
                                $seq="";
                            }else{
                                $seq = $seq . $line;
                            }
                        }
                        if( length($name)>0 && length($seq)>0  ){
                            $seq=~s/\s//g;
                            $seq=~s/-//g;
                            $seq = uc($seq);
                            $before_msa_seqs{$name} = $seq;
                        }
						close INPUT;
						$name="";
						$seq="";
                        open INPUT, "$outfile";
                        while( my $line=<INPUT> ){
                            if( $line=~/>(\S+)/  ){
                                if( length($name)>0 && length($seq)>0  ){
                                    $seq=~s/\s//g;
                                    $seq=~s/-//g;
                                    $seq = uc($seq);
                                    $after_msa_seqs{$name} = $seq;
                                }
                                $name=$1;
                                $seq="";
                            }else{
                                $seq = $seq . $line;
                            }
                        }
                        if( length($name)>0 && length($seq)>0  ){
                            $seq=~s/\s//g;
                            $seq=~s/-//g;
                            $seq = uc($seq);
                            $after_msa_seqs{$name} = $seq;
                        }
						close INPUT;
						while (my ($key, $value) = each(%before_msa_seqs)) {
                            if( (!exists( $after_msa_seqs{$key} )) || ( $after_msa_seqs{$key} ne $value ) ){
                                $iftoberun = 1;
                            }
                        }
					}
				}
				
				if(1== $iftoberun){
					print "mafft --auto $inputfile > $outfile\n";
				}
			}
		}
		close DIR;
	}
}
