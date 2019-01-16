#!perl -w
use strict;
use warnings FATAL => 'all';

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
#      VERSION:  1.0
#      CREATED:  09/17/2017 06:15:20 AM
#     REVISION:  ---
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
				my $iftoberun = 0;
				if(-e $outfile){
					my @args = stat ($outfile);
					my $size = $args[7];
					if( $size > 1 ){
                        open INPUTNON, "$inputfile";
                        my $name_non="";
                        my $seq_non="";
                        my %seqs_non;
                        while( my $line_non=<INPUTNON>  ){
                            if( $line_non=~/>(\S+)/  ){
                                if( length($name_non)>0 && length($seq_non)>0  ){
                                    $seq_non=~s/\s//g;
                                    $seq_non=~s/-//g;
                                    $seq_non = uc($seq_non);
                                    $seqs_non{$name_non} = $seq_non;
                                }
                                $name_non=$1;
                                $seq_non="";
                            }else{
                                $seq_non = $seq_non . $line_non;
                            }
                        }
                        if( length($name_non)>0 && length($seq_non)>0  ){
                            $seq_non=~s/\s//g;
                            $seq_non=~s/-//g;
                            $seq_non = uc($seq_non);
                            $seqs_non{$name_non} = $seq_non;
                        }
                        close INPUTNON;
                        
                        open INPUTMSA, "$outfile";
                        my $name_msa="";
                        my $seq_msa="";
                        my %seqs_msa;
                        while( my $line_msa=<INPUTMSA>  ){
                            if( $line_msa=~/>(\S+)/  ){
                                if( length($name_msa)>0 && length($seq_msa)>0  ){
                                    $seq_msa=~s/\s//g;
                                    $seq_msa=~s/-//g;
                                    $seq_msa = uc($seq_msa);
                                    $seqs_msa{$name_msa} = $seq_msa;
                                }
                                $name_msa=$1;
                                $seq_msa="";
                            }else{
                                $seq_msa = $seq_msa . $line_msa;
                            }
                        }
                        if( length($name_msa)>0 && length($seq_msa)>0  ){
                            $seq_msa=~s/\s//g;
                            $seq_msa=~s/-//g;
                            $seq_msa = uc($seq_non);
                            $seqs_msa{$name_msa} = $seq_msa;
                        }
                        close INPUTNON;
                        while( my ($k_non, $v_non)=each %seqs_non  ){
                            if( (exists $seqs_msa{$k_non}) && ($v_non eq $seqs_msa{$k_non})  ){

                            }else{
                                $iftoberun = 1;
                            }
                        }
                    }else{
                        $iftoberun = 1;
                    }
                }else{
                    $iftoberun = 1;
                }
				
				if(1== $iftoberun){
					print "mafft --auto $inputfile > $outfile\n";
				}
			}
		}
		close DIR;
	}
}
