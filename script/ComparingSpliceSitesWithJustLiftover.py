#!python
import numpy as np
import sys
import re

myargs = sys.argv

#
# This script takes two known-splicesite-infile and novel-splicesite-outfile as input,
# try to figure out which known-splicesite-infile is more similar with novel-splicesite-outfile
# known-splicesite-infile was used as input for RNA-seq mapping using hisat2.
# novel-splicesite-outfile is the output of hisat2.
# Here I implemented this script to compare the output of function gffCoordinateLiftOver and the output of function annotationAndExonerateAndNovo
# For this aim I only consider the splice site of CDS elements
#

print ("usage: python ComparingSpliceSitesWithJustLiftover.py known-splicesite-infile1 known-splicesite-infile2 novel-splicesite-outfile");

print ("you want to check the overlap between " + myargs[0] + " and " + myargs[2] + ", and check the overlap between " +myargs[1] + " and "+myargs[2])

def strandToInt(strand):
    if( strand == "+" ):
        return 1
    elif (strand == "-"):
        return 0
    else:
        return 2 # there is something wrong

def intToStrand(strand):
    if( strand == 1 ):
        return "+"
    elif (strand == 0):
        return "-"
    else:
        return "#" # there is something wrong

def chrToInt(preFix, chr):
    chr = re.sub(preFix, r"", chr)
    return chr

def readSpliceSitesFile(spliceSitesFile):
    ss = np.empty([0,4], int)
    with open(spliceSitesFile) as f:
        for line in f:
            line = re.sub(r"\n", "", line)
            elements = line.split("\t")
            chr = elements[0]
            chr = re.sub(r"Chr", "", chr)
            intStrand = strandToInt(elements[3])
            chr = int(chr)
            intStrand = int(intStrand)
            ss = np.append(ss, np.array([[chr, elements[1], elements[2], intStrand]]), axis=0)
    return ss

def compareTwoList(list1, list2):
    overlap = 0
    for l1 in list1:
        for l2 in list2:
            if ( l1[0]==l2[0] and l1[1]==l2[1] and l1[2]==l2[2] and l1[3]==l2[3] ):
                ++overlap
    return overlap



file1 = myargs[0]
file2 = myargs[1]
nfile = myargs[2]

splice1 = readSpliceSitesFile("/Users/song/PA10000.novelSpliceSite")
splice2 = readSpliceSitesFile("/Users/song/PA10000.novelSpliceSite")
spliceNovel = readSpliceSitesFile("/Users/song/PA10000.novelSpliceSite")

size1 = len(np.atleast_1d(splice1))
size2 = len(np.atleast_1d(splice2))

overlap1 = compareTwoList(splice1, spliceNovel)
overlap2 = compareTwoList(splice2, spliceNovel)

print ("size1: " , splice1.shape)
print ("size2: " , splice2.shape)
print ("overlap1: " , overlap1)
print ("overlap2: " , overlap2)
