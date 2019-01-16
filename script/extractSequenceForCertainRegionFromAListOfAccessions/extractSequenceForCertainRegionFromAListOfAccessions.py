#!python3
import subprocess
import re

# songbaoxing168@163.com

class Fasta:
    name = ""
    seq = ""
    def __init__(self, name, seq):
        self.name = name
        self.seq = seq

def readFastaFile(fastaFile):
    fastas = {}
    chromosome_names = []
    name = ""
    seq = []
    with open(fastaFile) as f:
        for line in f:
            m = re.search('^>(\S+)', line)
            if (m != None):
                if (len(name) > 0) & (len(seq) > 0):
                    s = ''.join(seq)
                    s = re.sub("\\s", "", s)
                    s = s.upper()
                    fasta = Fasta(name, s)
                    fastas[name]=(fasta)
                    chromosome_names.append(name)
                name = m.group(1)
                seq = []
            else:
                seq.append(line)
        if (len(name) > 0) & (len(seq) > 0):
            s = ''.join(seq)
            s = re.sub("\\s", "", s)
            s = s.upper()
            fasta = Fasta(name, s)
            fastas[name] = (fasta)
            chromosome_names.append(name)

    return chromosome_names, fastas



def getReverseComplementary(sequence):
    reversecomplementary=[]
    for c in  sequence[::-1]:
        if ('A' == c):
            c = 'T'
        elif ('T' == c):
            c = 'A'
        elif ('U' == c):
            c = 'A'
        elif ('C' == c):
            c = 'G'
        elif ('G' == c):
            c = 'C'
        elif ('R' == c):
            c = 'Y'
        elif ('Y' == c):
            c = 'R'
        elif ('K' == c):
            c = 'M'
        elif ('M' == c):
            c = 'K'
        elif ('B' == c):
            c = 'V'
        elif ('V' == c):
            c = 'B'
        elif ('D' == c):
            c = 'H'
        elif ('H' == c):
            c = 'D'
        reversecomplementary.append(c)

    return ''.join(reversecomplementary)


def getSubSequence(fastas, name, start, end, strand):
    # get a sequence fragment from fasta records
    start = start -1
    if start > len(fastas[name].seq):
        return ""
    if end  >  len(fastas[name].seq):
        end = len(fastas[name].seq)

    seq = (fastas[name].seq)[start:end]
    if( "+" == strand ):
        return seq
    else:
        return getReverseComplementary(seq)

fileoutput= open('/netscratch/dep_tsiantis/grp_gan/song/tsianst/GWAS/FRIORF/frigida.seq', 'w')
refStart = 4648158
refEnd = 4650582
length = refEnd-refStart+1
refStart=refStart-length
refEnd=refEnd+length
if refStart<0:
    refStart=0

chr = 'Chr8'
strand="-"
with open("list") as f:
    for line in f:
        elements = line.split()
        p = subprocess.Popen(['/home/song/zsdp/zsdp', 'lift', '-r', '/netscratch/dep_tsiantis/grp_gan/song/tsianst/GWAS/prepareForGenotyping/chi_v1.fa', '-v', elements[1], '-c', chr, '-p', str(refStart)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print (p)
        out, err = p.communicate()
        start = out.split()[1]

        p = subprocess.Popen(['/home/song/zsdp/zsdp', 'lift', '-r', '/netscratch/dep_tsiantis/grp_gan/song/tsianst/GWAS/prepareForGenotyping/chi_v1.fa', '-v', elements[1], '-c', chr, '-p', str(refEnd)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
        end = out.split()[1]
        chromosome_names, fastas=readFastaFile(elements[2])
        seq = getSubSequence(fastas, chr, int(start), int(end), strand)
        fileoutput.write(">")
        fileoutput.write(elements[0])
        fileoutput.write(" ")
        fileoutput.write(str(start))
        fileoutput.write(" ")
        fileoutput.write(str(end))
        fileoutput.write("\n")
        fileoutput.write(seq)
        fileoutput.write("\n")
fileoutput.close()
