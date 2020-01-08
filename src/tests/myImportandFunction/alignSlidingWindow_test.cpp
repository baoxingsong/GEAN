//
// Created by baoxing on 2/25/18.
//

#include "../../../googletest/googletest/include/gtest/gtest.h"
#include "../../impl/impl.h"
#include <string>
#include "../../myImportandFunction/alignSlidingWindow.h"

#include <iostream>
#include <ctime>
#include <ratio>
#include <chrono>


TEST(alignSlidingWindow, c1){
    std::cout << std::endl;
    std::string parameterFile = "configure";
    std::map<std::string, std::string> parameters = initialize_paramters(parameterFile, "./");
    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix(parameters);

    std::map<std::string, Fasta> querySequences;
    std::map<std::string, Fasta> databaseSequences;
    //readFastaFile("/home/song/wholeGenomeAlignment/ler.fa", querySequences);
    //readFastaFile("/home/song/wholeGenomeAlignment/tair10.fa", databaseSequences);

    readFastaFile("/home/bs674/genomes/ler.fa", querySequences);
    readFastaFile("/home/bs674/genomes/tair10.fa", databaseSequences);

    std::string dna_q = getSubsequence(querySequences, "000005F|arrow", 1, 263360);
    std::string dna_d = getSubsequence(databaseSequences, "Chr1", 7132, 270470);

    std::string _alignment_q;
    std::string _alignment_d;
    std::cout << "begin to align" << std::endl;
    alignSlidingWindow(dna_q, dna_d, _alignment_q, _alignment_d, 50, parameters, nucleotideCodeSubstitutionMatrix);

    size_t linewidth = 100;
    std::size_t pos = 0;
    while (pos < _alignment_q.length()) {
        std::cout << _alignment_q.substr(pos, linewidth) << std::endl;
        std::cout << _alignment_d.substr(pos, linewidth) << std::endl << std::endl;
        pos += linewidth;
    }
    ASSERT_EQ(0, 0);
}



TEST(alignSlidingWindow, c2){ // it takes about 25mins to finish the alignment
    std::cout << std::endl;
    std::string parameterFile = "configure";
    std::map<std::string, std::string> parameters = initialize_paramters(parameterFile, "./");
    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix(parameters);

    std::map<std::string, Fasta> querySequences;

    std::map<std::string, Fasta> databaseSequences;
    //readFastaFile("/home/song/wholeGenomeAlignment/ler.fa", querySequences);
    //readFastaFile("/home/song/wholeGenomeAlignment/tair10.fa", databaseSequences);

    readFastaFile("/home/bs674/genomes/ler.fa", querySequences);
    readFastaFile("/home/bs674/genomes/tair10.fa", databaseSequences);

    std::map<std::string, std::vector<AlignmentMatch>> alignmentMatchsMap;
    nucmerRead("/home/bs674/align.tbl", alignmentMatchsMap);
    std::cout << "begin to align" << std::endl;

    std::ofstream ofile;
    ofile.open( "/home/bs674/chr1_aln");
    for(std::map<std::string, std::vector<AlignmentMatch>>::iterator it1=alignmentMatchsMap.begin();
        it1!=alignmentMatchsMap.end(); ++it1) {
        if( databaseSequences.find(it1->first) != databaseSequences.end() ) {
            for (AlignmentMatch alignmentMatch : it1->second) {
                std::string dna_q = getSubsequence(querySequences, alignmentMatch.getQueryChr(), alignmentMatch.getQueryStart(), alignmentMatch.getQueryEnd(), alignmentMatch.getQueryStrand());
                std::string dna_d = getSubsequence(databaseSequences, alignmentMatch.getDatabaseChr(), alignmentMatch.getDatabaseStart(), alignmentMatch.getDatabaseEnd());

                std::string _alignment_q;
                std::string _alignment_d;
                alignSlidingWindow(dna_q, dna_d, _alignment_q, _alignment_d, 1000, parameters, nucleotideCodeSubstitutionMatrix);

                size_t linewidth = 100;
                std::size_t pos = 0;
                while (pos < _alignment_q.length()) {
                    ofile << _alignment_q.substr(pos, linewidth) << std::endl;
                    ofile << _alignment_d.substr(pos, linewidth) << std::endl << std::endl;
                    pos += linewidth;
                }
                ofile << std::endl << std::endl;
            }
        }
    }
    ofile.close();
    ASSERT_EQ(0, 0);
}




TEST(bandAlign, c1){ // it takes about 25mins to finish the alignment
    std::cout << std::endl;
    std::string parameterFile = "configure";
    std::map<std::string, std::string> parameters = initialize_paramters(parameterFile, "./");
    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix(parameters);

    std::map<std::string, Fasta> querySequences;

    std::map<std::string, Fasta> databaseSequences;
    //readFastaFile("/home/song/wholeGenomeAlignment/ler.fa", querySequences);
    //readFastaFile("/home/song/wholeGenomeAlignment/tair10.fa", databaseSequences);

    readFastaFile("/home/bs674/genomes/ler.fa", querySequences);
    readFastaFile("/home/bs674/genomes/tair10.fa", databaseSequences);

    std::map<std::string, std::vector<AlignmentMatch>> alignmentMatchsMap;
    nucmerRead("/home/bs674/align.tbl", alignmentMatchsMap);
    std::cout << "begin to align" << std::endl;

    std::ofstream ofile;
    ofile.open( "/home/bs674/chr1_aln_banded");
    for(std::map<std::string, std::vector<AlignmentMatch>>::iterator it1=alignmentMatchsMap.begin();
        it1!=alignmentMatchsMap.end(); ++it1) {
        if( databaseSequences.find(it1->first) != databaseSequences.end() ) {
            for (AlignmentMatch alignmentMatch : it1->second) {
                std::string dna_q = getSubsequence(querySequences, alignmentMatch.getQueryChr(), alignmentMatch.getQueryStart(), alignmentMatch.getQueryEnd(), alignmentMatch.getQueryStrand());
                std::string dna_d = getSubsequence(databaseSequences, alignmentMatch.getDatabaseChr(), alignmentMatch.getDatabaseStart(), alignmentMatch.getDatabaseEnd());

                std::string _alignment_q;
                std::string _alignment_d;
                bandAlign(dna_q, dna_d, _alignment_q, _alignment_d, 1000, parameters, nucleotideCodeSubstitutionMatrix);


                size_t linewidth = 100;
                std::size_t pos = 0;
                while (pos < _alignment_q.length()) {
                    ofile << _alignment_q.substr(pos, linewidth) << std::endl;
                    ofile << _alignment_d.substr(pos, linewidth) << std::endl << std::endl;
                    pos += linewidth;
                }
                ofile << std::endl << std::endl;
            }
        }
    }
    ofile.close();
    ASSERT_EQ(0, 0);
}




TEST(alignSlidingWindow, c3){
    std::cout << std::endl;
    std::string parameterFile = "configure";
    std::map<std::string, std::string> parameters = initialize_paramters(parameterFile, ".");
    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix(parameters);
    std::string dna_q = "TAAACCT";
    std::string dna_d = "CCTATCTGA";
    std::string _alignment_q;
    std::string _alignment_d;
    alignSlidingWindow(dna_q, dna_d, _alignment_q, _alignment_d, 5, parameters, nucleotideCodeSubstitutionMatrix);
    std::cout << _alignment_q << std::endl;
    std::cout << _alignment_d << std::endl;
    ASSERT_EQ(0, 0);
}




TEST(alignSlidingWindow, c5){
    std::cout << std::endl;
    std::string parameterFile = "configure";
    std::map<std::string, std::string> parameters = initialize_paramters(parameterFile, "/Users/song/Dropbox/gean");
    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix(parameters);

    std::map<std::string, Fasta> querySequences;
    std::map<std::string, Fasta> databaseSequences;

    readFastaFile("/Users/song/Dropbox/A1033_01_h50.fasta", querySequences);
    readFastaFile("/Users/song/Dropbox/adapter.fa", databaseSequences);


    std::string dna_q;
    std::string dna_d;
    std::cout << "begin to alignment" << std::endl;
    for( std::map<std::string, Fasta>::iterator it=querySequences.begin(); it != querySequences.end(); it++ ){
        for( std::map<std::string, Fasta>::iterator it1=databaseSequences.begin(); it1 != databaseSequences.end(); it1++ ) {
            dna_q = it->second.getSequence();
            dna_d = it1->second.getSequence();
            std::string _alignment_q="";
            std::string _alignment_d="";

            alignSlidingWindow(dna_q, dna_d, _alignment_q, _alignment_d, 60, parameters,
                               nucleotideCodeSubstitutionMatrix);
            std::cout << "sequence" << std::endl;
            std::cout << dna_q << std::endl;
            std::cout << dna_d << std::endl;
            std::cout << "alignment" << std::endl;
            std::cout << _alignment_q << std::endl;
            std::cout << _alignment_d << std::endl << std::endl;
        }
    }
    ASSERT_EQ(0, 0);
}
