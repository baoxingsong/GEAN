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

    readFastaFile("/netscratch/dep_tsiantis/grp_gan/song/wholeGenomeAlignment/ler.fa", querySequences);
    readFastaFile("/netscratch/dep_tsiantis/grp_gan/song/wholeGenomeAlignment/tair10.fa", databaseSequences);

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

    readFastaFile("/Users/song/Dropbox/tair11/ler.fa", querySequences);
    readFastaFile("/Users/song/Dropbox/tair11/Col.fa", databaseSequences);

    std::map<std::string, std::vector<AlignmentMatch>> alignmentMatchsMap;
    nucmerRead("/Users/song/Dropbox/tair11/align.tbl", alignmentMatchsMap);
    std::cout << "begin to align" << std::endl;

    std::ofstream ofile;
    ofile.open( "/Users/song/Dropbox/tair11/chr1_aln");
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

TEST(alignSlidingWindow, c3){
    std::cout << std::endl;
    std::string parameterFile = "configure";
    std::map<std::string, std::string> parameters = initialize_paramters(parameterFile, "./");
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
