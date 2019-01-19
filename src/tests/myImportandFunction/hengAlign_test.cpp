//
// Created by baoxing on 2/25/18.
//

#include "../../../googletest/googletest/include/gtest/gtest.h"
#include "../../impl/impl.h"
#include <string>
#include "../../myImportandFunction/hengAlign.h"

#include <iostream>
#include <ctime>
#include <ratio>
#include <chrono>


TEST(hengAlign, c1){
    std::cout << std::endl;
    std::string parameterFile = "configure";
    std::map<std::string, std::string> parameters = initialize_paramters(parameterFile, "./");
    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix(parameters);

    std::map<std::string, Fasta> querySequences;
    std::map<std::string, Fasta> databaseSequences;
    //readFastaFile("/Users/song/tair11/Col.fa", querySequences);
    //readFastaFile("/home/song/wholeGenomeAlignment/tair10.fa", databaseSequences);
/*
    mm_tbuf_t *tbuf = mm_tbuf_init();
    ksw_extz_t ez;
    memset(&ez, 0, sizeof(ksw_extz_t));
*/

    std::string dna_d = "TATATATATATACCTTG";
    std::string dna_q = "CCTTGCGCGCGC";

    std::string _alignment_q;
    std::string _alignment_d;
    std::cout << "begin to align" << std::endl;
    int8_t _open_gap_penalty=stoi(get_parameters("alignmentOpenGapP", parameters));
    int8_t _extend_gap_penalty=stoi(get_parameters("alignmentExtendGapP", parameters));
    std::stack<char> SQ;
    std::stack<char> SD;
    hengAlign(dna_q, dna_d, SQ, SD, _open_gap_penalty,
              _extend_gap_penalty, nucleotideCodeSubstitutionMatrix);
//              _extend_gap_penalty, nucleotideCodeSubstitutionMatrix, &(tbuf->km), &ez);
    //mm_tbuf_destroy(tbuf);
    while (!SQ.empty()) {
        _alignment_q += SQ.top();
        _alignment_d += SD.top();
        SQ.pop();
        SD.pop();
    }
    std::cout << _alignment_q << std::endl;
    std::cout << _alignment_d << std::endl;
    ASSERT_EQ(0, 0);
}
