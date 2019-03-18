//
// Created by song on 9/19/18.
//


#include "include/gtest/gtest.h"
#include "../../impl/impl.h"
#include <string>
#include <iostream>
#include <sstream>
#include "../../impl/deNovoGenomeVariantCalling.h"

TEST(orf, c1){
    std::string parameterFile = "/Users/bs674/Dropbox/gean/configure";
    std::map<std::string, std::string> parameters = initialize_paramters(parameterFile, "/Users/bs674/Dropbox/gean/");
    std::string gffFile = "/Users/bs674/M017.gff";
    std::string genomeFile = "/Users/bs674/Mo17.fa";
    std::string outputGffFile="/Users/bs674/M017_ORF.gff";
    int minIntron = 5;
    outPutORFConserveredTranscripts( genomeFile, gffFile, outputGffFile, minIntron, parameters);
    ASSERT_EQ(0, 0);
}
