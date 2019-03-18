//
// Created by song on 9/19/18.
//


#include "include/gtest/gtest.h"
#include "../../impl/impl.h"
#include <string>
#include <iostream>
#include <sstream>
#include "../../impl/deNovoGenomeVariantCalling.h"

TEST(longestPath, c1){
//    std::string parameterFile = "/Users/bs674/Dropbox/gean/configure";
//    std::map<std::string, std::string> parameters = initialize_paramters(parameterFile, "/Users/bs674/Dropbox/gean/");

    std::string parameterFile = "/Users/bs674/Dropbox/gean/configure";
    std::map<std::string, std::string> parameters = initialize_paramters(parameterFile, "/Users/bs674/Dropbox/gean/");

//    std::cout<< std::endl;
//    std::cout<< std::endl;
//    std::cout << "haha" << std::endl;
//    for( std::map<std::string, std::string>::iterator it=parameters.begin();
//        it!=parameters.end(); ++it){
//        std::cout << it->first << "\t" << it->second << std::endl;
//    }

    std::string referenceGffFile = "/Users/bs674/Zea_mays.AGPv3.31.gff3";
    std::string queryNewGffFile = "/Users/bs674/M017.gff";
    std::string queryGenomeFile = "/Users/bs674/Mo17.fa";
    std::string outputGffFile="/Users/bs674/M017_GEAN.gff";
    int minIntron = 5;
    bool keepTandemDuplication=true;
    double scoreThreshold=12;
    double score = 3.0;
    double penalty = -4.0;
    double syntenicScore=1.0;
    double orfScore=1.5;
    double dropLengthThredshold=0.2;

    generateLongestOutput( referenceGffFile,  queryNewGffFile, queryGenomeFile, outputGffFile,  minIntron,
                                score, penalty,  scoreThreshold,  keepTandemDuplication, parameters, syntenicScore,
                                orfScore, dropLengthThredshold);
    ASSERT_EQ(0, 0);
}
