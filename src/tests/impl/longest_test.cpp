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

    std::string referenceGffFile = "/Users/bs674/Zea_mays.AGPv3.31.gff3";
    std::string queryNewGffFile = "/Users/bs674/M017.gff";
    std::string queryGenomeFile = "/Users/bs674/Mo17.fa";
    std::string outputGffFile="/Users/bs674/M017_longestPath_true.gff";
    int minIntron = 5;
    bool keepTandemDuplication=true;
    double scoreThreshold=12;
    double score = 3.0;
    double penalty = -4.0;
    double syntenicScore=1.0;
    double orfScore=1.5;
    double dropLengthThredshold=0.2;
    bool onlySyntenic = true;
    generateLongestOutput( referenceGffFile,  queryNewGffFile, queryGenomeFile, outputGffFile,  minIntron,
                                score, penalty,  scoreThreshold,  keepTandemDuplication, parameters, syntenicScore,
                                orfScore, dropLengthThredshold, onlySyntenic);
    ASSERT_EQ(0, 0);
}


TEST(longestPath, c2){
//    std::string parameterFile = "/Users/bs674/Dropbox/gean/configure";
//    std::map<std::string, std::string> parameters = initialize_paramters(parameterFile, "/Users/bs674/Dropbox/gean/");

    std::string parameterFile = "/Users/bs674/Dropbox/gean/configure";
    std::map<std::string, std::string> parameters = initialize_paramters(parameterFile, "/Users/bs674/Dropbox/gean/");

    std::string referenceGffFile = "/Users/bs674/Zea_mays.AGPv3.31.gff3";
    std::string queryNewGffFile = "/Users/bs674/M017.gff";
    std::string queryGenomeFile = "/Users/bs674/Mo17.fa";
    std::string outputGffFile="/Users/bs674/M017_longestPath_false.gff";
    int minIntron = 5;
    bool keepTandemDuplication=true;
    double scoreThreshold=12;
    double score = 3.0;
    double penalty = -4.0;
    double syntenicScore=1.0;
    double orfScore=1.5;
    double dropLengthThredshold=0.2;
    bool onlySyntenic = false;
    generateLongestOutput( referenceGffFile,  queryNewGffFile, queryGenomeFile, outputGffFile,  minIntron,
                           score, penalty,  scoreThreshold,  keepTandemDuplication, parameters, syntenicScore,
                           orfScore, dropLengthThredshold, onlySyntenic);
    ASSERT_EQ(0, 0);
}


TEST(dagChainer, c1){
//    std::string parameterFile = "/Users/bs674/Dropbox/gean/configure";
//    std::map<std::string, std::string> parameters = initialize_paramters(parameterFile, "/Users/bs674/Dropbox/gean/");

    std::string parameterFile = "/Users/bs674/Dropbox/gean/configure";
    std::map<std::string, std::string> parameters = initialize_paramters(parameterFile, "/Users/bs674/Dropbox/gean/");

    std::string referenceGffFile = "/Users/bs674/Zea_mays.AGPv3.31_chr1.gff";
    std::string queryNewGffFile = "/Users/bs674/M017_chr1.gff";
    std::string queryGenomeFile = "/Users/bs674/Mo17.fa";
    std::string outputGffFile="/Users/bs674/M017_dagChainer_true.gff";
    int minIntron = 5;
    bool keepTandemDuplication=true;
    double syntenicScore=1.0;
    double orfScore=1.5;
    double dropLengthThredshold=0.2;
    int MAX_DIST_BETWEEN_MATCHES=200;
    int BP_GAP_SIZE=1;
    double INDEL_SCORE=-0.05;
    double GAP_OPEN_PENALTY=-0.1;
    double MIN_ALIGNMENT_SCORE = 6.0;
    bool onlySyntenic = true;
    generateDagChainerOutput( referenceGffFile, queryNewGffFile, queryGenomeFile, outputGffFile, minIntron, keepTandemDuplication,
                                   parameters, syntenicScore, orfScore, dropLengthThredshold, MAX_DIST_BETWEEN_MATCHES /*max gap in the term of number of genes*/,
                                   BP_GAP_SIZE, INDEL_SCORE, GAP_OPEN_PENALTY, MIN_ALIGNMENT_SCORE, onlySyntenic, false);
    ASSERT_EQ(0, 0);
}


TEST(dagChainer, c2){
//    std::string parameterFile = "/Users/bs674/Dropbox/gean/configure";
//    std::map<std::string, std::string> parameters = initialize_paramters(parameterFile, "/Users/bs674/Dropbox/gean/");

    std::string parameterFile = "/Users/bs674/Dropbox/gean/configure";
    std::map<std::string, std::string> parameters = initialize_paramters(parameterFile, "/Users/bs674/Dropbox/gean/");

    std::string referenceGffFile = "/Users/bs674/Zea_mays.AGPv3.31.gff3";
    std::string queryNewGffFile = "/Users/bs674/M017.gff";
    std::string queryGenomeFile = "/Users/bs674/Mo17.fa";
    std::string outputGffFile="/Users/bs674/M017_dagChainer_false.gff";
    int minIntron = 5;
    bool keepTandemDuplication=true;
    double syntenicScore=1.0;
    double orfScore=1.5;
    double dropLengthThredshold=0.2;
    int MAX_DIST_BETWEEN_MATCHES=20;
    int BP_GAP_SIZE=1;
    double INDEL_SCORE=-0.05;
    double GAP_OPEN_PENALTY=-0.1;
    double MIN_ALIGNMENT_SCORE = 6.0;
    bool onlySyntenic = false;
    generateDagChainerOutput( referenceGffFile, queryNewGffFile, queryGenomeFile, outputGffFile, minIntron, keepTandemDuplication,
                              parameters, syntenicScore, orfScore, dropLengthThredshold, MAX_DIST_BETWEEN_MATCHES /*max gap in the term of number of genes*/,
                              BP_GAP_SIZE, INDEL_SCORE, GAP_OPEN_PENALTY, MIN_ALIGNMENT_SCORE, onlySyntenic, false);
    ASSERT_EQ(0, 0);
}



TEST(longestQuotaOutput, c1){
    std::string parameterFile = "/Users/bs674/Dropbox/gean/configure";
    std::map<std::string, std::string> parameters = initialize_paramters(parameterFile, "/Users/bs674/Dropbox/gean/");
    std::string referenceGffFile = "/Users/bs674/Zea_mays.AGPv3.31.gff3";
    std::string queryNewGffFile = "/Users/bs674/M017_chr1.gff";
    std::string queryGenomeFile = "/Users/bs674/Mo17.fa";
    std::string outputGffFile="/Users/bs674/M017_quota_longestpath_quota_true.gff";
    int minIntron = 5;
    bool keepTandemDuplication=true;
    double syntenicScore=1.0;
    double orfScore=1.5;
    double dropLengthThredshold=0.2;
    int MAX_DIST_BETWEEN_MATCHES=20;
    double INDEL_SCORE=-0.05;
    double GAP_OPEN_PENALTY=-0.1;
    double MIN_ALIGNMENT_SCORE = 6.0;
    int refMaximumTimes=2;  // if there is a query duplication, then the reference gene could appear twice, so this value should be set as 2
    int queryMaximumTimes=2;
    bool onlySyntenic = true;
    generateLongestQuotaOutput( referenceGffFile, queryNewGffFile, queryGenomeFile, outputGffFile, minIntron, keepTandemDuplication,
        parameters, syntenicScore, orfScore, dropLengthThredshold, INDEL_SCORE, GAP_OPEN_PENALTY, MIN_ALIGNMENT_SCORE, MAX_DIST_BETWEEN_MATCHES,
        refMaximumTimes, queryMaximumTimes, onlySyntenic , false);
    ASSERT_EQ(0, 0);
}



TEST(longestQuotaOutput, c2){
    std::string parameterFile = "/Users/bs674/Dropbox/gean/configure";
    std::map<std::string, std::string> parameters = initialize_paramters(parameterFile, "/Users/bs674/Dropbox/gean/");
    std::string referenceGffFile = "/Users/bs674/Zea_mays.AGPv3.31.gff3";
    std::string queryNewGffFile = "/Users/bs674/M017_chr1.gff";
    std::string queryGenomeFile = "/Users/bs674/Mo17.fa";
    std::string outputGffFile="/Users/bs674/M017_quota_longestpath_quota_false.gff";
    int minIntron = 5;
    bool keepTandemDuplication=true;
    double syntenicScore=1.0;
    double orfScore=1.5;
    double dropLengthThredshold=0.2;
    int MAX_DIST_BETWEEN_MATCHES=20;
    double INDEL_SCORE=-0.05;
    double GAP_OPEN_PENALTY=-0.1;
    double MIN_ALIGNMENT_SCORE = 6.0;
    int refMaximumTimes=2;  // if there is a query duplication, then the reference gene could appear twice, so this value should be set as 2
    int queryMaximumTimes=2;
    bool onlySyntenic = false;
    generateLongestQuotaOutput( referenceGffFile, queryNewGffFile, queryGenomeFile, outputGffFile, minIntron, keepTandemDuplication,
                                parameters, syntenicScore, orfScore, dropLengthThredshold, INDEL_SCORE, GAP_OPEN_PENALTY, MIN_ALIGNMENT_SCORE, MAX_DIST_BETWEEN_MATCHES,
                                refMaximumTimes, queryMaximumTimes, onlySyntenic, false );
    ASSERT_EQ(0, 0);
}
