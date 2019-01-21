//
// Created by Baoxing song on 2019-01-20.
//

#include "../../../googletest/googletest/include/gtest/gtest.h"
#include "../../util/util.h"
#include <string>
#include <iostream>
#include <sstream>
#include "../../service/TransferGffWithNucmerResult.h"
TEST(TransferAllExonWithSpliceAlignmentResult, c1){
    std::string parameterFile = "/Users/song/Dropbox/GEAN/configure";
    std::map<std::string, std::string> parameters = initialize_paramters(parameterFile, "/Users/song/Dropbox/GEAN/");
    std::string gffFilePath = "/Users/song/CM4.0.noZero.gff3";
    std::string databaseFastaFilePath = "/Users/song/CM3.6.1_pseudomol.fa";
    std::string queryFastaFilePath="/Users/song/CM3.6.1_pseudomol.fa";
    std::string samFilePath="/Users/song/CM361_vs_CM361.sam";
    std::string outPutFilePath="/Users/song/ler.gff";
//    std::cout << " go to the function" << std::endl;
    size_t maxLengthForStructureAlignment=600;
    int minIntron = 5;
    int windowWidth = 60;
    TransferAllExonWithSpliceAlignmentResult( gffFilePath, databaseFastaFilePath, queryFastaFilePath, samFilePath, parameters, outPutFilePath, minIntron, windowWidth, maxLengthForStructureAlignment);
    ASSERT_EQ(0, 0);
}
