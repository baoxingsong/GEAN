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
//    std::string gffFilePath = "/Users/song/CM4.0.noZero.gff3";
    std::string gffFilePath = "/Users/song/Zea_mays.AGPv3.31.gff3";
    std::string databaseFastaFilePath = "/Users/song/Zea_mays.AGPv3.31.dna.genome.fa";
    std::string queryFastaFilePath="/Users/song/Zm-CML247-DRAFT-PANZEA-1.0.fasta";
//    std::string databaseFastaFilePath = "/Users/song/CM3.6.1_pseudomol.fa";
//    std::string queryFastaFilePath="/Users/song/CM3.6.1_pseudomol.fa";
    //std::string samFilePath="/Users/song/MELO3C031356.2.1.sam";
//    std::string samFilePath="/Users/song/MELO3C035475.2.1.sam";
    //std::string samFilePath="/Users/song/MELO3C003750.2.1.sam";
    std::string samFilePath="/Users/song/minimap2ab";
    std::string outPutFilePath="/Users/song/ler.gff";
//    std::cout << " go to the function" << std::endl;
    size_t maxLengthForStructureAlignment=60000;
    int minIntron = 5;
    int windowWidth = 20;
    TransferAllExonWithSpliceAlignmentResult( gffFilePath, databaseFastaFilePath, queryFastaFilePath, samFilePath, parameters, outPutFilePath, minIntron, windowWidth, maxLengthForStructureAlignment, 1, 1);
    ASSERT_EQ(0, 0);
}
