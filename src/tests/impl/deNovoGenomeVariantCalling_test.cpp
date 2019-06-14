//
// Created by song on 9/19/18.
//


#include "include/gtest/gtest.h"
#include "../../impl/impl.h"
#include <string>
#include <iostream>
#include <sstream>
#include "../../impl/deNovoGenomeVariantCalling.h"

TEST(deNovoGenomeVariantCalling, c1){
    std::string parameterFile = "/home/bs674/Dropbox/gean/configure";
    std::map<std::string, std::string> parameters = initialize_paramters(parameterFile, "/home/bs674/Dropbox/gean/");
    std::string refGffFilePath = "/media/bs674/1_8t/protein_ml/b73.gff";
    std::string referenceGenomeSequence = "/media/bs674/1_8t/protein_ml/b73.fa";
    std::string targetGffFilePath = "/media/bs674/1_8t/protein_ml/Zm-CML247-DRAFT-PANZEA-1.0/cml247.gff";
    std::string targetGenomeSequence = "/media/bs674/1_8t/protein_ml/Zm-CML247-DRAFT-PANZEA-1.0/cml247.fa";
    std::string outPutFilePath = "/home/bs674/temp_out";

    int minIntron = 5;
    int minGene = 9;
    size_t widownWidth = 60;
    deNovoGenomeVariantCalling(refGffFilePath, referenceGenomeSequence, targetGffFilePath, targetGenomeSequence,
                               minIntron, minGene, parameters, widownWidth, outPutFilePath);
    ASSERT_EQ(0, 0);
}
