//
// Created by song on 9/16/18.
//


#include "../../../googletest/googletest/include/gtest/gtest.h"
#include "../../control/myControl.h"
#include "../../control/myCoordinateLiftOver.h"
#include <string>

#include <iostream>
#include <ctime>
#include <ratio>
#include <chrono>

TEST(myGffCoordinateLiftOver, c1){
    std::string referenceGenomeSequence = "/biodata/dep_tsiantis/grp_gan/song/flygenomes/ref.fa";
    std::string inputGffFile = "/biodata/dep_tsiantis/grp_gan/song/flygenomes/modified_dmel-all-r6.11_v2.gtf";
    std::string variants = "/biodata/dep_tsiantis/grp_gan/song/flygenomes/PD1191.sdi";
    std::string outputGffFile = "/home/song/Desktop/PD1191.sdi";
    std::string prefix="";
    myGffCoordinateLiftOver(variants, inputGffFile, outputGffFile, prefix, referenceGenomeSequence);
}
//void myGffCoordinateLiftOver( std::string& sdiFile, std::string& gffFile, std::string& outputFile, std::string& vcfFix, const std::string& referenceGenomeFastaFile );
