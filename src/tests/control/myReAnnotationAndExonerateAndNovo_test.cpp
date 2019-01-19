//
// Created by song on 7/31/18.
//


#include "../../../googletest/googletest/include/gtest/gtest.h"
#include "../../control/myControl.h"
#include <string>
#include "../../myImportandFunction/alignNeedlemanWunsch_simd.h"
#include "../../myImportandFunction/alignNeedlemanWunsch.h"

#include <iostream>
#include <ctime>
#include <ratio>
#include <chrono>

TEST(myReAnnotationAndExonerateAndNovo, c1){
    std::string referenceGenomeSequence = "/biodata/dep_tsiantis/grp_gan/song/flygenomes/ref.fa";
    std::string inputGffFile = "/biodata/dep_tsiantis/grp_gan/song/flygenomes/modified_dmel-all-r6.11_v2.gtf";
    std::string novoGffFilePath = "/biodata/dep_tsiantis/grp_gan/song/flygenomes/augu.gff";
    std::string variants = "/biodata/dep_tsiantis/grp_gan/song/flygenomes/PD1125.sdi";
    std::string outputGffFile = "/netscratch/dep_tsiantis/grp_gan/song/geneStructureAnnotation/fruitFly/newGffFile/PD1125.gff";
    int threads = 4;
    int lengthThreadhold=5000;
    int minIntron = 5;
    bool remove_reference_orf_shift = true;
    std::string vcfFix = "";
    std::map<std::string, std::string> parameters = initialize_paramters("/netscratch/dep_tsiantis/grp_gan/song/geneStructureAnnotation/fruitFly/newGffFile/configure", "/netscratch/dep_tsiantis/grp_gan/song/geneStructureAnnotation/fruitFly/newGffFile/");

    myReAnnotationAndExonerateAndNovo( referenceGenomeSequence, inputGffFile, novoGffFilePath,
                                       variants, outputGffFile, threads, lengthThreadhold, vcfFix, parameters, minIntron, remove_reference_orf_shift);
}
