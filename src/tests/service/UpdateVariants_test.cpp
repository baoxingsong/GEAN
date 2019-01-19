//
// Created by song on 8/25/18.
//

#include "../../../googletest/googletest/include/gtest/gtest.h"
#include "../../service/UpdateVariants.h"
TEST(UpdateVariants, c1){
    std::string parameterFile = "/home/song/Dropbox/zsdp/configure";
    std::map<std::string, std::string> parameters = initialize_paramters(parameterFile, "/home/song/Dropbox/zsdp/");
    std::string gffFilePath = "/biodata/dep_tsiantis/grp_gan/song/gennepredication/TAIR10annotationclassification/unmodifyed/TAIR10_GFF3_genes.gff";
    std::string databaseFastaFilePath = "/netscratch/dep_tsiantis/grp_gan/song/wholeGenomeAlignment/tair10.fa";
    std::string sdiFile="/biodata/dep_tsiantis/grp_gan/song/rdINDELallHere/inputData/fullSdiFile/PA7213.sdi";
    std::string outPutFilePath="/netscratch/dep_tsiantis/grp_gan/song/wholeGenomeAlignment/ler.sdi";
    std::cout << std::endl;

    updateVariants( gffFilePath, databaseFastaFilePath, sdiFile, "", 5, parameters, outPutFilePath);
    ASSERT_EQ(0, 0);
}

TEST(UpdateVariants, c2){
    std::string parameterFile = "/home/song/Dropbox/zsdp/configure";
    std::map<std::string, std::string> parameters = initialize_paramters(parameterFile, "/home/song/Dropbox/zsdp/");
    std::string gffFilePath = "/biodata/dep_tsiantis/grp_gan/song/gennepredication/TAIR10annotationclassification/unmodifyed/TAIR10_GFF3_genes.gff";
    std::string databaseFastaFilePath = "/netscratch/dep_tsiantis/grp_gan/song/wholeGenomeAlignment/tair10.fa";
    std::string sdiFile="/biodata/dep_tsiantis/grp_gan/song/rdINDELallHere/inputData/fullSdiFile//PA1074.sdi";
    std::string outPutFilePath="/netscratch/dep_tsiantis/grp_gan/song/wholeGenomeAlignment/PA1074.sdi";
    std::cout << std::endl;

    updateVariants( gffFilePath, databaseFastaFilePath, sdiFile, "", 5, parameters, outPutFilePath);
    ASSERT_EQ(0, 0);
}
