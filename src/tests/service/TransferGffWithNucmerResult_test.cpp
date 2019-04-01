//
// Created by song on 8/5/18.
//

#include "../../../googletest/googletest/include/gtest/gtest.h"
#include "../../util/util.h"
#include <string>
#include <iostream>
#include <sstream>
#include "../../service/TransferGffWithNucmerResult.h"
TEST(TransferGffWithNucmerResult, c1){
    std::string parameterFile = "/home/song/Dropbox/zsdp/configure";
    std::map<std::string, std::string> parameters = initialize_paramters(parameterFile, "/home/song/Dropbox/zsdp/");
    std::string gffFilePath = "/biodata/dep_tsiantis/grp_gan/song/gennepredication/TAIR10annotationclassification/unmodifyed/TAIR10_GFF3_genes.gff";
    std::string databaseFastaFilePath = "/netscratch/dep_tsiantis/grp_gan/song/wholeGenomeAlignment/tair10.fa";
    std::string queryFastaFilePath="/netscratch/dep_tsiantis/grp_gan/song/wholeGenomeAlignment/ler.fa";
    std::string nucmerFilePath="/netscratch/dep_tsiantis/grp_gan/song/wholeGenomeAlignment/dtql.tab";
    std::string outPutFilePath="/netscratch/dep_tsiantis/grp_gan/song/wholeGenomeAlignment/ler.gff";
    std::cout << std::endl;
    size_t maxLengthForStructureAlignment=10000;
    TransferGffWithNucmerResult( gffFilePath, databaseFastaFilePath, queryFastaFilePath, nucmerFilePath, parameters, outPutFilePath, maxLengthForStructureAlignment);
    ASSERT_EQ(0, 0);
}

TEST(TransferAllExonWithNucmerResult, c1){
    std::string parameterFile = "/home/song/Dropbox/zsdp/configure";
    std::map<std::string, std::string> parameters = initialize_paramters(parameterFile, "/home/song/Dropbox/zsdp/");
    std::string gffFilePath = "/biodata/dep_tsiantis/grp_gan/song/gennepredication/TAIR10annotationclassification/unmodifyed/TAIR10_GFF3_genes.gff";
    std::string databaseFastaFilePath = "/netscratch/dep_tsiantis/grp_gan/song/wholeGenomeAlignment/tair10.fa";
    std::string queryFastaFilePath="/netscratch/dep_tsiantis/grp_gan/song/wholeGenomeAlignment/ler.fa";
    std::string nucmerFilePath="/netscratch/dep_tsiantis/grp_gan/song/wholeGenomeAlignment/test.tab";
    std::string outPutFilePath="/netscratch/dep_tsiantis/grp_gan/song/wholeGenomeAlignment/ler.gff";
    std::cout << std::endl;
    size_t maxLengthForStructureAlignment=10000;
    TransferAllExonWithNucmerResult( gffFilePath, databaseFastaFilePath, queryFastaFilePath, nucmerFilePath, parameters, outPutFilePath, 5 , false, 60, maxLengthForStructureAlignment);
    ASSERT_EQ(0, 0);
}

TEST(TransferAllExonWithNucmerResult, c1_1){
    std::string parameterFile = "/home/song/Dropbox/zsdp/configure";
    std::map<std::string, std::string> parameters = initialize_paramters(parameterFile, "/home/song/Dropbox/zsdp/");
    std::string gffFilePath = "/biodata/dep_tsiantis/grp_gan/song/gennepredication/TAIR10annotationclassification/unmodifyed/TAIR10_GFF3_genes.gff";
    std::string databaseFastaFilePath = "/netscratch/dep_tsiantis/grp_gan/song/wholeGenomeAlignment/tair10.fa";
    std::string queryFastaFilePath="/netscratch/dep_tsiantis/grp_gan/song/wholeGenomeAlignment/ler.fa";
    std::string nucmerFilePath="/netscratch/dep_tsiantis/grp_gan/song/wholeGenomeAlignment/test.tab";
    std::string outPutFilePath="/netscratch/dep_tsiantis/grp_gan/song/wholeGenomeAlignment/ler.gff";
    std::cout << std::endl;
    size_t maxLengthForStructureAlignment=10000;
    TransferAllExonWithNucmerResult( gffFilePath, databaseFastaFilePath, queryFastaFilePath, nucmerFilePath, parameters, outPutFilePath, 5 , true, 60, maxLengthForStructureAlignment);
    ASSERT_EQ(0, 0);
}


TEST(TransferAllExonWithNucmerResult, c2){
    std::string parameterFile = "/home/song/Dropbox/zsdp/configure";
    std::map<std::string, std::string> parameters = initialize_paramters(parameterFile, "/home/song/Dropbox/zsdp/");
    std::string gffFilePath = "/netscratch/dep_tsiantis/grp_gan/song/zbl/haiwang/Zea_mays.AGPv4.38.gtf";
    std::string databaseFastaFilePath = "/netscratch/dep_tsiantis/grp_gan/song/zbl/haiwang/Zea_mays.AGPv4.dna.toplevel.fa";
    std::string queryFastaFilePath="/netscratch/dep_tsiantis/grp_gan/song/zbl/haiwang/Mo17.fasta";
    std::string nucmerFilePath="/netscratch/dep_tsiantis/grp_gan/song/zbl/haiwang/dtql.tab";
    std::string outPutFilePath="/netscratch/dep_tsiantis/grp_gan/song/zbl/haiwang/Mo17.gff";
    std::cout << std::endl;
    size_t maxLengthForStructureAlignment=10000;
    TransferAllExonWithNucmerResult( gffFilePath, databaseFastaFilePath, queryFastaFilePath, nucmerFilePath, parameters, outPutFilePath, 5, false, 60, maxLengthForStructureAlignment);
    ASSERT_EQ(0, 0);
}

TEST(TransferAllExonWithNucmerResult, c3){
    std::string parameterFile = "/home/song/Dropbox/zsdp/configure";
    std::map<std::string, std::string> parameters = initialize_paramters(parameterFile, "/home/song/Dropbox/zsdp/");
    std::string gffFilePath = "/netscratch/dep_tsiantis/grp_gan/song/wholeGenomeAlignment/TAIR10_GFF3_genes.gff";
    std::string databaseFastaFilePath = "/netscratch/dep_tsiantis/grp_gan/song/wholeGenomeAlignment/tair10.fa";
    std::string queryFastaFilePath="/netscratch/dep_tsiantis/grp_gan/song/wholeGenomeAlignment/chi_v1.fa";
    std::string nucmerFilePath="/netscratch/dep_tsiantis/grp_gan/song/wholeGenomeAlignment/test.tab";
    std::string outPutFilePath="/netscratch/dep_tsiantis/grp_gan/song/wholeGenomeAlignment/test_chi.gff";
    std::cout << std::endl;
    size_t maxLengthForStructureAlignment=10000;
    TransferAllExonWithNucmerResult( gffFilePath, databaseFastaFilePath, queryFastaFilePath, nucmerFilePath, parameters, outPutFilePath, 5, false, 60, maxLengthForStructureAlignment);
    ASSERT_EQ(0, 0);
}

TEST(TransferAllExonWithNucmerResult, c4){
    std::string parameterFile = "/home/song/Dropbox/zsdp/configure";
    std::map<std::string, std::string> parameters = initialize_paramters(parameterFile, "/home/song/Dropbox/zsdp/");
    std::string gffFilePath = "/netscratch/dep_tsiantis/grp_gan/song/wholeGenomeAlignment/TAIR10_GFF3_genes.gff";
    std::string databaseFastaFilePath = "/netscratch/dep_tsiantis/grp_gan/song/wholeGenomeAlignment/tair10.fa";
    std::string queryFastaFilePath="/netscratch/dep_tsiantis/grp_gan/song/wholeGenomeAlignment/chi_v1.fa";
    std::string nucmerFilePath="/netscratch/dep_tsiantis/grp_gan/song/wholeGenomeAlignment/chi.tab";
    std::string outPutFilePath="/netscratch/dep_tsiantis/grp_gan/song/wholeGenomeAlignment/chi.gff";
    std::cout << std::endl;
    size_t maxLengthForStructureAlignment=10000;
    TransferAllExonWithNucmerResult( gffFilePath, databaseFastaFilePath, queryFastaFilePath, nucmerFilePath, parameters, outPutFilePath, 5, false, 60, maxLengthForStructureAlignment);
    ASSERT_EQ(0, 0);
}

TEST(TransferAllExonWithNucmerResult, c5){
    std::string parameterFile = "/home/bs674/Dropbox/gean/configure";
    std::map<std::string, std::string> parameters = initialize_paramters(parameterFile, "/home/bs674/Dropbox/gean/");
    std::string gffFilePath = "/media/bs674/genomeSequence/maize/Zea_mays.AGPv4.34.gff3";
    std::string databaseFastaFilePath = "/media/bs674/genomeSequence/maize/Zea_mays.AGPv4.dna.toplevel.fa";
    std::string queryFastaFilePath="/media/bs674/A1013-0002_annotation/gapclosed.fasta";
    std::string nucmerFilePath="/media/bs674/A1013-0002_annotation/Zm00001d020160.sam";
    std::string outPutFilePath="/media/bs674/A1013-0002_annotation/Zm00001d020160.gff";
    std::cout << std::endl;
    size_t maxLengthForStructureAlignment=10000;
    TransferAllExonWithSpliceAlignmentResult( gffFilePath, databaseFastaFilePath, queryFastaFilePath, nucmerFilePath, parameters, outPutFilePath, 5, 6000, maxLengthForStructureAlignment, 1);

//    TransferAllExonWithNucmerResult( gffFilePath, databaseFastaFilePath, queryFastaFilePath, nucmerFilePath, parameters, outPutFilePath, 5, false, 60, maxLengthForStructureAlignment);
    ASSERT_EQ(0, 0);
}


