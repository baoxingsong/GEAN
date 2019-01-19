//
// Created by baoxing on 10/10/17.
//

#ifndef ANNOTATIONLIFTOVER_RESDIFROMMSA_H
#define ANNOTATIONLIFTOVER_RESDIFROMMSA_H

#include "../../impl/impl.h"
#include <thread>
#include <atomic>
#include <mutex>

//void msaFileReadPversion(MsaFileRecord msaFileRecord, std::string fileLocation, std::set<std::string>& accessionNames, std::vector<MsaFileRecord>& msaFileRecords, std::atomic_int& number_of_runing_threads);

void trimWindow(std::vector<TwoSeqOfMsaResult>& twoSeqOfMsaResults, std::string& accessionName, std::vector<MsaFileRecord>& msaFileRecords,
                std::map<std::string, Fasta>& referenceGenome, std::map<std::string, Fasta>& targetGenome, std::string & chrName,
                int& thisTargetChromosomeLength);

//void newSdiFileForOneAccession(std::string& accessionName, std::map<std::string, std::string>& sdiFiles,
//                               std::string & vcfFix, std::map<std::string, Fasta>& referenceGenome,
//                               std::map<std::string, std::string>& parameters, std::string & chrName,
//                               std::vector<MsaFileRecord>& msaFileRecords, std::atomic_int& number_of_runing_threads);

void constructSdiFromMsa(std::vector<std::string>& chromosomes,  std::string& folder, std::string& outputFolder, std::string & referenceGenomeFilePath,
                         std::map<std::string, std::string>& sdiFiles, std::string & vcfFix, std::map<std::string, std::string>& parameters, int & maxThread);

#endif //ANNOTATIONLIFTOVER_RESDIFROMMSA_H
