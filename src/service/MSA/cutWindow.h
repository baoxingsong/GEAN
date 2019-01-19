//
// Created by baoxing on 10/10/17.
//

#ifndef ANNOTATIONLIFTOVER_CUTWINDOW_H
#define ANNOTATIONLIFTOVER_CUTWINDOW_H

#include <thread>
#include "../../util/util.h"
#include "../../impl/impl.h"
#include "transcriptRealignmentPversion.h"
#include <atomic>
#include <mutex>
//
//#include <pthread.h>
//#include <sched.h>
#include <unistd.h>




//void writeToFileParallelVersion_ref(Window& window, int& thisChrLength, int& msaWindowOverlap, std::string& chr, std::map<std::string, Fasta>& sequences, std::string& accessionName, std::string& folder, bool& append, std::atomic_int & number_of_runing_threads);
//
//// this function is used to output the sequence for non-reference accessions
//void writeToFileParallelVersion_non_ref(Window& window, int start, int end, std::string& chr, std::map<std::string, Fasta>& sequences, std::string& accessionName, std::string& folder, bool& append, std::atomic_int & number_of_runing_threads);

void checkOrfPversion(Transcript& transcript, std::map<std::string, Fasta>& genome, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix, int& minIntron, std::atomic_int & number_of_runing_threads);

void prepareForMsa( std::string& referenceGenomeFilePath, std::string& referenceGffFilePath, std::map<std::string, std::string>& sdiFiles,
                    int& maxThread, int & lengthThread, std::string & vcfFix, std::map<std::string, std::string>& parameters,
                    bool& append, int& msaWindowSize, int& msaWindowOverlap, int& minIntron, int & outputPoolSize);

void constructSdiFromMsa_v_beta(std::vector<std::string>& chromosomes, std::set<std::string>& accessionNames, std::string& folder, std::string& outputFolder, std::string & referenceGenomeFilePath,
                                std::map<std::string, std::string>& sdiFiles, std::string & vcfFix, std::map<std::string, std::string>& parameters);



#endif //ANNOTATIONLIFTOVER_CUTWINDOW_H
