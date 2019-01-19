//
// Created by baoxing on 10/10/17.
//

#ifndef ANNOTATIONLIFTOVER_TRANSCRIPTREALIGNMENTANDEXONERATE_H
#define ANNOTATIONLIFTOVER_TRANSCRIPTREALIGNMENTANDEXONERATE_H

#include "../../impl/impl.h"
#include "../../util/util.h"
#include "runExonerate.h"
#include <thread>
#include <atomic>
#include <mutex>
#include "../../myImportandFunction/myImportantFunction.h"

bool transcriptRealignment( Transcript& targetTranscript, int& startTarget, int & endTarget, Transcript& referenceTranscript,
                            NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix,
                            std::map<std::string, Fasta>& targetGenome, std::string& refGenomeSequence,
                            std::map<std::string, Fasta>& referenceGenome, std::string& chromosomeName,
                            std::map<std::string, Transcript>& targetTranscriptsHashMap, int & lengthThread,
                            std::map<std::string, std::string>& parameters, int & minIntron );


void transcriptRealignmentAndExonerate( Transcript tartgetTranscript, Transcript referenceTranscript,
                                        NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix,
                                        std::map<std::string, Fasta>& targetGenome,
                                        std::map<std::string, Fasta>& referenceGenome, std::string chromosomeName, std::map<std::string, Transcript>& targetTranscriptsHashMap,
                                        std::atomic_int & number_of_runing_threads , std::string & prefixUuid, int & lengthThread,
                                        std::map<std::string, std::string>& parameters, int& minIntron );

#endif //ANNOTATIONLIFTOVER_TRANSCRIPTREALIGNMENTANDEXONERATE_H
