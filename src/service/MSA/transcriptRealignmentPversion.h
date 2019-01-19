//
// Created by baoxing on 10/10/17.
//

#ifndef ANNOTATIONLIFTOVER_TRANSCRIPTREALIGNMENTPVERSION_H
#define ANNOTATIONLIFTOVER_TRANSCRIPTREALIGNMENTPVERSION_H
#include "../../myImportandFunction/myImportantFunction.h"
#include "../reannotation/transcriptRealignmentAndExonerate.h"
#include "../../impl/impl.h"
#include "../../util/util.h"
#include <atomic>
#include <thread>
#include <mutex>

void transcriptRealignmentPversion( Transcript& tartgetTranscript, Transcript& referenceTranscript,
                                    NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix,
                                    std::map<std::string, Fasta>& targetGenome, std::atomic_int & number_of_runing_threads,
                                    std::map<std::string, Fasta>& referenceGenome, std::string chromosomeName, std::map<std::string,
                                    Transcript>& targetTranscriptsHashMap, int & lengthThread,
                                    std::map<std::string, std::string>& parameters, int & minIntron );


#endif //ANNOTATIONLIFTOVER_TRANSCRIPTREALIGNMENTPVERSION_H
