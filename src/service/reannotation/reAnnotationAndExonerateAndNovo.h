//
// Created by baoxing on 10/10/17.
//

#ifndef ANNOTATIONLIFTOVER_REANNOTATIONANDEXONERATEANDNOVO_H
#define ANNOTATIONLIFTOVER_REANNOTATIONANDEXONERATEANDNOVO_H

#include "../../model/model.h"
#include "../../util/util.h"
#include "../../impl/impl.h"
#include "TranscriptsTogenes.h"
#include "runExonerate.h"
#include <atomic>
#include <thread>
#include "transcriptRealignmentAndExonerate.h"

void reAnnotationAndExonerateAndNovo( std::string& referenceGenomeFilePath, std::string referenceGffFilePath, std::string novoGffFilePath, std::string sdiFile,
                                      std::map<std::string, std::vector<Gene> >& genes, std::map<std::string, Transcript > & targetTranscriptsHashMap, int & maxThread,
                                      std::string & outputGffFile, int & lengthThread, std::string & vcfFix,
                                      std::map<std::string, std::string>& parameters, int& minIntron, bool & remove_reference_orf_shift);
#endif //ANNOTATIONLIFTOVER_REANNOTATIONANDEXONERATEANDNOVO_H
