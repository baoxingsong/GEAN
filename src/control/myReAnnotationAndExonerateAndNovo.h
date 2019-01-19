//
// Created by baoxing on 10/10/17.
//

#ifndef ANNOTATIONLIFTOVER_MYREANNOTATIONANDEXONERATEANDNOVO_H
#define ANNOTATIONLIFTOVER_MYREANNOTATIONANDEXONERATEANDNOVO_H

#include "../util/util.h"
#include "../impl/impl.h"
#include "../service/service.h"

void myReAnnotationAndExonerateAndNovo( std::string& referenceGenomeFile, std::string& inputGffFile,std::string& novoGffFilePath,
                                        std::string& variantsFile, std::string& outputGffFile, int &maxThread,
                                        int& lengthThread, std::string & vcfFix, std::map<std::string, std::string>& parameters,
                                        int & minIntron, bool & remove_reference_orf_shift);
void hereOutPutLiftOrOrthologousResult(std::map<std::string, std::vector<Gene> >& genes,
        std::map<std::string, Transcript> & targetTranscriptsHashMap, std::string& outputGffFile );

#endif //ANNOTATIONLIFTOVER_MYREANNOTATIONANDEXONERATEANDNOVO_H
