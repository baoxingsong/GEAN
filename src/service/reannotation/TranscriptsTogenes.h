//
// Created by baoxing on 10/10/17.
//

#ifndef ANNOTATIONLIFTOVER_TRANSCRIPTSTOGENES_H
#define ANNOTATIONLIFTOVER_TRANSCRIPTSTOGENES_H

#include "../../model/model.h"
#include <climits>
void TranscriptsTogenes(std::map<std::string, std::string >& transcript_to_gene_map, std::map<std::string, std::vector<Gene> >& genes,
                        std::map<std::string, Transcript>& targetTranscriptsHashMap);

#endif //ANNOTATIONLIFTOVER_TRANSCRIPTSTOGENES_H
