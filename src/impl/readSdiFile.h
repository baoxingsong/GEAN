//
// Created by baoxing on 10/10/17.
//

#ifndef ANNOTATIONLIFTOVER_READSDIFILE_H
#define ANNOTATIONLIFTOVER_READSDIFILE_H

#include "../model/model.h"
#include "../util/util.h"

void readSdiFile (const std::string& filePath, std::map<std::string, std::vector<Variant> >& variantsMap, const std::string & vcfFix, std::map<std::string, Fasta>& referenceGenome);
void readSdiFile (const std::string& filePath, std::map<std::string, std::vector<Variant> >& variantsMap, const std::string& chromosome, const std::string & vcfFix, std::map<std::string, Fasta>& referenceGenome);


#endif //ANNOTATIONLIFTOVER_READSDIFILE_H
