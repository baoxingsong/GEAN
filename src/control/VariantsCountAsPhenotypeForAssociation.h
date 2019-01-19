//
// Created by baoxing on 10/10/17.
//

#ifndef ANNOTATIONLIFTOVER_VARIANTSCOUNTASPHENOTYPEFORASSOCIATION_H
#define ANNOTATIONLIFTOVER_VARIANTSCOUNTASPHENOTYPEFORASSOCIATION_H
#include "../impl/impl.h"
#include "../model/model.h"
#include "../util/util.h"

void countNumberOfTwoneighborSNP( std::string& sdiFile, std::string & outputPrefix, int & rangeLength, std::string& vcfFix, std::map<std::string, std::string>& parameters, const std::string& referenceGenomeFastaFile);
void countNumberSNPAndIndel( std::string& sdiFile, std::string & outputPrefix, int & rangeLength, std::string& vcfFix, std::map<std::string, std::string>& parameters, const std::string& referenceGenomeFastaFile);
void generateRandomSdi( std::string& sdiFile, std::string & outputPrefix, std::string& vcfFix, const std::string& referenceGenomeFastaFile  );
#endif //ANNOTATIONLIFTOVER_VARIANTSCOUNTASPHENOTYPEFORASSOCIATION_H
