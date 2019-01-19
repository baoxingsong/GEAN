//
// Created by baoxing on 10/10/17.
//

#ifndef ANNOTATIONLIFTOVER_MYCOORDINATELIFTOVER_H
#define ANNOTATIONLIFTOVER_MYCOORDINATELIFTOVER_H
#include "../model/model.h"
#include "../util/util.h"
#include "../impl/impl.h"

int myGetChangedFromReference( std::string& sdiFile, std::string& chromosome, int& position, std::string& vcfFix, const std::string& referenceGenomeFastaFile );
int myGetReferenceFromChanged( std::string& sdiFile, std::string& chromosome, int& position, std::string& vcfFix, const std::string& referenceGenomeFastaFile );

void myGffCoordinateLiftOver( std::string& sdiFile, std::string& gffFile, std::string& outputFile, std::string& vcfFix, const std::string& referenceGenomeFastaFile );
void myRevGffCoordinateLiftOver( std::string& sdiFile, std::string& gffFile, std::string& outputFile, std::string& vcfFix, const std::string& referenceGenomeFastaFile );
#endif //ANNOTATIONLIFTOVER_MYCOORDINATELIFTOVER_H
