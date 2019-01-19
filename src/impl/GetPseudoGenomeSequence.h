//
// Created by baoxing on 10/10/17.
//

#ifndef ANNOTATIONLIFTOVER_GETPSEUDIGENOMESEQUENCE_H
#define ANNOTATIONLIFTOVER_GETPSEUDIGENOMESEQUENCE_H

#include "../model/model.h"
#include "readSdiFile.h"
#include "WriteFasta.h"
#include "readFastaFile.h"


int getPseudoGenomeSequence(std::map<std::string, Fasta>& referenceSequences, const std::string& sdiFile,
                            std::map<std::string, Fasta>& targetSequences, const std::string& vcfFix, std::map<std::string, std::string>& parameters);

int getPseudoGenomeSequence(const std::string& referenceGenomeFastaFile,
                            const std::string& sdiFile, std::map<std::string, Fasta>& targetSequences, const std::string& vcfFix, std::map<std::string, std::string>& parameters);

int getPseudoGenomeSequence(const std::string& referenceGenomeFastaFile,
                            const std::string& sdiFile, const std::string& outputFile, const std::string& vcfFix, std::map<std::string, std::string>& parameters);

int getPseudoGenomeSequence(std::map<std::string, Fasta>& referenceSequences, std::map<std::string, std::vector<Variant> >& variantsMap,
                            std::map<std::string, Fasta>& targetSequences, std::map<std::string, std::string>& parameters);

int getPseudoGenomeSequence(std::map<std::string, Fasta>& referenceSequences, std::map<std::string, std::vector<Variant> >& variantsMap,
                            std::map<std::string, Fasta>& targetSequences, std::string& chromosome, std::map<std::string, std::string>& parameters);

int getPseudoGenomeSequence(std::map<std::string, Fasta>& referenceSequences, std::map<std::string, std::vector<Variant> >& variantsMap,
                            std::map<std::string, Fasta>& targetSequences, std::set<std::string>& chromosomes, std::map<std::string, std::string>& parameters);

#endif //ANNOTATIONLIFTOVER_GETPSEUDIGENOMESEQUENCE_H
