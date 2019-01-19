//
// Created by baoxing on 10/10/17.
//

#ifndef ANNOTATIONLIFTOVER_RUNEXONERATE_H
#define ANNOTATIONLIFTOVER_RUNEXONERATE_H

#include "../../impl/impl.h"
#include <thread>
#include <mutex>

void runExonerateEst(std::string& transcriptName, std::string& cdsSequence, std::string& targetSequence,
                     NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix,
                     std::map<std::string, Transcript>& targetTranscriptsHashMap, int& startTarget, int& endTarget,
                     STRAND strand, std::string& tchromeSomeName, std::string& prefixUuid, std::map<std::string, Fasta>& targetGenome, std::map<std::string, std::string>& parameters, int& minIntron);
void readExonerateEstResult( std::string& transcriptName, std::string& cdsSequence, std::string& targetSequence,
                             NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix,
                             std::map<std::string, Transcript>& targetTranscriptsHashMap, int& startTarget, int& endTarget,
                             STRAND& strand, std::string& tchromeSomeName, std::string& fileLocation, std::map<std::string, Fasta>& targetGenome,
                             std::map<std::string, std::string>& parameters, int& minIntron);

void runExonerateProtein(std::string& transcriptName, std::string& protenSequene, std::string& targetSequence,
                         NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix,
                         std::map<std::string, Transcript>& targetTranscriptsHashMap, int& startTarget, int& endTarget,
                         STRAND& strand, std::string& tchromeSomeName, std::string& fileLocation , std::map<std::string, Fasta>& targetGenome,
                         std::map<std::string, std::string>& parameters, int& minIntron );

void readExonerateProteinResult( std::string& fileLocation, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix,
                                 int& startTarget, int& endTarget, STRAND& strand, std::map<std::string, Transcript>& targetTranscriptsHashMap,
                                 std::string& transcriptName, std::string & tchromeSomeName, std::map<std::string, Fasta>& targetGenome, int& minIntron);



#endif //ANNOTATIONLIFTOVER_RUNEXONERATE_H
