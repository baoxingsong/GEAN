//
// Created by baoxing on 10/10/17.
//

#ifndef ANNOTATIONLIFTOVER_ANNOTATIONLIFTOVER_H
#define ANNOTATIONLIFTOVER_ANNOTATIONLIFTOVER_H
#include "../model/model.h"
#include "../util/util.h"
#include "coordinateLiftOver.h"
#include "checkOrfState.h"
#include "TranscriptUpdateInformation.h"
#include <thread>
#include <chrono>
#include <time.h>

void annotationLiftOver(std::map<std::string, std::vector<Transcript> >& refTranscriptHashSet,
                        std::map<std::string, Transcript >& targetTranscriptHashMap,
                        std::map<std::string, std::vector<Variant> >& variantsMap,
                        std::map<std::string, Fasta>& targetGenome, std::map<std::string, Fasta>& referenceGenome,
                        std::string chromosome, const int& minIntron);

void annotationLiftOver(std::map<std::string, std::vector<Transcript> >& refTranscriptHashSet,
                        std::map<std::string, Transcript >& targetTranscriptHashMap,
                        std::map<std::string, std::vector<Variant> >& variantsMap,
                        std::map<std::string, Fasta>& targetGenome, std::map<std::string, Fasta>& referenceGenome,
                        std::map<std::string, std::string>& parameters, const int& minIntron);

void annotationLiftOver(std::map<std::string, std::vector<Transcript> >& refTranscriptHashSet,
                        std::map<std::string, Transcript >& targetTranscriptHashMap,
                        std::map<std::string, std::vector<Variant> >& variantsMap,
                        std::map<std::string, Fasta>& targetGenome, std::map<std::string, Fasta>& referenceGenome,
                        std::set<std::string>& chromosomes, std::map<std::string, std::string>& parameters, const int& minIntron);

void annotationLiftOverPversion(std::map<std::string, std::vector<Transcript> >& refTranscriptHashSet,
                                std::map<std::string, Transcript >& targetTranscriptHashMap,
                                std::map<std::string, std::vector<Variant> >& variantsMap,
                                std::map<std::string, Fasta>& targetGenome, std::map<std::string, Fasta>& referenceGenome,
                                std::map<std::string, std::string>& parameters, int& minIntron, const int & maxThread);
#endif //ANNOTATIONLIFTOVER_ANNOTATIONLIFTOVER_H
