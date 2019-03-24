//
// Created by Baoxing Song on 2019-03-13.
//

#ifndef GEAN_LONGESTPATH_H
#define GEAN_LONGESTPATH_H

#include "../model/model.h"
#include "readGffFileWithEverything.h"
#include "readFastaFile.h"
#include "TranscriptUpdateInformation.h"
#include "checkOrfState.h"

void generateLongestOutput( const std::string & referenceGffFile, const std::string & queryNewGffFile,
                            const std::string & queryGenomeFile, const std::string & outputGffFile, const int & minIntron,
                            double score, double penalty, double scoreThreshold, const bool & keepTandemDuplication,
                            std::map<std::string, std::string>& parameters, const double & syntenicScore, const double & orfScore, const double & dropLengthThredshold );



void generateLongestQuotaOutput( const std::string & referenceGffFile, const std::string & queryNewGffFile,
                                 const std::string & queryGenomeFile, const std::string & outputGffFile, const int & minIntron,
                                 const bool & keepTandemDuplication,
                                 std::map<std::string, std::string>& parameters, const double & syntenicScore,
                                 const double & orfScore, const double & dropLengthThredshold,
                                 double & INDEL_SCORE, double & GAP_OPEN_PENALTY, double & MIN_ALIGNMENT_SCORE,
                                 int & refMaximumTimes, int & queryMaximumTimes );

void generateDagChainerOutput( const std::string & referenceGffFile, const std::string & queryNewGffFile,
                               const std::string & queryGenomeFile, const std::string & outputGffFile, const int & minIntron,
                               const bool & keepTandemDuplication,
                               std::map<std::string, std::string>& parameters, const double & syntenicScore, const double & orfScore,
                               const double & dropLengthThredshold, int MAX_DIST_BETWEEN_MATCHES /*max gap in the term of number of genes*/,
                               int BP_GAP_SIZE, double INDEL_SCORE, double GAP_OPEN_PENALTY, double MIN_ALIGNMENT_SCORE);

#endif //GEAN_LONGESTPATH_H
