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
#include "math.h"
void generateLongestOutput( const std::string & referenceGffFile, const std::string & queryNewGffFile,
                            const std::string & queryGenomeFile, const std::string & outputGffFile, const int & minIntron,
                            double score, double penalty, double scoreThreshold, const bool & keepTandemDuplication,
                            std::map<std::string, std::string>& parameters, const double & syntenicScore,
                            const double & orfScore, const double & dropLengthThredshold, const bool & onlySyntenic );



void generateLongestQuotaOutput( const std::string & referenceGffFile, const std::string & queryNewGffFile,
                                 const std::string & queryGenomeFile, const std::string & outputGffFile, const int & minIntron,
                                 const bool & keepTandemDuplication,
                                 std::map<std::string, std::string>& parameters, const double & syntenicScore,
                                 const double & orfScore, const double & dropLengthThredshold,
                                 double & INDEL_SCORE, double & GAP_OPEN_PENALTY, double & MIN_ALIGNMENT_SCORE, int & MAX_DIST_BETWEEN_MATCHES,
                                 int & refMaximumTimes, int & queryMaximumTimes, const bool & onlySyntenic, const bool & sortOutPutGffBycoordinate );

void generateDagChainerOutput( const std::string & referenceGffFile, const std::string & queryNewGffFile,
                               const std::string & queryGenomeFile, const std::string & outputGffFile, const int & minIntron,
                               const bool & keepTandemDuplication,
                               std::map<std::string, std::string>& parameters, const double & syntenicScore, const double & orfScore,
                               const double & dropLengthThredshold, int MAX_DIST_BETWEEN_MATCHES /*max gap in the term of number of genes*/,
                               int BP_GAP_SIZE, double INDEL_SCORE, double GAP_OPEN_PENALTY, double MIN_ALIGNMENT_SCORE, const bool & onlySyntenic, const bool & sortOutPutGffBycoordinate);

#endif //GEAN_LONGESTPATH_H
