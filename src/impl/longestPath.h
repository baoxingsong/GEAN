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

#endif //GEAN_LONGESTPATH_H
