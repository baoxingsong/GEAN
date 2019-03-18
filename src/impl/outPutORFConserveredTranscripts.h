//
// Created by Baoxing Song on 2019-03-17.
//

#ifndef GEAN_OUTPUTORFCONSERVEREDTRANSCRIPTS_H
#define GEAN_OUTPUTORFCONSERVEREDTRANSCRIPTS_H

#include "../model/model.h"
#include "readGffFileWithEverything.h"
#include "readFastaFile.h"
#include "TranscriptUpdateInformation.h"
#include "checkOrfState.h"

void outPutORFConserveredTranscripts( const std::string & genomeFile, const std::string & gffFile,
                                      const std::string & outputGffFile, const int & minIntron,
                                      std::map<std::string, std::string>& parameters );

#endif //GEAN_OUTPUTORFCONSERVEREDTRANSCRIPTS_H
