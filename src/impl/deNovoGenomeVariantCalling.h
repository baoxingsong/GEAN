//
// Created by Baoxing song on 20.10.18.
//

#ifndef ZSDP_DENOVOGENOMEVARIANTCALLING_H
#define ZSDP_DENOVOGENOMEVARIANTCALLING_H

#include "../model/model.h"
#include "../util/util.h"
#include "./readFastaFile.h"
#include "./readGffFileWithEverything.h"
#include "./organizeGffRecords.h"
#include "../myImportandFunction/myImportantFunction.h"

void deNovoGenomeVariantCalling(const std::string & refGffFilePath, const std::string & refFastaFilePath,
                           const std::string & targetGffFilePath, const std::string & targetFastaFilePath,
                           const size_t & minIntron, const size_t & minGene,
                           std::map<std::string, std::string>& parameters, const size_t & widownWidth, const std::string & outPutFilePath);

#endif //ZSDP_DENOVOGENOMEVARIANTCALLING_H
