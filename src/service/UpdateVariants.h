//
// Created by song on 8/22/18.
//

#ifndef ZSDP_UPDATEVARIANTS_H
#define ZSDP_UPDATEVARIANTS_H

#include "../impl/impl.h"
#include "../myImportandFunction/myImportantFunction.h"
#include <map>


void updateVariants( const std::string & referenceGffFilePath, const std::string & referenceGenomeFilePath,
                     const std::string & sdiFile, const std::string & vcfFix, const size_t & minIntron,
                     std::map<std::string, std::string>& parameters, const std::string & outPutFilePath  );

#endif //ZSDP_UPDATEVARIANTS_H
