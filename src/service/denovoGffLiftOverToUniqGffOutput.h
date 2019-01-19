//
// Created by Baoxing song on 19.10.18.
//

#ifndef ZSDP_DENOVOGFFLIFTOVERTOUNIQGFFOUTPUT_H
#define ZSDP_DENOVOGFFLIFTOVERTOUNIQGFFOUTPUT_H

#include "../impl/impl.h"
#include "../myImportandFunction/myImportantFunction.h"


void denovoGffLiftOverToUniqGffOutput(const std::string & gffFilePath, const std::string & outPutFilePath,
                                      const std::string & queryFastaFilePath, const size_t & minIntron, const size_t & minGene, std::map<std::string, std::string>& parameters);

#endif //ZSDP_DENOVOGFFLIFTOVERTOUNIQGFFOUTPUT_H
