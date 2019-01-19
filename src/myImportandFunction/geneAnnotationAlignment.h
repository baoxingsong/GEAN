//
// Created by Baoxing song on 09.10.18.
//

#ifndef ZSDP_GENEANNOTATIONALIGNMENT_H
#define ZSDP_GENEANNOTATIONALIGNMENT_H

#include "../model/model.h"
#include "../util/util.h"
#include <stack>

int32_t globalAlignment( std::vector<Gene> & g1s, std::vector<Gene> & g2s, std::vector<std::string> & _alignment_q,
    std::vector<std::string> & _alignment_d);

#endif //ZSDP_GENEANNOTATIONALIGNMENT_H
