//
// Created by song on 7/27/18.
//

#ifndef ANNOTATIONLIFTOVER_ALIGNNEEDLEMANFORTRANSCRIPT_SIMD_AVX2INT32_H
#define ANNOTATIONLIFTOVER_ALIGNNEEDLEMANFORTRANSCRIPT_SIMD_AVX2INT32_H

#include "alignNeedlemanForTranscript.h"
#include <immintrin.h>
#include <stack>
#include <bitset>



#include <iostream>
#include <map>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "../util/parameters.h"
#include <algorithm>

void alignNeedlemanForTranscript_simd_avx2int32(std::string& dna_d,
                                                std::string& dna_q, int & startCodonPosition, int & stopCodonPosition,
                                                std::vector<SpliceSitePosition>& spliceSitePositions,
                                                std::map<std::string, std::string>& parameters,
                                                NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix,
                                                std::string & alignment_q, std::string & alignment_d, std::string & infor);

#endif //ANNOTATIONLIFTOVER_ALIGNNEEDLEMANFORTRANSCRIPT_SIMD_AVX2INT32_H
