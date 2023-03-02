#ifndef GEAN_ALIGN_TRANSCRIPT_SIMD32_H
#define GEAN_ALIGN_TRANSCRIPT_SIMD32_H

#include "alignNeedlemanForTranscript.h"

void alignNeedlemanForTranscript_simd_avx2int32(
        std::string& reference,
        std::string& query,
        int& startCodonPosition,
        int& stopCodonPosition,
        std::vector<SpliceSitePosition>& spliceSitePositions,
        std::map<std::string, std::string>& parameters,
        NucleotideCodeSubstitutionMatrix& substitutionMatrix,
        std::string& alignedQuery,
        std::string& alignedReference,
        std::string& alignmentInformation);

#endif
