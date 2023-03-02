#include "alignNeedlemanForTranscript_simd_avx2int16.h"

#include <iostream>

alignNeedlemanForTranscript_simd_avx2int16::alignNeedlemanForTranscript_simd_avx2int16(
        std::string& reference,
        std::string& query,
        int& startCodonPosition,
        int& stopCodonPosition,
        std::vector<SpliceSitePosition>& spliceSitePositions,
        std::map<std::string, std::string>& parameters,
        NucleotideCodeSubstitutionMatrix& substitutionMatrix) {
    alignTranscriptGotohExact(reference, query, startCodonPosition, stopCodonPosition,
                              spliceSitePositions, parameters, substitutionMatrix,
                              alignmentQuery, alignmentReference, alignmentInformation);
}

std::string alignNeedlemanForTranscript_simd_avx2int16::getAlignment_q() {
    return alignmentQuery;
}

std::string alignNeedlemanForTranscript_simd_avx2int16::getAlignment_d() {
    return alignmentReference;
}

std::string alignNeedlemanForTranscript_simd_avx2int16::getInfor() {
    return alignmentInformation;
}

void alignNeedlemanForTranscript_simd_avx2int16::print_results() {
    std::cout << alignmentQuery << '\n' << alignmentReference << '\n';
}
