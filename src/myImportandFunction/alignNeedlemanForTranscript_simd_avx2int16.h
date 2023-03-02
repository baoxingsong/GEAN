#ifndef GEAN_ALIGN_TRANSCRIPT_SIMD16_H
#define GEAN_ALIGN_TRANSCRIPT_SIMD16_H

#include "alignNeedlemanForTranscript.h"

class alignNeedlemanForTranscript_simd_avx2int16 {
private:
    std::string alignmentQuery;
    std::string alignmentReference;
    std::string alignmentInformation;

public:
    alignNeedlemanForTranscript_simd_avx2int16(
            std::string& reference,
            std::string& query,
            int& startCodonPosition,
            int& stopCodonPosition,
            std::vector<SpliceSitePosition>& spliceSitePositions,
            std::map<std::string, std::string>& parameters,
            NucleotideCodeSubstitutionMatrix& substitutionMatrix);
    ~alignNeedlemanForTranscript_simd_avx2int16() = default;
    std::string getAlignment_q();
    std::string getAlignment_d();
    std::string getInfor();
    void print_results();
};

#endif
