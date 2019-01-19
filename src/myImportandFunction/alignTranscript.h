//
// Created by song on 8/26/18.
//

#ifndef ZSDP_ALIGNTRANSCRIPT_H
#define ZSDP_ALIGNTRANSCRIPT_H

#include "alignNeedlemanForTranscript_simd_avx2int32.h"
#include "alignNeedlemanForTranscript_simd_avx2int16.h"


class AlignTranscript {
    private:
        alignNeedlemanForTranscript_simd_avx2int32 * avx32 = NULL;
        alignNeedlemanForTranscript_simd_avx2int16 * avx16 = NULL;
        std::string model;
    public:
        AlignTranscript(std::string& dna_d, std::string& dna_q, int & startCodonPosition, int & stopCodonPosition,
                        std::vector<SpliceSitePosition>& splitSitePositions, std::map<std::string, std::string>& parameters,
                        NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix);
        ~AlignTranscript();
        std::string getAlignment_q();
        std::string getAlignment_d();
        void print_results();
};


#endif //ZSDP_ALIGNTRANSCRIPT_H
