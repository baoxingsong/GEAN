//
// Created by song on 8/26/18.
//

#include "alignTranscript.h"

AlignTranscript::AlignTranscript(std::string& dna_d, std::string& dna_q, int & startCodonPosition, int & stopCodonPosition,
        std::vector<SpliceSitePosition>& splitSitePositions, std::map<std::string, std::string>& parameters, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix){
#ifdef __AVX2__
    if( dna_d.length() > 3270 || dna_q.length()>3270 ){
        alignNeedlemanForTranscript_simd_avx2int32(dna_d, dna_q, startCodonPosition, stopCodonPosition, splitSitePositions, parameters, nucleotideCodeSubstitutionMatrix,alignment_q, alignment_d, infor);
    }else{
        avx16 = new alignNeedlemanForTranscript_simd_avx2int16(dna_d, dna_q, startCodonPosition, stopCodonPosition, splitSitePositions, parameters, nucleotideCodeSubstitutionMatrix);
        alignment_q = avx16->getAlignment_q();
        alignment_d = avx16->getAlignment_d();
        infor = avx16->getInfor();
    }
#else
    alignNeedlemanForTranscript(dna_d, dna_q, startCodonPosition, stopCodonPosition, splitSitePositions, parameters, nucleotideCodeSubstitutionMatrix,alignment_q, alignment_d, infor);
#endif
}

AlignTranscript::~AlignTranscript(){
    if( avx16 != NULL ){
        delete(avx16);
    }
}

std::string AlignTranscript::getAlignment_q(){
    return this->alignment_q;
}

std::string AlignTranscript::getAlignment_d(){
    return this->alignment_d;
}
