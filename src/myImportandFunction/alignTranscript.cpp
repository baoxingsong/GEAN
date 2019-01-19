//
// Created by song on 8/26/18.
//

#include "alignTranscript.h"

AlignTranscript::AlignTranscript(std::string& dna_d, std::string& dna_q, int & startCodonPosition, int & stopCodonPosition,
        std::vector<SpliceSitePosition>& splitSitePositions, std::map<std::string, std::string>& parameters, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix){
    if( dna_d.length() > 3270 || dna_q.length()>3270 ){
        avx32 = new alignNeedlemanForTranscript_simd_avx2int32(dna_d, dna_q, startCodonPosition, stopCodonPosition, splitSitePositions, parameters, nucleotideCodeSubstitutionMatrix);
    }else{
        avx16 = new alignNeedlemanForTranscript_simd_avx2int16(dna_d, dna_q, startCodonPosition, stopCodonPosition, splitSitePositions, parameters, nucleotideCodeSubstitutionMatrix);
    }
}

AlignTranscript::~AlignTranscript(){
    if( avx32 != NULL ){
        delete(avx32);
    }
    if( avx16 != NULL ){
        delete(avx16);
    }
}

std::string AlignTranscript::getAlignment_q(){
    if( avx32 != NULL ){
        return avx32->getAlignment_q();
    }
    if( avx16 != NULL ){
        return avx16->getAlignment_q();
    }
    return "";
}

std::string AlignTranscript::getAlignment_d(){
    if( avx32 != NULL ){
        return avx32->getAlignment_d();
    }
    if( avx16 != NULL ){
        return avx16->getAlignment_d();
    }
    return "";
}

void AlignTranscript::print_results(){
    if( avx32 != NULL ){
        avx32->print_results();
    }else if( avx16 != NULL ){
        avx16->print_results();
    }
}

