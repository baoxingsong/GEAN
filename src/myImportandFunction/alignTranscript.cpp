#include "alignTranscript.h"

AlignTranscript::AlignTranscript(
        std::string& reference,
        std::string& query,
        int& startCodonPosition,
        int& stopCodonPosition,
        std::vector<SpliceSitePosition>& spliceSitePositions,
        std::map<std::string, std::string>& parameters,
        NucleotideCodeSubstitutionMatrix& substitutionMatrix) {
    alignTranscriptGotohExact(reference, query, startCodonPosition, stopCodonPosition,
                              spliceSitePositions, parameters, substitutionMatrix,
                              alignment_q, alignment_d, infor);
}

std::string AlignTranscript::getAlignment_q() {
    return alignment_q;
}

std::string AlignTranscript::getAlignment_d() {
    return alignment_d;
}
