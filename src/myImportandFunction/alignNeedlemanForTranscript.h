#ifndef GEAN_ALIGN_NEEDLEMAN_FOR_TRANSCRIPT_H
#define GEAN_ALIGN_NEEDLEMAN_FOR_TRANSCRIPT_H

#include <map>
#include <string>
#include <vector>

#include "../util/nucleotideCodeSubstitutionMatrix.h"

enum ELEMENTS {
    START,
    STOP,
    EXON,
    INTRON,
    SPLICEDONOR,
    SPLICEACCEPTOR
};

class SpliceSitePosition {
private:
    int donorSpliceSitePosition;
    int acceptorSpliceSitePosition;

public:
    SpliceSitePosition(int donorSpliceSitePosition, int acceptorSpliceSitePosition);
    int getDonorSpliceSitePosition() const;
    int getAcceptorSpliceSitePosition() const;
};

void alignTranscriptGotohExact(
        const std::string& reference,
        const std::string& query,
        int startCodonPosition,
        int stopCodonPosition,
        const std::vector<SpliceSitePosition>& spliceSitePositions,
        std::map<std::string, std::string>& parameters,
        NucleotideCodeSubstitutionMatrix& substitutionMatrix,
        std::string& alignedQuery,
        std::string& alignedReference,
        std::string& alignmentInformation);

class NeedlemanWunschForTranscript {
private:
    std::string alignmentQuery;
    std::string alignmentReference;

public:
    NeedlemanWunschForTranscript(
            std::string& reference,
            std::string& query,
            int& startCodonPosition,
            int& stopCodonPosition,
            std::vector<SpliceSitePosition>& spliceSitePositions,
            std::map<std::string, std::string>& parameters,
            NucleotideCodeSubstitutionMatrix& substitutionMatrix);
    ~NeedlemanWunschForTranscript() = default;
    std::string getAlignment_q();
    std::string getAlignment_d();
    void print_results();
};

void alignNeedlemanForTranscript(
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
