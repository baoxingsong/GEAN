#include "alignNeedlemanForTranscript.h"
#include "gotohAlignment.h"

#include <array>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>

#include "../util/parameters.h"

namespace {

int32_t parameterAsInt32(const std::string& name,
                         std::map<std::string, std::string>& parameters) {
    const long long value = std::stoll(get_parameters(name, parameters));
    if (value < std::numeric_limits<int32_t>::min() ||
        value > std::numeric_limits<int32_t>::max()) {
        throw std::out_of_range(name + " is outside the int32 range");
    }
    return static_cast<int32_t>(value);
}

ELEMENTS classifyReferencePosition(
        int position,
        int startCodonPosition,
        int stopCodonPosition,
        const std::vector<SpliceSitePosition>& spliceSitePositions) {
    if (position < startCodonPosition || position > stopCodonPosition + 2) {
        return INTRON;
    }
    for (const SpliceSitePosition& spliceSite : spliceSitePositions) {
        if (position > spliceSite.getDonorSpliceSitePosition() + 1 &&
            position < spliceSite.getAcceptorSpliceSitePosition() - 1) {
            return INTRON;
        }
        if (position == spliceSite.getAcceptorSpliceSitePosition() ||
            position == spliceSite.getAcceptorSpliceSitePosition() - 1) {
            return SPLICEACCEPTOR;
        }
        if (position == spliceSite.getDonorSpliceSitePosition() ||
            position == spliceSite.getDonorSpliceSitePosition() + 1) {
            return SPLICEDONOR;
        }
    }
    if (position == startCodonPosition || position == startCodonPosition + 1 ||
        position == startCodonPosition + 2) {
        return START;
    }
    if (position == stopCodonPosition || position == stopCodonPosition + 1 ||
        position == stopCodonPosition + 2) {
        return STOP;
    }
    return EXON;
}

size_t scoreClass(ELEMENTS element) {
    switch (element) {
        case EXON:
            return 0;
        case INTRON:
            return 1;
        case SPLICEDONOR:
        case SPLICEACCEPTOR:
            return 2;
        case START:
        case STOP:
            return 3;
    }
    throw std::runtime_error("unknown reference element");
}

char informationCharacter(ELEMENTS element) {
    switch (element) {
        case EXON:
            return 'E';
        case INTRON:
            return 'I';
        case SPLICEDONOR:
        case SPLICEACCEPTOR:
            return 'P';
        case START:
            return 'S';
        case STOP:
            return 'T';
    }
    return '?';
}

} // namespace

SpliceSitePosition::SpliceSitePosition(int donorPosition, int acceptorPosition)
        : donorSpliceSitePosition(donorPosition),
          acceptorSpliceSitePosition(acceptorPosition) {}

int SpliceSitePosition::getDonorSpliceSitePosition() const {
    return donorSpliceSitePosition;
}

int SpliceSitePosition::getAcceptorSpliceSitePosition() const {
    return acceptorSpliceSitePosition;
}

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
        std::string& alignmentInformation) {
    const std::array<int32_t, 4> openGapPenalties = {{
            parameterAsInt32("alignmentExonOpenGapP", parameters),
            parameterAsInt32("alignmentIntronOpenGapP", parameters),
            parameterAsInt32("alignmentSpliceSitesOpenGapP", parameters),
            parameterAsInt32("alignmentStartStopCodonOpenGapP", parameters)}};
    const std::array<int32_t, 4> extendGapPenalties = {{
            parameterAsInt32("alignmentExonExtendGapP", parameters),
            parameterAsInt32("alignmentIntronExtendGapP", parameters),
            parameterAsInt32("alignmentSpliceSitesExtendGapP", parameters),
            parameterAsInt32("alignmentStartStopCodonExtendGapP", parameters)}};

    std::vector<ELEMENTS> referenceElements(reference.size());
    for (size_t i = 0; i < reference.size(); ++i) {
        referenceElements[i] = classifyReferencePosition(
                static_cast<int>(i + 1), startCodonPosition, stopCodonPosition,
                spliceSitePositions);
    }

    const auto score = [&](size_t referenceIndex, size_t queryIndex) {
        const int referenceCode = substitutionMatrix.get_dna_acid_map(reference[referenceIndex]);
        const int queryCode = substitutionMatrix.get_dna_acid_map(query[queryIndex]);
        switch (scoreClass(referenceElements[referenceIndex])) {
            case 0:
                return substitutionMatrix.get_exon_subsitition_matrix()[referenceCode][queryCode];
            case 1:
                return substitutionMatrix.get_intron_subsitition_matrix()[referenceCode][queryCode];
            case 2:
                return substitutionMatrix.get_splice_sites_subsitition_matrix()[referenceCode][queryCode];
            default:
                return substitutionMatrix.get_start_stop_codon_subsitition_matrix()[referenceCode][queryCode];
        }
    };
    const auto open = [&](size_t referenceIndex) {
        return openGapPenalties[scoreClass(referenceElements[referenceIndex])];
    };
    const auto extend = [&](size_t referenceIndex) {
        return extendGapPenalties[scoreClass(referenceElements[referenceIndex])];
    };

    const GotohAlignmentResult result = gotohAlign(
            query, reference, score, open, extend,
            GotohBoundaryMode::ReferenceGlobal);
    alignedQuery += result.query;
    alignedReference += result.reference;

    size_t referencePosition = 0;
    for (size_t i = 0; i < result.reference.size(); ++i) {
        if (result.reference[i] == '-') {
            alignmentInformation.push_back('-');
        } else {
            alignmentInformation.push_back(
                    informationCharacter(referenceElements[referencePosition++]));
        }
    }
}

NeedlemanWunschForTranscript::NeedlemanWunschForTranscript(
        std::string& reference,
        std::string& query,
        int& startCodonPosition,
        int& stopCodonPosition,
        std::vector<SpliceSitePosition>& spliceSitePositions,
        std::map<std::string, std::string>& parameters,
        NucleotideCodeSubstitutionMatrix& substitutionMatrix) {
    std::string information;
    alignTranscriptGotohExact(reference, query, startCodonPosition, stopCodonPosition,
                              spliceSitePositions, parameters, substitutionMatrix,
                              alignmentQuery, alignmentReference, information);
}

std::string NeedlemanWunschForTranscript::getAlignment_q() {
    return alignmentQuery;
}

std::string NeedlemanWunschForTranscript::getAlignment_d() {
    return alignmentReference;
}

void NeedlemanWunschForTranscript::print_results() {
    std::cout << alignmentQuery << '\n' << alignmentReference << '\n';
}

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
        std::string& alignmentInformation) {
    alignTranscriptGotohExact(reference, query, startCodonPosition, stopCodonPosition,
                              spliceSitePositions, parameters, substitutionMatrix,
                              alignedQuery, alignedReference, alignmentInformation);
}
