#include "alignSlidingWindow.h"
#include "gotohAlignment.h"
#include "../util/parameters.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <stack>
#include <stdexcept>

namespace {

GotohAlignmentResult alignNucleotides(
        const std::string& query,
        const std::string& reference,
        int32_t openGapPenalty,
        int32_t extendGapPenalty,
        NucleotideCodeSubstitutionMatrix& substitutionMatrix,
        GotohBoundaryMode boundaryMode,
        int bandWidth = -1) {
    const auto score = [&](size_t referenceIndex, size_t queryIndex) {
        const int referenceCode = substitutionMatrix.get_dna_acid_map(reference[referenceIndex]);
        const int queryCode = substitutionMatrix.get_dna_acid_map(query[queryIndex]);
        return substitutionMatrix.nucleotide_substitution_matrix[referenceCode][queryCode];
    };
    const auto open = [openGapPenalty](size_t) { return openGapPenalty; };
    const auto extend = [extendGapPenalty](size_t) { return extendGapPenalty; };
    return gotohAlign(query, reference, score, open, extend, boundaryMode, bandWidth);
}

void pushAlignment(const GotohAlignmentResult& result,
                   std::stack<char>& queryStack,
                   std::stack<char>& referenceStack) {
    for (size_t i = result.query.size(); i > 0; --i) {
        queryStack.push(result.query[i - 1]);
        referenceStack.push(result.reference[i - 1]);
    }
}

int32_t parameterAsInt32(const std::string& name,
                         std::map<std::string, std::string>& parameters) {
    const long long value = std::stoll(get_parameters(name, parameters));
    if (value < std::numeric_limits<int32_t>::min() ||
        value > std::numeric_limits<int32_t>::max()) {
        throw std::out_of_range(name + " is outside the int32 range");
    }
    return static_cast<int32_t>(value);
}

void appendSlidingAlignment(
        const std::string& query,
        const std::string& reference,
        std::string& alignedQuery,
        std::string& alignedReference,
        int slidingWindowSize,
        int firstWindowExtra,
        int32_t openGapPenalty,
        int32_t extendGapPenalty,
        NucleotideCodeSubstitutionMatrix& substitutionMatrix) {
    if (slidingWindowSize <= 0) {
        throw std::invalid_argument("sliding window size must be positive");
    }
    if (firstWindowExtra < 0) {
        firstWindowExtra = 0;
    }

    size_t queryStart = 0;
    size_t referenceStart = 0;
    size_t currentWindowSize = static_cast<size_t>(slidingWindowSize) +
                               static_cast<size_t>(firstWindowExtra);

    while (queryStart < query.size() && referenceStart < reference.size()) {
        const size_t queryLength = std::min(currentWindowSize, query.size() - queryStart);
        const size_t referenceLength = std::min(currentWindowSize, reference.size() - referenceStart);
        const std::string queryWindow = query.substr(queryStart, queryLength);
        const std::string referenceWindow = reference.substr(referenceStart, referenceLength);

        const GotohAlignmentResult result = alignNucleotides(
                queryWindow, referenceWindow, openGapPenalty, extendGapPenalty,
                substitutionMatrix, GotohBoundaryMode::Overlap);
        if (result.query.size() != result.reference.size()) {
            throw std::runtime_error("alignment strings have different lengths");
        }

        size_t consumedQuery = 0;
        size_t consumedReference = 0;
        for (size_t i = 0; i < result.query.size(); ++i) {
            consumedQuery += result.query[i] != '-';
            consumedReference += result.reference[i] != '-';
        }
        if (consumedQuery == 0 && consumedReference == 0) {
            throw std::runtime_error("sliding alignment made no progress");
        }

        alignedQuery += result.query;
        alignedReference += result.reference;
        queryStart += consumedQuery;
        referenceStart += consumedReference;
        currentWindowSize = static_cast<size_t>(slidingWindowSize);
    }

    while (referenceStart < reference.size()) {
        alignedQuery.push_back('-');
        alignedReference.push_back(reference[referenceStart++]);
    }
    while (queryStart < query.size()) {
        alignedQuery.push_back(query[queryStart++]);
        alignedReference.push_back('-');
    }
}

} // namespace

// These entry points retain their historical names for source compatibility.
// All three now use the same exact-width Gotoh implementation.
void alignment_avx_int32(const std::string& query, const std::string& reference,
                         std::stack<char>& queryStack, std::stack<char>& referenceStack,
                         const int32_t& openGapPenalty, const int32_t& extendGapPenalty,
                         NucleotideCodeSubstitutionMatrix& substitutionMatrix) {
    pushAlignment(alignNucleotides(query, reference, openGapPenalty, extendGapPenalty,
                                   substitutionMatrix, GotohBoundaryMode::Overlap),
                  queryStack, referenceStack);
}

void alignment_avx_int16(const std::string& query, const std::string& reference,
                         std::stack<char>& queryStack, std::stack<char>& referenceStack,
                         const int32_t& openGapPenalty, const int32_t& extendGapPenalty,
                         NucleotideCodeSubstitutionMatrix& substitutionMatrix) {
    alignment_avx_int32(query, reference, queryStack, referenceStack,
                        openGapPenalty, extendGapPenalty, substitutionMatrix);
}

void alignment_avx_int8(const std::string& query, const std::string& reference,
                        std::stack<char>& queryStack, std::stack<char>& referenceStack,
                        const int32_t& openGapPenalty, const int32_t& extendGapPenalty,
                        NucleotideCodeSubstitutionMatrix& substitutionMatrix) {
    alignment_avx_int32(query, reference, queryStack, referenceStack,
                        openGapPenalty, extendGapPenalty, substitutionMatrix);
}

void alignSlidingWindow(const std::string& query, const std::string& reference,
                        std::string& alignedQuery, std::string& alignedReference,
                        const int& slidingWindowSize,
                        std::map<std::string, std::string>& parameters,
                        NucleotideCodeSubstitutionMatrix& substitutionMatrix) {
    appendSlidingAlignment(query, reference, alignedQuery, alignedReference,
                           slidingWindowSize, 0,
                           parameterAsInt32("alignmentOpenGapP", parameters),
                           parameterAsInt32("alignmentExtendGapP", parameters),
                           substitutionMatrix);
}

void alignSlidingWindow(const std::string& query, const std::string& reference,
                        std::string& alignedQuery, std::string& alignedReference,
                        const int& slidingWindowSize,
                        const int& startShiftDistance, const int& endShiftDistance,
                        std::map<std::string, std::string>& parameters,
                        NucleotideCodeSubstitutionMatrix& substitutionMatrix) {
    appendSlidingAlignment(query, reference, alignedQuery, alignedReference,
                           slidingWindowSize,
                           std::max(startShiftDistance, endShiftDistance),
                           parameterAsInt32("alignmentOpenGapP", parameters),
                           parameterAsInt32("alignmentExtendGapP", parameters),
                           substitutionMatrix);
}

void bandAlign(const std::string& query, const std::string& reference,
               std::string& alignedQuery, std::string& alignedReference,
               const int& slidingWindowSize,
               const int& startShiftDistance, const int& endShiftDistance,
               std::map<std::string, std::string>& parameters,
               NucleotideCodeSubstitutionMatrix& substitutionMatrix) {
    const int bandWidth = std::max(slidingWindowSize,
                                   std::max(startShiftDistance, endShiftDistance));
    bandAlign(query, reference, alignedQuery, alignedReference, bandWidth,
              parameters, substitutionMatrix);
}

void bandAlign(const std::string& query, const std::string& reference,
               std::string& alignedQuery, std::string& alignedReference,
               const int& slidingWindowSize,
               std::map<std::string, std::string>& parameters,
               NucleotideCodeSubstitutionMatrix& substitutionMatrix) {
    if (slidingWindowSize < 0) {
        throw std::invalid_argument("band width must not be negative");
    }
    const GotohAlignmentResult result = alignNucleotides(
            query, reference,
            parameterAsInt32("alignmentOpenGapP", parameters),
            parameterAsInt32("alignmentExtendGapP", parameters),
            substitutionMatrix, GotohBoundaryMode::Global, slidingWindowSize);
    alignedQuery += result.query;
    alignedReference += result.reference;
}
