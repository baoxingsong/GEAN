#ifndef GEAN_ALIGN_SLIDING_WINDOW_H
#define GEAN_ALIGN_SLIDING_WINDOW_H

#include <map>
#include <string>

#include "../util/nucleotideCodeSubstitutionMatrix.h"

void alignSlidingWindow(
        const std::string& query,
        const std::string& reference,
        std::string& alignedQuery,
        std::string& alignedReference,
        const int& slidingWindowSize,
        const int& startShiftDistance,
        const int& endShiftDistance,
        std::map<std::string, std::string>& parameters,
        NucleotideCodeSubstitutionMatrix& substitutionMatrix);

void alignSlidingWindow(
        const std::string& query,
        const std::string& reference,
        std::string& alignedQuery,
        std::string& alignedReference,
        const int& slidingWindowSize,
        std::map<std::string, std::string>& parameters,
        NucleotideCodeSubstitutionMatrix& substitutionMatrix);

void bandAlign(
        const std::string& query,
        const std::string& reference,
        std::string& alignedQuery,
        std::string& alignedReference,
        const int& slidingWindowSize,
        std::map<std::string, std::string>& parameters,
        NucleotideCodeSubstitutionMatrix& substitutionMatrix);

void bandAlign(
        const std::string& query,
        const std::string& reference,
        std::string& alignedQuery,
        std::string& alignedReference,
        const int& slidingWindowSize,
        const int& startShiftDistance,
        const int& endShiftDistance,
        std::map<std::string, std::string>& parameters,
        NucleotideCodeSubstitutionMatrix& substitutionMatrix);

#endif
