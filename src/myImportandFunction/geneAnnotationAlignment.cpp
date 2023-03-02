#include "geneAnnotationAlignment.h"

#include <algorithm>
#include <limits>
#include <stack>
#include <stdexcept>

namespace {

enum TraceState : uint8_t {
    TRACE_M = 0,
    TRACE_F = 1,
    TRACE_E = 2,
    TRACE_INVALID = 3
};

const uint8_t M_MASK = 0x03;
const uint8_t F_EXTENDED = 0x04;
const uint8_t E_EXTENDED = 0x08;
const int64_t NEGATIVE_INFINITY = std::numeric_limits<int64_t>::lowest() / 4;

size_t traceIndex(size_t i, size_t j, size_t secondLength) {
    return i * (secondLength + 1) + j;
}

void setMPredecessor(uint8_t& cell, uint8_t state) {
    cell = static_cast<uint8_t>((cell & ~M_MASK) | state);
}

int64_t addScore(int64_t score, int32_t delta) {
    return score == NEGATIVE_INFINITY
                   ? NEGATIVE_INFINITY
                   : score + static_cast<int64_t>(delta);
}

int64_t bestStateScore(int64_t matchScore,
                       int64_t verticalGapScore,
                       int64_t horizontalGapScore,
                       uint8_t& state) {
    int64_t best = matchScore;
    state = TRACE_M;
    if (verticalGapScore > best) {
        best = verticalGapScore;
        state = TRACE_F;
    }
    if (horizontalGapScore > best) {
        best = horizontalGapScore;
        state = TRACE_E;
    }
    return best;
}

} // namespace

int32_t globalAlignment(
        std::vector<Gene>& firstGenes,
        std::vector<Gene>& secondGenes,
        std::vector<std::string>& alignedFirst,
        std::vector<std::string>& alignedSecond) {
    const size_t firstLength = firstGenes.size();
    const size_t secondLength = secondGenes.size();
    const int32_t openGapPenalty = -2;
    const int32_t extendGapPenalty = -1;

    if (firstLength + 1 > std::numeric_limits<size_t>::max() / (secondLength + 1)) {
        throw std::length_error("gene alignment traceback matrix is too large");
    }
    std::vector<uint8_t> trace((firstLength + 1) * (secondLength + 1),
                               TRACE_INVALID);

    std::vector<int64_t> previousM(secondLength + 1, NEGATIVE_INFINITY);
    std::vector<int64_t> currentM(secondLength + 1, NEGATIVE_INFINITY);
    std::vector<int64_t> previousF(secondLength + 1, NEGATIVE_INFINITY);
    std::vector<int64_t> currentF(secondLength + 1, NEGATIVE_INFINITY);
    std::vector<int64_t> previousE(secondLength + 1, NEGATIVE_INFINITY);
    std::vector<int64_t> currentE(secondLength + 1, NEGATIVE_INFINITY);

    previousM[0] = 0;
    setMPredecessor(trace[0], TRACE_M);
    for (size_t j = 1; j <= secondLength; ++j) {
        const size_t index = traceIndex(0, j, secondLength);
        if (j == 1) {
            previousE[j] = addScore(previousM[0], openGapPenalty);
            trace[index] &= static_cast<uint8_t>(~E_EXTENDED);
        } else {
            previousE[j] = addScore(previousE[j - 1], extendGapPenalty);
            trace[index] |= E_EXTENDED;
        }
    }

    for (size_t i = 1; i <= firstLength; ++i) {
        std::fill(currentM.begin(), currentM.end(), NEGATIVE_INFINITY);
        std::fill(currentF.begin(), currentF.end(), NEGATIVE_INFINITY);
        std::fill(currentE.begin(), currentE.end(), NEGATIVE_INFINITY);

        const size_t leftIndex = traceIndex(i, 0, secondLength);
        const int64_t openedLeft = addScore(previousM[0], openGapPenalty);
        const int64_t extendedLeft = addScore(previousF[0], extendGapPenalty);
        if (extendedLeft >= openedLeft) {
            currentF[0] = extendedLeft;
            trace[leftIndex] |= F_EXTENDED;
        } else {
            currentF[0] = openedLeft;
            trace[leftIndex] &= static_cast<uint8_t>(~F_EXTENDED);
        }

        for (size_t j = 1; j <= secondLength; ++j) {
            const size_t index = traceIndex(i, j, secondLength);

            uint8_t matchPredecessor = TRACE_INVALID;
            const int64_t previousBest = bestStateScore(
                    previousM[j - 1], previousF[j - 1], previousE[j - 1],
                    matchPredecessor);
            const int32_t substitution =
                    firstGenes[i - 1].getName() == secondGenes[j - 1].getName()
                            ? 1 : -1;
            currentM[j] = addScore(previousBest, substitution);
            setMPredecessor(trace[index], matchPredecessor);

            const int64_t openedF = addScore(previousM[j], openGapPenalty);
            const int64_t extendedF = addScore(previousF[j], extendGapPenalty);
            if (extendedF >= openedF) {
                currentF[j] = extendedF;
                trace[index] |= F_EXTENDED;
            } else {
                currentF[j] = openedF;
                trace[index] &= static_cast<uint8_t>(~F_EXTENDED);
            }

            const int64_t openedE = addScore(currentM[j - 1], openGapPenalty);
            const int64_t extendedE = addScore(currentE[j - 1], extendGapPenalty);
            if (extendedE >= openedE) {
                currentE[j] = extendedE;
                trace[index] |= E_EXTENDED;
            } else {
                currentE[j] = openedE;
                trace[index] &= static_cast<uint8_t>(~E_EXTENDED);
            }
        }

        previousM.swap(currentM);
        previousF.swap(currentF);
        previousE.swap(currentE);
    }

    uint8_t state = TRACE_INVALID;
    const int64_t finalScore = bestStateScore(
            previousM[secondLength], previousF[secondLength],
            previousE[secondLength], state);
    if (finalScore < std::numeric_limits<int32_t>::min() ||
        finalScore > std::numeric_limits<int32_t>::max()) {
        throw std::overflow_error("gene alignment score exceeds the int32 range");
    }

    std::stack<std::string> firstStack;
    std::stack<std::string> secondStack;
    size_t i = firstLength;
    size_t j = secondLength;
    while (i > 0 || j > 0) {
        if (state == TRACE_M) {
            if (i == 0 || j == 0) {
                throw std::runtime_error("invalid gene M traceback boundary");
            }
            const uint8_t predecessor =
                    trace[traceIndex(i, j, secondLength)] & M_MASK;
            if (predecessor > TRACE_E) {
                throw std::runtime_error("invalid gene M traceback state");
            }
            firstStack.push(firstGenes[i - 1].getName());
            secondStack.push(secondGenes[j - 1].getName());
            --i;
            --j;
            state = predecessor;
        } else if (state == TRACE_F) {
            if (i == 0) {
                throw std::runtime_error("invalid gene F traceback boundary");
            }
            const uint8_t predecessor =
                    (trace[traceIndex(i, j, secondLength)] & F_EXTENDED)
                            ? TRACE_F : TRACE_M;
            firstStack.push(firstGenes[i - 1].getName());
            secondStack.push("-");
            --i;
            state = predecessor;
        } else if (state == TRACE_E) {
            if (j == 0) {
                throw std::runtime_error("invalid gene E traceback boundary");
            }
            const uint8_t predecessor =
                    (trace[traceIndex(i, j, secondLength)] & E_EXTENDED)
                            ? TRACE_E : TRACE_M;
            firstStack.push("-");
            secondStack.push(secondGenes[j - 1].getName());
            --j;
            state = predecessor;
        } else {
            throw std::runtime_error("invalid gene traceback state");
        }
    }

    while (!firstStack.empty()) {
        alignedFirst.push_back(firstStack.top());
        alignedSecond.push_back(secondStack.top());
        firstStack.pop();
        secondStack.pop();
    }
    return static_cast<int32_t>(finalScore);
}
