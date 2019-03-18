//
// Created by Baoxing Song on 2019-03-13.
//

#include "OrthologPair.h"

OrthologPair::OrthologPair(const int & queryIndex, const uint32_t & refStartPos, const uint32_t & refEndPos,
        const uint32_t & queryStartPos, const uint32_t & queryEndPos, const double & score,
        const STRAND & strand) : queryIndex(queryIndex), refStartPos(refStartPos), refEndPos(refEndPos),
        queryStartPos(queryStartPos), queryEndPos(queryEndPos), score(score), strand(strand) {
    refMiddlePos = (refStartPos+refEndPos)/2;
    queryMiddlePos = (queryStartPos+queryEndPos)/2;
}

uint32_t OrthologPair::getRefStartPos() const {
    return refStartPos;
}

void OrthologPair::setRefStartPos(uint32_t refStartPos) {
    OrthologPair::refStartPos = refStartPos;
}

uint32_t OrthologPair::getRefEndPos() const {
    return refEndPos;
}

void OrthologPair::setRefEndPos(uint32_t refEndPos) {
    OrthologPair::refEndPos = refEndPos;
}

uint32_t OrthologPair::getRefMiddlePos() const {
    return refMiddlePos;
}

void OrthologPair::setRefMiddlePos(uint32_t refMiddlePos) {
    OrthologPair::refMiddlePos = refMiddlePos;
}

uint32_t OrthologPair::getQueryStartPos() const {
    return queryStartPos;
}

void OrthologPair::setQueryStartPos(uint32_t queryStartPos) {
    OrthologPair::queryStartPos = queryStartPos;
}

uint32_t OrthologPair::getQueryEndPos() const {
    return queryEndPos;
}

void OrthologPair::setQueryEndPos(uint32_t queryEndPos) {
    OrthologPair::queryEndPos = queryEndPos;
}

uint32_t OrthologPair::getQueryMiddlePos() const {
    return queryMiddlePos;
}

void OrthologPair::setQueryMiddlePos(uint32_t queryMiddlePos) {
    OrthologPair::queryMiddlePos = queryMiddlePos;
}

double OrthologPair::getScore() const {
    return score;
}

void OrthologPair::setScore(double score) {
    OrthologPair::score = score;
}

STRAND OrthologPair::getStrand() const {
    return strand;
}

void OrthologPair::setStrand(STRAND & strand) {
    OrthologPair::strand = strand;
}

int OrthologPair::getQueryIndex() const {
    return queryIndex;
}

void OrthologPair::setQueryIndex(int queryIndex) {
    OrthologPair::queryIndex = queryIndex;
}
