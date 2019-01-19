//
// Created by song on 8/4/18.
//

#include "Range.h"

Range::Range(const std::string &chr, const size_t & start, const size_t & end, const STRAND & strand) : chr(chr), start(start), end(end),
                                                                                strand(strand) {}
Range::Range(){
    chr="";
    start=0;
    end=0;
    strand=POSITIVE;
}

const std::string &Range::getChr() const {
    return chr;
}

void Range::setChr(const std::string &chr) {
    Range::chr = chr;
}

const size_t & Range::getStart() const {
    return start;
}

void Range::setStart(const size_t & start) {
    Range::start = start;
}

const size_t & Range::getEnd() const {
    return end;
}

void Range::setEnd(const size_t & end) {
    Range::end = end;
}

const STRAND & Range::getStrand() const {
    return strand;
}

void Range::setStrand(const STRAND & strand) {
    Range::strand = strand;
}
