//
// Created by song on 8/4/18.
//

#ifndef ANNOTATIONLIFTOVER_RANGE_H
#define ANNOTATIONLIFTOVER_RANGE_H

#include <string>
#include "STRAND.h"
class Range {
private:
    std::string chr;
    size_t start;
    size_t end;
    STRAND strand;
public:
    Range(const std::string &chr, const size_t & start, const size_t & end, const STRAND & strand);
    Range();
    const std::string &getChr() const;
    void setChr(const std::string &chr);
    const size_t & getStart() const;
    void setStart(const size_t & start);
    const size_t & getEnd() const;
    void setEnd(const size_t & end);
    const STRAND & getStrand() const;
    void setStrand(const STRAND & strand);
};


#endif //ANNOTATIONLIFTOVER_RANGE_H
