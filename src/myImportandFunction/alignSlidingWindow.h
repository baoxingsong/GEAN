//
// Created by song on 8/5/18.
//

#ifndef ANNOTATIONLIFTOVER_ALIGNSLIDINGWINDOW_H
#define ANNOTATIONLIFTOVER_ALIGNSLIDINGWINDOW_H

#include <string>
#include <stack>
#include <immintrin.h>
#include <map>
#include "../impl/impl.h"
#include "../util/nucleotideCodeSubstitutionMatrix.h"
#include "../util/parameters.h"
#include "hengAlign.h"
void alignSlidingWindow( const std::string& dna_q, const std::string& dna_d,
                         std::string & _alignment_q, std::string & _alignment_d, const int & slidingWindowSize,
                         std::map<std::string, std::string>& parameters, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix );

#endif //ANNOTATIONLIFTOVER_ALIGNSLIDINGWINDOW_H
