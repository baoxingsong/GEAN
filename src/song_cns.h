//
// Created by Baoxing song on 2018-12-28.
//

#ifndef GEAN_SONG_CNS_H
#define GEAN_SONG_CNS_H

#include <string>
#include <map>
#include <set>
#include "./util/nucleotideCodeSubstitutionMatrix.h"
#include <stack>
#include<iostream>
#include <cstring>
#include <fstream>
#include <cstdlib>
#include <limits>
#include <sstream>
#include <iomanip>


class SimilarBlocks {
private:
    std::string _dna_q, _dna_d;
    size_t _length_of_q, _length_of_d, _n;
    int32_t  _open_gap_penalty;
    int32_t  _extend_gap_penalty;
    int32_t  **_similarity_matrix;
    int32_t  **_similarity_matrix_E;
    int32_t  **_similarity_matrix_F;
    int32_t  **_substitute_matrix;

    std::string _alignment_q;
    std::string _alignment_d;
public:
    SimilarBlocks(const std::string & dna_q, const std::string & dna_d, const int32_t & match_score, const int32_t & mis_match_score,
                    const int32_t & open_gap_penalty, const int32_t & extend_gap_penalty,
                    NucleotideCodeSubstitutionMatrix& _nucleotideCodeSubstitutionMatrix);
    ~SimilarBlocks();

    void calculate_similarity( NucleotideCodeSubstitutionMatrix& _nucleotideCodeSubstitutionMatrix );
    void get_optimal_alignment();
    std::string getAlignment_q();
    std::string getAlignment_d();
    void print_results();
};


#endif //GEAN_SONG_CNS_H
