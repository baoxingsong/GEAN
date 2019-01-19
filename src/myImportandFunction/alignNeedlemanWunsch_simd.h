//
// Created by Baoxing song on 25.04.18.
//

#ifndef ANNOTATIONLIFTOVER_ALIGNNEEDLEMANWUNSCH_SIMD_H
#define ANNOTATIONLIFTOVER_ALIGNNEEDLEMANWUNSCH_SIMD_H
#include <string>
#include <map>
#include <cmath>
#include <stack>
#include <iostream>
#include <cstring>
#include "../util/nucleotideCodeSubstitutionMatrix.h"
#include <immintrin.h>


#include <fstream>
#include <cstdlib>
#include <limits>
#include <sstream>
#include <iomanip>


class NeedlemanWunsch_simd_Slow {
    private:
        std::string _dna_q, _dna_d;
        size_t _length_of_q, _length_of_d, _segLen, _n;
        int32_t  _open_gap_penalty;
        int32_t  _extend_gap_penalty;
        int32_t  **_similarity_matrix;
        int32_t  **_substitute_matrix;
        int32_t **_track_match;
        int32_t **_track_del;
        int32_t **_track_ins;

        int32_t bias; // here bias should be 0, it is used to score the part beyond the validated sequence

        std::string _alignment_q;
        std::string _alignment_d;
    public:
        NeedlemanWunsch_simd_Slow(const std::string & dna_q, const std::string & dna_d, const int8_t & match_score,
                             const int8_t & mis_match_score, const int8_t & open_gap_penalty,
                             const int8_t & extend_gap_penalty,  NucleotideCodeSubstitutionMatrix & _nucleotideCodeSubstitutionMatrix);
        ~NeedlemanWunsch_simd_Slow();
        __m256i* query_profile_avx2_byte(  NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix);
        void ssw_avx2( __m256i* vProfile,  NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix );
        int32_t ** getSimilarity_matrix();
        void get_optimal_alignment( );
        void print_results();
        std::string getAlignment_q();
        std::string getAlignment_d();
};


class NeedlemanWunsch_simd_Fast {
private:
    std::string _dna_q, _dna_d;
    size_t _length_of_q, _length_of_d, _segLen, _n;
    int8_t  _open_gap_penalty;
    int8_t  _extend_gap_penalty;
    int32_t  **_similarity_matrix;
    int32_t  **_substitute_matrix;
    int8_t bias; // here bias should be 0, it is used to score the part beyond the validated sequence

    std::string _alignment_q;
    std::string _alignment_d;

public:
    NeedlemanWunsch_simd_Fast(const std::string & dna_q, const std::string & dna_d, const int8_t & match_score,
                         const int8_t & mis_match_score, const int8_t & open_gap_penalty,
                         const int8_t & extend_gap_penalty, const NucleotideCodeSubstitutionMatrix & _nucleotideCodeSubstitutionMatrix);
    ~NeedlemanWunsch_simd_Fast();
    __m256i* query_profile_avx2_byte( NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix);
    void ssw_avx2( __m256i* vProfile, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix );
    int32_t ** getSimilarity_matrix();
    void get_optimal_alignment(NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix);
    void print_results();
    std::string get_alignment_q();
    std::string get_alignment_d();
};

#endif //ANNOTATIONLIFTOVER_ALIGNNEEDLEMANWUNSCH_SIMD_H
