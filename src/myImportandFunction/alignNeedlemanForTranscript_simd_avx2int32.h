//
// Created by song on 7/27/18.
//

#ifndef ANNOTATIONLIFTOVER_ALIGNNEEDLEMANFORTRANSCRIPT_SIMD_AVX2INT32_H
#define ANNOTATIONLIFTOVER_ALIGNNEEDLEMANFORTRANSCRIPT_SIMD_AVX2INT32_H

#include "alignNeedlemanForTranscript.h"
#include <immintrin.h>
#include <stack>
#include <bitset>



#include <iostream>
#include <map>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "../util/parameters.h"
#include <algorithm>


class alignNeedlemanForTranscript_simd_avx2int32 {
    private:
        std::string _alignment_q;
        std::string _alignment_d;
        std::string _infor;
    //    std::string _signs;

        std::string _dna_q, _dna_d;

        int32_t _exon_match_score;
        int32_t _exon_mismatch_score;

        //std::map<char, <std::map<char, double> > _exon_subsitition_matrix;

        int32_t _intron_match_score;
        int32_t _intron_mismatch_score;

        int32_t _start_stop_codon_match_score;
        int32_t _start_stop_codon_mismatch_score;

        int32_t _splice_sites_match_score;
        int32_t _splice_sites_mismatch_score;

        int32_t *_open_gap_penalty;
        int32_t *_extend_gap_penalty;

        int8_t *_ref_elements;

        int32_t *_similarity_matrix;
        std::vector<bool> *_trace_match;
        std::vector<bool> *_trace_del;
        std::vector<bool> *_trace_ins;
        //std::vector<bool> *_trace_ins; //std::vector<bool>  is a special class, in which each element use 1bit RAM space

    //    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix;
        int32_t _startCodonPosition;
        int32_t _stopCodonPosition;
        std::vector<SpliceSitePosition> _spliceSitePositions;

        size_t _length_of_q, _length_of_d, _segLen, _n;
        int8_t bias;
        __m256i** query_profile_avx2_byte( NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix);// dna_q is the query, (the non-reference sequence, with gene structure unknown)
        void calculate_similarity_not_very_good( __m256i** vProfile, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix);
        void calculate_similarity( __m256i** vProfile, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix);
        void get_optimal_alignment( NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix );
        int checkElements( int position );
public:
        alignNeedlemanForTranscript_simd_avx2int32(std::string& dna_d, std::string& dna_q, int & startCodonPosition, int & stopCodonPosition,
                                     std::vector<SpliceSitePosition>& splitSitePositions, std::map<std::string, std::string>& parameters, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix);
        ~alignNeedlemanForTranscript_simd_avx2int32();
        std::string getAlignment_q();
        std::string getAlignment_d();
        void print_results();
};

#endif //ANNOTATIONLIFTOVER_ALIGNNEEDLEMANFORTRANSCRIPT_SIMD_AVX2INT32_H
