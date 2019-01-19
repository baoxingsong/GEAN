//
// Created by song on 8/26/18.
//

#ifndef ZSDP_ALIGNNEEDLEMANFORTRANSCRIPT_SIMD_AVX2INT16_H
#define ZSDP_ALIGNNEEDLEMANFORTRANSCRIPT_SIMD_AVX2INT16_H



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


class alignNeedlemanForTranscript_simd_avx2int16 {
private:
    std::string _alignment_q;
    std::string _alignment_d;
    std::string _infor;
    //    std::string _signs;

    std::string _dna_q, _dna_d;

    int16_t _exon_match_score;
    int16_t _exon_mismatch_score;
    int16_t _exon_open_gap_penalty;
    int16_t _exon_extend_gap_penalty;
    //std::map<char, <std::map<char, double> > _exon_subsitition_matrix;

    int16_t _intron_match_score;
    int16_t _intron_mismatch_score;
    int16_t _intron_open_gap_penalty;
    int16_t _intron_extend_gap_penalty;

    int16_t _start_stop_codon_match_score;
    int16_t _start_stop_codon_mismatch_score;
    int16_t _start_stop_codon_open_gap_penalty;
    int16_t _start_stop_codon_extend_gap_penalty;

    int16_t _splice_sites_match_score;
    int16_t _splice_sites_mismatch_score;
    int16_t _splice_sites_open_gap_penalty;
    int16_t _splice_sites_extend_gap_penalty;

    int16_t *_similarity_matrix;
    std::vector<bool> *_track_match;
    std::vector<bool> *_track_del;
    std::vector<bool> *_track_ins; //std::vector<bool>  is a special class, in which each element use 1bit RAM space

    //    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix;
    int16_t _startCodonPosition;
    int16_t _stopCodonPosition;
    std::vector<SpliceSitePosition> _spliceSitePositions;

    size_t _length_of_q, _length_of_d, _segLen, _n;
    int8_t bias;
    std::map<std::string, __m256i*> query_profile_avx2_byte( NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix);// dna_q is the query, (the non-reference sequence, with gene structure unknown)
    void calculate_similarity( std::map<std::string, __m256i*> & vProfile, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix);
    void get_optimal_alignment( NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix );
    ELEMENTS checkElements( int position );

public:
    alignNeedlemanForTranscript_simd_avx2int16(std::string& dna_d, std::string& dna_q, int & startCodonPosition, int & stopCodonPosition,
                                               std::vector<SpliceSitePosition>& splitSitePositions, std::map<std::string, std::string>& parameters, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix);
    ~alignNeedlemanForTranscript_simd_avx2int16();
    std::string getAlignment_q();
    std::string getAlignment_d();
    void print_results();
};


#endif //ZSDP_ALIGNNEEDLEMANFORTRANSCRIPT_SIMD_AVX2INT16_H
