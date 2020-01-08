/*
 * =====================================================================================
 *
 *       Filename:  alignNeedlemanForTranscript.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/02/2017 15:11:41
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Baoxing Song (songbx.me), songbaoxing168@163.com
 *
 * =====================================================================================
 */
/*************************************************************************




 ************************************************************************/
#ifndef _ALIGNNEEDLEMANWUNSCHFORTRANSCRIPT_H
#define _ALIGNNEEDLEMANWUNSCHFORTRANSCRIPT_H

#include <stdlib.h>
#include <string>
#include <map>
#include <vector>
#include "../util/nucleotideCodeSubstitutionMatrix.h"
#include "../model/model.h"
#include <stack>


enum ELEMENTS {
    START, STOP, EXON, INTRON, SPLICEDONOR, SPLICEACCEPTOR
};

class SpliceSitePosition{
    private:
        int _donorSpliceSitePosition;
        int _acceptorSpliceSitePosition;
    public:
        SpliceSitePosition(int donorsSpliceSitePosition, int acceptorSpliceSitePosition);
        int getDonorSpliceSitePosition();
        int getAcceptorSpliceSitePosition();
};

class NeedlemanWunschForTranscript {
    private:
    std::string _alignment_q;
    std::string _alignment_d;
//    std::string _signs;

    std::string _dna_q, _dna_d;
    size_t _length_of_q, _length_of_d;

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

    int32_t **_similarity_matrix;
    VARIANTCATEGORY **_track_matrix;

//    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix;
    int32_t _startCodonPosition;
    int32_t _stopCodonPosition;
    std::vector<SpliceSitePosition> _spliceSitePositions;

    public:
        NeedlemanWunschForTranscript(std::string& dna_d, std::string & dna_q, int & startCodonPosition, int & stopCodonPosition,
            std::vector<SpliceSitePosition>& splitSitePositions, std::map<std::string, std::string>& parameters, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix);
        ~NeedlemanWunschForTranscript();
        std::string getAlignment_q();
        std::string getAlignment_d();
        void setScore(int32_t & match, int32_t & insert, int32_t & del, int & i, int & j);
        void calculate_similarity(NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix);
        void get_optimal_alignment();
        ELEMENTS checkElements( int & position );
        void print_results();
};

void alignNeedlemanForTranscript(std::string& dna_d,
                                 std::string& dna_q, int & startCodonPosition, int & stopCodonPosition,
                                 std::vector<SpliceSitePosition>& spliceSitePositions,
                                 std::map<std::string, std::string>& parameters,
                                 NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix,
                                 std::string & alignment_q, std::string & alignment_d, std::string & infor);
#endif
