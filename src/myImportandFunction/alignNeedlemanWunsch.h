/*
 * =====================================================================================
 *
 *       Filename:  alignNeedlemanWunsch.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/02/2017 15:11:41
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Baoxing Song (songbx.me), songbaoxing168@163.com
 *
 * =====================================================================================
 */
/*************************************************************************




 ************************************************************************/

#ifndef _ALIGNNEEDLEMANWUNSCH_H
#define _ALIGNNEEDLEMANWUNSCH_H

#include <string>
#include <map>
#include <set>
#include "../util/nucleotideCodeSubstitutionMatrix.h"
#include <stack>

class NeedlemanWunsch {
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
        VARIANTCATEGORY **_track_matrix;
    public:
        NeedlemanWunsch(const std::string & dna_q, const std::string & dna_d, const int32_t & match_score, const int32_t & mis_match_score,
                    const int32_t & open_gap_penalty, const int32_t & extend_gap_penalty,
                    NucleotideCodeSubstitutionMatrix& _nucleotideCodeSubstitutionMatrix);
        ~NeedlemanWunsch();

        void calculate_similarity( NucleotideCodeSubstitutionMatrix& _nucleotideCodeSubstitutionMatrix );
    void get_optimal_alignment();
    std::string getAlignment_q();
    std::string getAlignment_d();
    void print_results();
};

#endif
