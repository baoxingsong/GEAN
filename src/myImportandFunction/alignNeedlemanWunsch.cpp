/*
 * =====================================================================================
 *
 *       Filename:  alignNeedlemanWunsch.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/02/2017 15:11:39
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Baoxing Song (songbx.me), songbaoxing168@163.com
 *
 * =====================================================================================
 */
/*************************************************************************




 ************************************************************************/

#include<iostream>
#include <map>
#include "alignNeedlemanWunsch.h"
#include "../util/nucleotideCodeSubstitutionMatrix.h"
#include <string>
#include <cstring>
#include <fstream>
#include <cstdlib>
#include <limits>
#include <sstream>
#include <iomanip>

NeedlemanWunsch::NeedlemanWunsch(const std::string & dna_q, const std::string & dna_d, const int32_t & match_score, const int32_t & mis_match_score,
                                 const int32_t & open_gap_penalty, const int32_t & extend_gap_penalty,
                                 NucleotideCodeSubstitutionMatrix& _nucleotideCodeSubstitutionMatrix) {
    this->_dna_q = dna_q;
    this->_dna_d = dna_d;
    this->_open_gap_penalty = open_gap_penalty;
    this->_extend_gap_penalty = extend_gap_penalty;

    this->_length_of_q = dna_q.length();
    this->_length_of_d = dna_d.length();
    this->_n = 5;

    this->_similarity_matrix = new int32_t *[this->_length_of_q + 1];

    size_t  i, j;
    for (i = 0; i < (this->_length_of_q + 1); ++i) {
        this->_similarity_matrix[i] = new int32_t [this->_length_of_d + 1];
    }

    i=0;
    for ( j = 0; j <= _length_of_d; ++j) {
        this->_similarity_matrix[i][j] = 0;
    }
    j=0;
    for (i = 0; i < (this->_length_of_q+1); ++i) {
        this->_similarity_matrix[i][j] = 0;
    }

    // this matrix is for set different penalty for open gap and extend gap begin
    // and the track also changed to use this matrix
    this->_track_matrix = new VARIANTCATEGORY*[this->_length_of_q + 1];
    for (i = 0; i < (this->_length_of_q + 1); ++i) {
        _track_matrix[i] = new VARIANTCATEGORY[this->_length_of_d + 1];
    }

    for ( i = 0; i <= _length_of_q; ++i) {
        for ( j = 0; j <= _length_of_d; ++j) {
            _track_matrix[i][j] = SNPORINSERTIONORDELETION;
        }
    }

    this->_substitute_matrix = new int32_t *[_n];
    for( i=0; i<_n; ++i){
        this->_substitute_matrix[i] = new int32_t [_n];
        for( j=0; j<_n; ++j){
            if( _nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[i][j] > 0 ){
                this->_substitute_matrix[i][j] = _nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[i][j] * match_score;
            }else{
                this->_substitute_matrix[i][j] = -(_nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[i][j] * mis_match_score);
            }
        }
    }
}

NeedlemanWunsch::~NeedlemanWunsch() {
    size_t i;
    for (i = 0; i <= _length_of_q; ++i) {
        delete[] this->_similarity_matrix[i];
        delete[] this->_track_matrix[i];
    }
    delete[] this->_similarity_matrix;
    delete[] this->_track_matrix;
    for (i = 0; i < _n; ++i) {
        delete[] this->_substitute_matrix[i];
    }
    delete[] this->_substitute_matrix;
}

// Calculating similarity matrix
void NeedlemanWunsch::calculate_similarity( NucleotideCodeSubstitutionMatrix& _nucleotideCodeSubstitutionMatrix ) {

    int match = 0, insert = 0, del = 0 , selected = 0, t;
    for (size_t i = 1; i < _length_of_q + 1; ++i) {
        for (size_t j = 1; j < _length_of_d + 1; ++j) {
            match = _similarity_matrix[i - 1][j - 1] + 
                _substitute_matrix[_nucleotideCodeSubstitutionMatrix.get_dna_acid_map(_dna_q[i - 1] ) ][_nucleotideCodeSubstitutionMatrix.get_dna_acid_map(_dna_d[j - 1]) ];
            insert = _similarity_matrix[i - 1][j] + this->_open_gap_penalty;
            del = _similarity_matrix[i][j - 1] + this->_open_gap_penalty;

            if ( _track_matrix[i - 1][j] == INSERTION || _track_matrix[i - 1][j] == SNPORINSERTION ||  _track_matrix[i - 1][j] == INSERTIONORDELETION
                || _track_matrix[i - 1][j] == SNPORINSERTIONORDELETION ){
                insert = _similarity_matrix[i - 1][j] + this->_extend_gap_penalty;
            }
            if ( _track_matrix[i][j - 1] == DELETION || _track_matrix[i][j - 1] == SNPORDELETION ||_track_matrix[i][j-1] == INSERTIONORDELETION
                || _track_matrix[i][j - 1] == SNPORINSERTIONORDELETION ){
                del = _similarity_matrix[i][j - 1] + this->_extend_gap_penalty;
            }
            if( del >insert && del==match  ){
                selected = del;
                _track_matrix[i][j] = SNPORDELETION;
            }else if( insert >del && insert == match  ){
                selected = match;
                _track_matrix[i][j] = SNPORINSERTION;
            }else if ( insert > match && insert > del){// prefer deletion
                int t = 1;
                while( i-t >=1 && (_track_matrix[i - t][j] == SNPORINSERTION || _track_matrix[i - t][j] == SNPORINSERTIONORDELETION || _track_matrix[i - t][j]==INSERTIONORDELETION ) ){
                    _track_matrix[i - t][j] = INSERTION;
                    ++t;
                }
                selected = insert;
                _track_matrix[i][j] = INSERTION;
            }else if( del > match && del > insert ){//prefer insertion, so that the INDELs could be put together
                t = 1;
                while( j-t >=1 && (_track_matrix[i][j-t] == SNPORDELETION || _track_matrix[i][j-t] == SNPORINSERTIONORDELETION || _track_matrix[i][j-t]==INSERTIONORDELETION) ){
                    _track_matrix[i][j-t] = DELETION;
                    ++t;
                }
                selected = del;
                _track_matrix[i][j] = DELETION;
            }else if (match > insert && match > del){
                t = 1;
                while( i-t >=1 && j-t>=1 && (_track_matrix[i-t][j-t] == SNPORINSERTION || _track_matrix[i-t][j-t] == SNPORINSERTIONORDELETION || _track_matrix[i-t][j-t]==SNPORDELETION ) ){
                    _track_matrix[i-t][j-t] = SNP;
                    ++t;
                }
                selected = match;
                _track_matrix[i][j] = SNP;
            }else if ( del >match && insert==del  ){
                selected = del;
                _track_matrix[i][j] = INSERTIONORDELETION;
            } else{
                selected = del;
                _track_matrix[i][j] = SNPORINSERTIONORDELETION;
            }
            _similarity_matrix[i][j] = selected;
        }
    }
}


void NeedlemanWunsch::get_optimal_alignment() {
    _alignment_q = "";
    _alignment_d = "";
    std::stack<char> SQ, SD;

    size_t k;
    size_t i = this->_length_of_q;
    size_t j = this->_length_of_d;
    int32_t highestScore = this->_similarity_matrix[i][j];
    for( k=_length_of_q; k >0; --k ){
        if( _similarity_matrix[k][j] > highestScore ){
            i = k;
            highestScore= _similarity_matrix[k][j];
        }
    }

    for( k=_length_of_d; k >0; --k ){
        if( _similarity_matrix[_length_of_q][k] > highestScore ){
            j = k;
            i = this->_length_of_q;
            highestScore= _similarity_matrix[_length_of_q][k];
        }
    }

    for( k=_length_of_q; k>i; --k ){
        SQ.push(_dna_q[k-1]);
        SD.push('-');
    }

    for( k=_length_of_d; k>j; --k ){
        SQ.push('-');
        SD.push(_dna_d[k-1]);
    }

    while (i > 0 || j > 0) {
        if (i == 0) {
            SQ.push('-');
            SD.push(_dna_d[j-1]);
            --j;
        } else if (j == 0) {
            SQ.push(_dna_q[i-1]);
            SD.push('-');
            --i;
        }else{
            if ( _track_matrix[i][j]==SNP || _track_matrix[i][j]==SNPORDELETION || _track_matrix[i][j]==SNPORINSERTION || _track_matrix[i][j]==SNPORINSERTIONORDELETION  ) {
                SQ.push(_dna_q[i - 1]);
                SD.push(_dna_d[j - 1]);
                --i;
                --j;
            } else if (_track_matrix[i][j]==INSERTION || _track_matrix[i][j]==INSERTIONORDELETION) {//insertion
                SQ.push(_dna_q[i - 1]);
                SD.push('-');
                --i;
            } else {
                SQ.push('-');
                SD.push(_dna_d[j - 1]);
                --j;
            }
        }
    }
    while (!SQ.empty()) {
        _alignment_q += SQ.top();
        _alignment_d += SD.top();
        SQ.pop();
        SD.pop();
    }
}

std::string NeedlemanWunsch::getAlignment_q(){
    return _alignment_q;
}
std::string NeedlemanWunsch::getAlignment_d(){
    return _alignment_d;
}
void NeedlemanWunsch::print_results() {
    std::cout << _alignment_q << std::endl;
    std::cout << _alignment_d << std::endl;
}
