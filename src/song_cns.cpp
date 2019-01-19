//
// Created by Baoxing song on 2018-12-28.
//

#include "song_cns.h"


SimilarBlocks::SimilarBlocks(const std::string & dna_q, const std::string & dna_d, const int32_t & match_score, const int32_t & mis_match_score,
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
        this->_similarity_matrix_E[i][j] = 0;
        this->_similarity_matrix_F[i][j] = -std::numeric_limits<int32_t >::infinity();
    }
    j=0;
    for (i = 0; i < (this->_length_of_q+1); ++i) {
        this->_similarity_matrix[i][j] = 0;
        this->_similarity_matrix_E[i][j] = -std::numeric_limits<int32_t >::infinity();;;
        this->_similarity_matrix_F[i][j] = 0;
    }

    // this matrix is for set different penalty for open gap and extend gap begin
    // and the track also changed to use this matrix

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

SimilarBlocks::~SimilarBlocks() {
    size_t i;
    for (i = 0; i <= _length_of_q; ++i) {
        delete[] this->_similarity_matrix[i];
    }
    delete[] this->_similarity_matrix;
    for (i = 0; i < _n; ++i) {
        delete[] this->_substitute_matrix[i];
    }
    delete[] this->_substitute_matrix;
}

// Calculating similarity matrix
void SimilarBlocks::calculate_similarity( NucleotideCodeSubstitutionMatrix& _nucleotideCodeSubstitutionMatrix ) {
    int match = 0;
    for (size_t i = 1; i < _length_of_q + 1; ++i) {
        for (size_t j = 1; j < _length_of_d + 1; ++j) {
            match = _similarity_matrix[i - 1][j - 1] +
                    _substitute_matrix[_nucleotideCodeSubstitutionMatrix.get_dna_acid_map(_dna_q[i - 1] ) ][_nucleotideCodeSubstitutionMatrix.get_dna_acid_map(_dna_d[j - 1]) ];

            _similarity_matrix_F[i][j] = _similarity_matrix[i-1][j] + this->_open_gap_penalty > _similarity_matrix_F[i-1][j]+this->_extend_gap_penalty ? _similarity_matrix[i-1][j] + this->_open_gap_penalty : _similarity_matrix_F[i-1][j]+this->_extend_gap_penalty;
            _similarity_matrix_E[i][j] = _similarity_matrix[i][j-1] + this->_open_gap_penalty > _similarity_matrix_E[i][j-1]+this->_extend_gap_penalty ? _similarity_matrix[i][j-1] + this->_open_gap_penalty : _similarity_matrix_E[i][j-1]+this->_extend_gap_penalty;

            if( match >= _similarity_matrix_E[i][j] && match >= _similarity_matrix_F[i][j] ){
                _similarity_matrix[i][j] = match;
            }else if ( _similarity_matrix_E[i][j] >= match && _similarity_matrix_E[i][j] >= _similarity_matrix_F[i][j] ){
                _similarity_matrix[i][j] = _similarity_matrix_E[i][j];
            }else{
                _similarity_matrix[i][j] = _similarity_matrix_F[i][j];
            }
        }
    }
}


void SimilarBlocks::get_optimal_alignment() {
    /*
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
     */
}

std::string SimilarBlocks::getAlignment_q(){
    return _alignment_q;
}
std::string SimilarBlocks::getAlignment_d(){
    return _alignment_d;
}
void SimilarBlocks::print_results() {
    std::cout << _alignment_q << std::endl;
    std::cout << _alignment_d << std::endl;
}
