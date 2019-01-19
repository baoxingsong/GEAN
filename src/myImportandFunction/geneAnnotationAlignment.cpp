//
// Created by Baoxing song on 09.10.18.
//

#include "geneAnnotationAlignment.h"

int32_t globalAlignment( std::vector<Gene> & g1s, std::vector<Gene> & g2s, std::vector<std::string> & _alignment_q,
    std::vector<std::string> & _alignment_d){
    int32_t ** _similarity_matrix = new int32_t *[g1s.size() + 1];
    int32_t  i, j;
    for (i = 0; i < (g1s.size() + 1); ++i) {
        _similarity_matrix[i] = new int32_t [g2s.size() + 1];
    }
    i=0;
    for ( j = 0; j <= g2s.size(); ++j) {
        _similarity_matrix[i][j] = -j;
    }
    j=0;
    for (i = 0; i < (g1s.size()+1); ++i) {
        _similarity_matrix[i][j] = -i;
    }

    // this matrix is for set different penalty for open gap and extend gap begin
    // and the track also changed to use this matrix
    VARIANTCATEGORY** _track_matrix = new VARIANTCATEGORY *[g1s.size() + 1];
    for (i = 0; i < (g1s.size() + 1); ++i) {
        _track_matrix[i] = new VARIANTCATEGORY[g2s.size() + 1];
    }

    for ( i = 0; i <= g1s.size(); ++i) {
        for ( j = 0; j <= g2s.size(); ++j) {
            _track_matrix[i][j] = SNPORINSERTIONORDELETION;
        }
    }
    i=0;
    for ( j = 1; j <= g2s.size(); ++j) {
        _track_matrix[i][j] = DELETION;
    }
    j=0;
    for (i = 1; i < (g1s.size()+1); ++i) {
        _track_matrix[i][j] = INSERTION;
    }

    int32_t _open_gap_penalty = -2;
    int32_t _extend_gap_penalty = -1;
    int match = 0, insert = 0, del = 0 , selected = 0, t;
    for (size_t i = 1; i < g1s.size() + 1; ++i) {
        for (size_t j = 1; j < g2s.size() + 1; ++j) {
            if( g1s[i-1].getName().compare(  g2s[j-1].getName() ) ==0 ){
                match = _similarity_matrix[i - 1][j - 1] +1;
            }else{
                match = _similarity_matrix[i - 1][j - 1] -1;
            }

            insert = _similarity_matrix[i - 1][j] + _open_gap_penalty;
            del = _similarity_matrix[i][j - 1] + _open_gap_penalty;

            if ( _track_matrix[i - 1][j] == INSERTION || _track_matrix[i - 1][j] == SNPORINSERTION ||  _track_matrix[i - 1][j] == INSERTIONORDELETION
                 || _track_matrix[i - 1][j] == SNPORINSERTIONORDELETION ){
                insert = _similarity_matrix[i - 1][j] + _extend_gap_penalty;
            }
            if ( _track_matrix[i][j - 1] == DELETION || _track_matrix[i][j - 1] == SNPORDELETION ||_track_matrix[i][j-1] == INSERTIONORDELETION
                 || _track_matrix[i][j - 1] == SNPORINSERTIONORDELETION ){
                del = _similarity_matrix[i][j - 1] + _extend_gap_penalty;
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
            } else {
                selected = del;
                _track_matrix[i][j] = SNPORINSERTIONORDELETION;
            }
            _similarity_matrix[i][j] = selected;
        }
    }

    std::stack<std::string> SQ, SD;

    size_t k;
    i = g1s.size();
    j = g2s.size();

    while (i > 0 || j > 0) {
        if (i == 0) {
            SQ.push("-");
            SD.push(g2s[j-1].getName());
            --j;
        } else if (j == 0) {
            SQ.push(g1s[i-1].getName());
            SD.push("-");
            --i;
        }else{
            if ( _track_matrix[i][j]==SNP || _track_matrix[i][j]==SNPORDELETION || _track_matrix[i][j]==SNPORINSERTION || _track_matrix[i][j]==SNPORINSERTIONORDELETION  ) {
                SQ.push(g1s[i - 1].getName());
                SD.push(g2s[j - 1].getName());
                --i;
                --j;
            } else if (_track_matrix[i][j]==INSERTION || _track_matrix[i][j]==INSERTIONORDELETION) {//insertion
                SQ.push(g1s[i - 1].getName());
                SD.push("-");
                --i;
            } else {
                SQ.push("-");
                SD.push(g2s[j - 1].getName());
                --j;
            }
        }
    }
    while (!SQ.empty()) {
        _alignment_q.push_back(SQ.top() );
        _alignment_d.push_back(SD.top());
        SQ.pop();
        SD.pop();
    }
    return _similarity_matrix[g1s.size()][g2s.size()];
}
