/*
 * =====================================================================================
 *
 *       Filename:  alignNeedlemanForTranscript.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/02/2017 15:11:39
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Baoxing Song (songbx.me), songbaoxing168@163.com
 *
 * =====================================================================================
 */
#include <stdlib.h>
/*************************************************************************




 ************************************************************************/

#include <iostream>
#include <map>
#include "alignNeedlemanForTranscript.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include "../util/parameters.h"
#include <algorithm>
SpliceSitePosition::SpliceSitePosition(int donorsSpliceSitePosition, int acceptorSpliceSitePosition) {
    this->_donorSpliceSitePosition=donorsSpliceSitePosition;
    this->_acceptorSpliceSitePosition=acceptorSpliceSitePosition;
}
int SpliceSitePosition::getDonorSpliceSitePosition(){
    return _donorSpliceSitePosition;
}

int SpliceSitePosition::getAcceptorSpliceSitePosition(){
    return _acceptorSpliceSitePosition;
}
NeedlemanWunschForTranscript::NeedlemanWunschForTranscript(std::string& dna_d, std::string& dna_q, int & startCodonPosition, int & stopCodonPosition,
                                                           std::vector<SpliceSitePosition>& splitSitePositions, std::map<std::string, std::string>& parameters, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix){
    this->_dna_d = dna_d;
    this->_dna_q = dna_q;
    this->_exon_match_score = stoi(get_parameters("alignmentExonMatchP", parameters));
    this->_exon_mismatch_score = stoi(get_parameters("alignmentExonMismatchP", parameters));
    this->_exon_open_gap_penalty = stoi(get_parameters("alignmentExonOpenGapP", parameters));
    this->_exon_extend_gap_penalty = stoi(get_parameters("alignmentExonExtendGapP", parameters));

    this->_intron_match_score = stoi(get_parameters("alignmentIntronMatchP", parameters));
    this->_intron_mismatch_score  = stoi(get_parameters("alignmentIntronMismatchP", parameters));
    this->_intron_open_gap_penalty = stoi(get_parameters("alignmentIntronOpenGapP", parameters));
    this->_intron_extend_gap_penalty = stoi(get_parameters("alignmentIntronExtendGapP", parameters));

    this->_splice_sites_match_score = stoi(get_parameters("alignmentSpliceSitesMatchP", parameters));
    this->_splice_sites_mismatch_score = stoi(get_parameters("alignmentSpliceSitesMismatchP", parameters));
    this->_splice_sites_open_gap_penalty = stoi(get_parameters("alignmentSpliceSitesOpenGapP", parameters));
    this->_splice_sites_extend_gap_penalty = stoi(get_parameters("alignmentSpliceSitesExtendGapP", parameters));

    this->_start_stop_codon_match_score = stoi(get_parameters("alignmentStartStopCodonMatchP", parameters));
    this->_start_stop_codon_mismatch_score = stoi(get_parameters("alignmentStartStopCodonMismatchP", parameters));
    this->_start_stop_codon_open_gap_penalty = stoi(get_parameters("alignmentStartStopCodonOpenGapP", parameters));
    this->_start_stop_codon_extend_gap_penalty = stoi(get_parameters("alignmentStartStopCodonExtendGapP", parameters));
    this->_startCodonPosition = startCodonPosition;
    this->_stopCodonPosition = stopCodonPosition;
    this->_spliceSitePositions = splitSitePositions;
    this->_length_of_q = this->_dna_q.length();
    this->_length_of_d = this->_dna_d.length();
    this->_similarity_matrix = new int32_t*[this->_length_of_q + 1];
    for (int i = 0; i < (this->_length_of_q + 1); i++) {
        this->_similarity_matrix[i] = new int32_t[this->_length_of_d + 1];
    }
    for (int i = 0; i<= _length_of_q; i++) {
        for (int j = 0; j<= 1; j++) {
            _similarity_matrix[i][j] = 0;
        }
    }
    for (int j = 0; j<= _length_of_d; j++){
        for (int i = 0; i<= 1; i++) {
            _similarity_matrix[i][j] = 0;
        }
    }
    // this matrix is for set different penalty for open gap and extend gap begin
    // 0 for match, 1 for deletion, 2 for insertation
    // and the track also changed to use this matrix
    this->_track_matrix = new VARIANTCATEGORY*[this->_length_of_q + 1];
    for (int i = 0; i < (this->_length_of_q + 1); i++) {
        _track_matrix[i] = new VARIANTCATEGORY[this->_length_of_d + 1];
    }

    for (int i = 0; i <= 1; i++) {
        for (int j = 0; j <= _length_of_d; j++) {
            _track_matrix[i][j] = SNPORINSERTIONORDELETION;
        }
    }
    for (int j = 0; j <= 1; j++){
        for (int i = 0; i <= _length_of_q; i++) {
            _track_matrix[i][j] = SNPORINSERTIONORDELETION;
        }
    }

    this->calculate_similarity(nucleotideCodeSubstitutionMatrix);
    this->get_optimal_alignment();
}

NeedlemanWunschForTranscript::~NeedlemanWunschForTranscript() {
    for (int i = 0; i <= _length_of_q; i++) {
        delete[] this->_similarity_matrix[i];
        delete[] this->_track_matrix[i];
    }
    delete[] this->_similarity_matrix;
    delete[] this->_track_matrix;
}

ELEMENTS NeedlemanWunschForTranscript::checkElements( int & position ){
    if( position < this->_startCodonPosition || position > this->_stopCodonPosition+2  ){
        return INTRON;
    }
    for(std::vector<SpliceSitePosition>::size_type i = 0; i != this->_spliceSitePositions.size(); i++ ){
    	SpliceSitePosition spliceSite = this->_spliceSitePositions[i];
        if(position > (spliceSite.getDonorSpliceSitePosition()+1) && position < (spliceSite.getAcceptorSpliceSitePosition()-1) ){
            return INTRON;
        }else if( position== spliceSite.getAcceptorSpliceSitePosition() || position == spliceSite.getAcceptorSpliceSitePosition()-1){
            return SPLICEACCEPTOR;
        }else if( position== spliceSite.getDonorSpliceSitePosition() || position == spliceSite.getDonorSpliceSitePosition()+1  ){
            return SPLICEDONOR;
        }
    }
    if (position == this->_startCodonPosition || position == this->_startCodonPosition+1 || position == this->_startCodonPosition+2 ){
        return START;
    }
    if (position == this->_stopCodonPosition || position == this->_stopCodonPosition+1 || position == this->_stopCodonPosition+2 ){
        return STOP;
    }
    return EXON;
}

// Calculating similarity matrix
void NeedlemanWunschForTranscript::calculate_similarity( NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix){

    int i;
    int match = 0, insert = 0, del = 0;
    for (int j=1; j < _length_of_d + 1; ++j){
        if( INTRON == this->checkElements(j) ) {
            for (i = 1; i < _length_of_q + 1; ++i) {
                match = _similarity_matrix[i - 1][j - 1] + nucleotideCodeSubstitutionMatrix.get_intron_subsitition_matrix()
                [nucleotideCodeSubstitutionMatrix.get_dna_acid_map(_dna_q[i - 1])]
                [nucleotideCodeSubstitutionMatrix.get_dna_acid_map(_dna_d[j - 1])];
                if (i==1 || _track_matrix[i - 1][j] == INSERTION || _track_matrix[i - 1][j] == SNPORINSERTION
                            || _track_matrix[i - 1][j] == SNPORINSERTIONORDELETION || _track_matrix[i - 1][j] == INSERTIONORDELETION ) { //deletion
                    insert = _similarity_matrix[i - 1][j] + this->_intron_extend_gap_penalty;
                } else {
                    insert = _similarity_matrix[i - 1][j] + this->_intron_open_gap_penalty;
                }
                if (j==1 || _track_matrix[i][j - 1] == DELETION || _track_matrix[i][j - 1] == SNPORDELETION ||
                            _track_matrix[i][j - 1] == SNPORINSERTIONORDELETION || _track_matrix[i][j-1] == INSERTIONORDELETION ) { //insertion
                    del = _similarity_matrix[i][j - 1] + this->_intron_extend_gap_penalty;
                } else {
                    del = _similarity_matrix[i][j - 1] + this->_intron_open_gap_penalty;
                }
                setScore( match, insert, del, i, j);
            }
        }else if( EXON == this->checkElements(j)  ){
            for (i = 1; i < _length_of_q + 1; ++i) {
                match = _similarity_matrix[i - 1][j - 1] + nucleotideCodeSubstitutionMatrix.get_exon_subsitition_matrix()
                [nucleotideCodeSubstitutionMatrix.get_dna_acid_map(_dna_q[i - 1])]
                [nucleotideCodeSubstitutionMatrix.get_dna_acid_map(_dna_d[j - 1])];
                if (i==1 || _track_matrix[i - 1][j] == INSERTION || _track_matrix[i - 1][j] == SNPORINSERTION
                            || _track_matrix[i - 1][j] == SNPORINSERTIONORDELETION || _track_matrix[i - 1][j] == INSERTIONORDELETION ) { //deletion
                    insert = _similarity_matrix[i - 1][j] + this->_exon_extend_gap_penalty;
                } else {
                    insert = _similarity_matrix[i - 1][j] + this->_exon_open_gap_penalty;
                }
                if (j==1 || _track_matrix[i][j - 1] == DELETION || _track_matrix[i][j - 1] == SNPORDELETION ||
                            _track_matrix[i][j - 1] == SNPORINSERTIONORDELETION || _track_matrix[i][j-1] == INSERTIONORDELETION ) { //insertion
                    del = _similarity_matrix[i][j - 1] + this->_exon_extend_gap_penalty;
                } else {
                    del = _similarity_matrix[i][j - 1] + this->_exon_open_gap_penalty;
                }
                setScore( match, insert, del, i, j);
            }
        }else if( START == this->checkElements(j) || STOP == this->checkElements(j)){
            for ( i = 1; i < _length_of_q + 1; ++i) {
                match = _similarity_matrix[i - 1][j - 1] + nucleotideCodeSubstitutionMatrix.get_start_stop_codon_subsitition_matrix()
                [nucleotideCodeSubstitutionMatrix.get_dna_acid_map(_dna_q[i - 1])]
                [nucleotideCodeSubstitutionMatrix.get_dna_acid_map(_dna_d[j - 1])];

                if (i==1 || _track_matrix[i - 1][j] == INSERTION || _track_matrix[i - 1][j] == SNPORINSERTION
                            || _track_matrix[i - 1][j] == SNPORINSERTIONORDELETION || _track_matrix[i - 1][j] == INSERTIONORDELETION ) { //deletion
                    insert = _similarity_matrix[i - 1][j] + this->_start_stop_codon_extend_gap_penalty;
                } else {
                    insert = _similarity_matrix[i - 1][j] + this->_start_stop_codon_open_gap_penalty;
                }
                if (j==1 || _track_matrix[i][j - 1] == DELETION || _track_matrix[i][j - 1] == SNPORDELETION ||
                            _track_matrix[i][j - 1] == SNPORINSERTIONORDELETION || _track_matrix[i][j-1] == INSERTIONORDELETION ) { //insertion
                    del = _similarity_matrix[i][j - 1] + this->_start_stop_codon_extend_gap_penalty;
                } else {
                    del = _similarity_matrix[i][j - 1] + this->_start_stop_codon_open_gap_penalty;
                }
                setScore( match, insert, del, i, j);
/*
                // test for new start and stop codon score method start
                if (j == _startCodonPosition + 2 && i > 2) {
                    std::stringstream targetthree;
                    targetthree << _dna_d[j - 3];
                    targetthree << _dna_d[j - 2];
                    targetthree << _dna_d[j - 1];
                    std::string startThree = targetthree.str();
                    if (nucleotideCodeSubstitutionMatrix.getPossibleStartCodons().find(startThree) !=
                        nucleotideCodeSubstitutionMatrix.getPossibleStartCodons().end()) {
                        _similarity_matrix[i - 2][j - 2] =
                                _similarity_matrix[i - 3][j - 3] + nucleotideCodeSubstitutionMatrix.get_start_stop_codon_subsitition_matrix()[0][0];
                        _track_matrix[i - 2][j - 2] = SNP;
                        _similarity_matrix[i - 1][j - 1] =
                                _similarity_matrix[i - 2][j - 2] + nucleotideCodeSubstitutionMatrix.get_start_stop_codon_subsitition_matrix()[0][0];
                        _track_matrix[i - 1][j - 1] = SNP;
                        _similarity_matrix[i][j] =
                                _similarity_matrix[i - 1][j - 1] + nucleotideCodeSubstitutionMatrix.get_start_stop_codon_subsitition_matrix()[0][0];
                        _track_matrix[i][j] = SNP;
                    }
                }

                if (j == _stopCodonPosition + 2 && i > 2) {
                    std::stringstream targetthree;
                    targetthree << _dna_d[j - 3];
                    targetthree << _dna_d[j - 2];
                    targetthree << _dna_d[j - 1];
                    std::string stopThree = targetthree.str();
                    // the the query has been the end, and the snp is not higher than insertion, do not change
                    if (nucleotideCodeSubstitutionMatrix.getPossibleStopCodons().find(stopThree) !=
                        nucleotideCodeSubstitutionMatrix.getPossibleStopCodons().end()
                        && !(i == _length_of_q &&
                             (_similarity_matrix[i - 3][j - 3] + 3 * nucleotideCodeSubstitutionMatrix.get_start_stop_codon_subsitition_matrix()[0][0]) <
                             _similarity_matrix[i][j])) {
                        _similarity_matrix[i - 2][j - 2] =
                                _similarity_matrix[i - 3][j - 3] + nucleotideCodeSubstitutionMatrix.get_start_stop_codon_subsitition_matrix()[0][0];
                        _track_matrix[i - 2][j - 2] = SNP;
                        _similarity_matrix[i - 1][j - 1] =
                                _similarity_matrix[i - 2][j - 2] + nucleotideCodeSubstitutionMatrix.get_start_stop_codon_subsitition_matrix()[0][0];
                        _track_matrix[i - 1][j - 1] = SNP;
                        _similarity_matrix[i][j] =
                                _similarity_matrix[i - 1][j - 1] + nucleotideCodeSubstitutionMatrix.get_start_stop_codon_subsitition_matrix()[0][0];
                        _track_matrix[i][j] = SNP;
                    }
                }
                // test for new start and stop codon score method end
                */
            }
        }else if ( SPLICEDONOR == this->checkElements(j) || SPLICEACCEPTOR == this->checkElements(j)  ){
            for (i = 1; i < _length_of_q + 1; ++i) {
                match = _similarity_matrix[i - 1][j - 1] + nucleotideCodeSubstitutionMatrix.get_splice_sites_subsitition_matrix()
                [nucleotideCodeSubstitutionMatrix.get_dna_acid_map(_dna_q[i - 1] ) ]
                [nucleotideCodeSubstitutionMatrix.get_dna_acid_map(_dna_d[j - 1] ) ];
                if(i==1 ||  _track_matrix[i - 1][j] == INSERTION || _track_matrix[i - 1][j] == SNPORINSERTION
                            || _track_matrix[i - 1][j] == SNPORINSERTIONORDELETION || _track_matrix[i - 1][j] == INSERTIONORDELETION){ //deletion
                	insert = _similarity_matrix[i - 1][j] + this->_splice_sites_extend_gap_penalty;
            	}else{
                    insert = _similarity_matrix[i - 1][j] + this->_splice_sites_open_gap_penalty;
                }
            	if(j==1 || _track_matrix[i][j - 1] == DELETION || _track_matrix[i][j - 1] == SNPORDELETION ||
                           _track_matrix[i][j - 1] == SNPORINSERTIONORDELETION || _track_matrix[i][j-1] == INSERTIONORDELETION ){ //insertion
                	del = _similarity_matrix[i][j - 1] + this->_splice_sites_extend_gap_penalty;
            	}else{
                    del = _similarity_matrix[i][j - 1] + this->_splice_sites_open_gap_penalty;
                }
                setScore( match, insert, del, i, j);
            }
        }
    }
}

void NeedlemanWunschForTranscript::setScore(int32_t & match, int32_t & insert, int32_t & del, int & i, int & j){
    int32_t selected=0;
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
        int t = 1;
        while( j-t >=1 && (_track_matrix[i][j-t] == SNPORDELETION || _track_matrix[i][j-t] == SNPORINSERTIONORDELETION || _track_matrix[i][j-t]==INSERTIONORDELETION) ){
            _track_matrix[i][j-t] = DELETION;
            ++t;
        }
        selected = del;
        _track_matrix[i][j] = DELETION;
    }else if (match > insert && match > del){
        int t = 1;
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

//// Trace back step.
void NeedlemanWunschForTranscript::get_optimal_alignment() {
    _alignment_q = "";
    _alignment_d = "";
    std::stack<char> SQ, SD;

    int32_t k;
    size_t i = this->_length_of_q;
    size_t j = this->_length_of_d;
    int32_t highestScore = this->_similarity_matrix[i][j];
    for( k=_length_of_q; k >0; --k ){
        if( _similarity_matrix[k][j] > highestScore ){
            i = k;
            highestScore= _similarity_matrix[k][j];
        }
    }
    for( k=_length_of_q; k>i; --k ){
        SQ.push(_dna_q[k-1]);
        SD.push('-');
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
            } else if (_track_matrix[i][j]==DELETION || _track_matrix[i][j]==INSERTIONORDELETION) {// Going to S(i, j-1) //deletion
                SQ.push('-');
                SD.push(_dna_d[j - 1]);
                --j;
            } else {        //insertion
                SQ.push(_dna_q[i - 1]);
                SD.push('-');
                --i;
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

std::string NeedlemanWunschForTranscript::getAlignment_q(){
    return _alignment_q;
}
std::string NeedlemanWunschForTranscript::getAlignment_d(){
    return _alignment_d;
}
void NeedlemanWunschForTranscript::print_results() {
    std::cout << _alignment_q << std::endl;
    std::cout << _alignment_d << std::endl;
}
