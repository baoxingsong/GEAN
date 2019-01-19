/*
 * =====================================================================================
 *
 *       Filename:  alignNeedlemanForTranscript_simd.cpp
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  03/August/2018 22:11:39
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Baoxing Song (songbx.me), songbaoxing168@163.com
 *   Organization:  The department of life science, Qiannan Normal University for Nationalities
 *
 * =====================================================================================
 */

/*************************************************************************




 ************************************************************************/

#include <iostream>
#include <map>
#include "alignNeedlemanForTranscript_simd_avx2int32.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include "../util/parameters.h"
#include <algorithm>

alignNeedlemanForTranscript_simd_avx2int32::alignNeedlemanForTranscript_simd_avx2int32(std::string& dna_d,
        std::string& dna_q, int & startCodonPosition, int & stopCodonPosition,
        std::vector<SpliceSitePosition>& splitSitePositions,
        std::map<std::string, std::string>& parameters, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix){
    //std::cout << "NeedlemanWunschForTranscript begin" << std::endl;
    this->_dna_q = dna_q;
    this->_dna_d = dna_d;
    this->_exon_match_score = stoi(get_parameters("alignmentExonMatchP", parameters));
    this->_exon_mismatch_score = stoi(get_parameters("alignmentExonMismatchP", parameters));
    this->_exon_open_gap_penalty = stoi(get_parameters("alignmentExonOpenGapP", parameters))*36;
    this->_exon_extend_gap_penalty = stoi(get_parameters("alignmentExonExtendGapP", parameters))*36;

    this->_intron_match_score = stoi(get_parameters("alignmentIntronMatchP", parameters));
    this->_intron_mismatch_score  = stoi(get_parameters("alignmentIntronMismatchP", parameters));
    this->_intron_open_gap_penalty = stoi(get_parameters("alignmentIntronOpenGapP", parameters))*36;
    this->_intron_extend_gap_penalty = stoi(get_parameters("alignmentIntronExtendGapP", parameters))*36;

    this->_splice_sites_match_score = stoi(get_parameters("alignmentSpliceSitesMatchP", parameters));
    this->_splice_sites_mismatch_score = stoi(get_parameters("alignmentSpliceSitesMismatchP", parameters));
    this->_splice_sites_open_gap_penalty = stoi(get_parameters("alignmentSpliceSitesOpenGapP", parameters))*36;
    this->_splice_sites_extend_gap_penalty = stoi(get_parameters("alignmentSpliceSitesExtendGapP", parameters))*36;

    this->_start_stop_codon_match_score = stoi(get_parameters("alignmentStartStopCodonMatchP", parameters));
    this->_start_stop_codon_mismatch_score = stoi(get_parameters("alignmentStartStopCodonMismatchP", parameters));
    this->_start_stop_codon_open_gap_penalty = stoi(get_parameters("alignmentStartStopCodonOpenGapP", parameters))*36;
    this->_start_stop_codon_extend_gap_penalty = stoi(get_parameters("alignmentStartStopCodonExtendGapP", parameters))*36;
    this->_startCodonPosition = startCodonPosition;
    this->_stopCodonPosition = stopCodonPosition;
    this->_spliceSitePositions = splitSitePositions;

    this->_length_of_q = this->_dna_q.length();
    this->_length_of_d = this->_dna_d.length();
    this->bias = 0;
    this->_n = 16;
    this->_segLen = ( _length_of_q+ 7)/8;

    //check potential overflow problem begin
    int32_t maxLength = _length_of_d;
    if( maxLength <  _length_of_q){
        maxLength = _length_of_q;
    }
    int32_t maxScore = _exon_match_score;
    if( maxScore < _intron_match_score ){
        maxScore = _intron_match_score;
    }
    if( maxScore < _splice_sites_match_score ){
        maxScore = _splice_sites_match_score;
    }
    if( maxScore < _start_stop_codon_match_score ){
        maxScore = _start_stop_codon_match_score;
    }
    maxScore = 36 * maxScore;
    if( (maxLength * maxScore) > pow(2,31)){
        std::cerr << "the sequence is too long, could not be optimized with this pipeline" << std::endl;
        return;
    }
    int32_t minimumPenalty = _exon_open_gap_penalty;
    if( minimumPenalty > _exon_mismatch_score*36 ){
        minimumPenalty = _exon_mismatch_score*36;
    }
    if( minimumPenalty > _intron_open_gap_penalty ){
        minimumPenalty = _intron_open_gap_penalty;
    }
    if( minimumPenalty > _intron_mismatch_score*36 ){
        minimumPenalty = _intron_mismatch_score*36;
    }
    if( minimumPenalty > _splice_sites_open_gap_penalty ){
        minimumPenalty = _splice_sites_open_gap_penalty;
    }
    if( minimumPenalty > _splice_sites_mismatch_score*36 ){
        minimumPenalty = _splice_sites_mismatch_score*36;
    }

    if( minimumPenalty > _start_stop_codon_open_gap_penalty ){
        minimumPenalty = _start_stop_codon_open_gap_penalty;
    }
    if( minimumPenalty > _start_stop_codon_mismatch_score*36 ){
        minimumPenalty = _start_stop_codon_mismatch_score*36;
    }
    minimumPenalty = 36 * minimumPenalty;
    if( (maxLength * minimumPenalty) < -pow(2,31)){
        std::cerr << "the sequence is too long, could not be optimized with this pipeline" << std::endl;
        return;
    } //check potential overflow problem end

    if( _exon_open_gap_penalty > _exon_extend_gap_penalty ){//warning message about score strategy begin
        std::cerr << "it is not reasonable to have open gap penalty larger than extend gap penalty" << std::endl;
    }
    if( _intron_open_gap_penalty > _intron_extend_gap_penalty ){
        std::cerr << "it is not reasonable to have open gap penalty larger than extend gap penalty" << std::endl;
    }
    if( _splice_sites_open_gap_penalty > _splice_sites_extend_gap_penalty ){
        std::cerr << "it is not reasonable to have open gap penalty larger than extend gap penalty" << std::endl;
    }
    if( _start_stop_codon_open_gap_penalty > _start_stop_codon_extend_gap_penalty ){
        std::cerr << "it is not reasonable to have open gap penalty larger than extend gap penalty" << std::endl;
    }//warning message about score strategy end

    int32_t i, j;
    this->_similarity_matrix = new int32_t*[this->_length_of_q + 1];
    for (i = 0; i < (this->_length_of_q + 1); ++i) {
        this->_similarity_matrix[i] = new int32_t[this->_length_of_d + 1];
    }
    j=0;
    for ( i = 0; i<= _length_of_q; ++i) {
        _similarity_matrix[i][j] = 0;
    }
    i=0;
    for (j = 0; j<= _length_of_d; ++j){
        _similarity_matrix[i][j] = 0;
    }

    this->_track_match = new int32_t*[this->_length_of_q + 1];
    this->_track_del = new int32_t*[this->_length_of_q + 1];
    this->_track_ins = new int32_t*[this->_length_of_q + 1];
    for ( i = 0; i < (this->_length_of_q + 1); i++) {
        _track_match[i] = new int32_t[this->_length_of_d + 1];
        _track_del[i] = new int32_t[this->_length_of_d + 1];
        _track_ins[i] = new int32_t[this->_length_of_d + 1];
    }
    for ( i = 0; i <= 1; i++) {
        for ( j = 0; j <= _length_of_d; j++) {
            _track_match[i][j] = -1;
            _track_del[i][j] = -1;
            _track_ins[i][j] = -1;
        }
    }
    for ( j = 0; j <= 1; j++){
        for ( i = 0; i <= _length_of_q; i++) {
            _track_match[i][j] = -1;
            _track_del[i][j] = -1;
            _track_ins[i][j] = -1;
        }
    }
    std::map<std::string, __m256i*> vProfile = query_profile_avx2_byte(nucleotideCodeSubstitutionMatrix);
    this->calculate_similarity(vProfile, nucleotideCodeSubstitutionMatrix);
    this->get_optimal_alignment(nucleotideCodeSubstitutionMatrix);
}

alignNeedlemanForTranscript_simd_avx2int32::~alignNeedlemanForTranscript_simd_avx2int32() {
    for (int i = 0; i < (_length_of_q + 1); i++) {
        delete[] this->_similarity_matrix[i];
        delete[] this->_track_match[i];
        delete[] this->_track_del[i];
        delete[] this->_track_ins[i];
    }
    delete[] this->_similarity_matrix;
    delete[] this->_track_match;
    delete[] this->_track_del;
    delete[] this->_track_ins;
}

ELEMENTS alignNeedlemanForTranscript_simd_avx2int32::checkElements( int position ){
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

std::map<std::string, __m256i*> alignNeedlemanForTranscript_simd_avx2int32::query_profile_avx2_byte(
        NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix){ // dna_q is the query, (the non-reference sequence, with gene structure unknown)
    __m256i* exon_profile = (__m256i*)malloc(this->_n * this->_segLen * sizeof(__m256i));
    __m256i* intron_profile = (__m256i*)malloc(this->_n * this->_segLen * sizeof(__m256i));
    __m256i* splice_sites_profile = (__m256i*)malloc(this->_n * this->_segLen * sizeof(__m256i));
    __m256i* start_stop_profile = (__m256i*)malloc(this->_n * this->_segLen * sizeof(__m256i));

    int32_t a, i, j, k;
    for( a =0; a<_n; ++a ){
        //vProfile[a] = new __m256i[segLen];
        int32_t* exon_s = (int32_t*)(exon_profile+a*this->_segLen);
        int32_t* intron_s = (int32_t*)(intron_profile+a*this->_segLen);
        int32_t* splice_sites_s = (int32_t*)(splice_sites_profile+a*this->_segLen);
        int32_t* start_stop_s = (int32_t*)(start_stop_profile+a*this->_segLen);

        for( i=0; i<this->_segLen; ++i ){
            j = i;
            for( k=0; k<8; ++k ){
                *exon_s++ = j>= _length_of_q ? bias : nucleotideCodeSubstitutionMatrix.get_exon_subsitition_matrix()[a][nucleotideCodeSubstitutionMatrix.get_dna_acid_map()[_dna_q[j] ]];
                *intron_s++ = j>= _length_of_q ? bias : nucleotideCodeSubstitutionMatrix.get_intron_subsitition_matrix()[a][nucleotideCodeSubstitutionMatrix.get_dna_acid_map()[_dna_q[j] ]];
                *splice_sites_s++ = j>= _length_of_q ? bias : nucleotideCodeSubstitutionMatrix.get_splice_sites_subsitition_matrix()[a][nucleotideCodeSubstitutionMatrix.get_dna_acid_map()[_dna_q[j] ]];
                *start_stop_s++ = j>= _length_of_q ? bias : nucleotideCodeSubstitutionMatrix.get_start_stop_codon_subsitition_matrix()[a][nucleotideCodeSubstitutionMatrix.get_dna_acid_map()[_dna_q[j] ]];
                j += this->_segLen;
            }
        }
    }
    std::map<std::string, __m256i*> vProfile;
    vProfile["exon"] = exon_profile;
    vProfile["intron"] = intron_profile;
    vProfile["splice_sites"] = splice_sites_profile;
    vProfile["start_stop"] = start_stop_profile;
    return vProfile;
}

//
//// Calculating similarity matrix
//void calculate_similarity_it_works( std::map<std::string, __m256i*> & vProfile, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix){
//
//    const __m256i vZero = _mm256_set1_epi32(0);
//    const __m256i vMinuteOne = _mm256_set1_epi32(-1);
//    const __m256i vOne = _mm256_set1_epi32(-1);
//
//    __m256i* pvHStore = (__m256i*) calloc(this->_segLen, sizeof(__m256i));
//    __m256i* pvHLoad = (__m256i*) calloc(this->_segLen, sizeof(__m256i));
//    __m256i* trackMatch = (__m256i*) calloc(this->_segLen, sizeof(__m256i));
//    __m256i* pvE = (__m256i*) calloc(this->_segLen, sizeof(__m256i)); //vE is the score from left side
//    __m256i* pv;
//    __m256i* vP;
//    __m256i vcmp, vGapO, vGapE;
//
//    int32_t i, j, z;
//    int16_t * open_gap_penalty;
//    int16_t * extend_gap_penalty;
//
//       //for outer loop
//    __m256i vE, vH, vHsave;
//    double deletionScore;
//    int32_t* a;
//    int32_t* b;
//    int lastOne;
//    for (i = 0; (i != _length_of_d); ++i) {//important
//        if( INTRON == this->checkElements(i) ) {
//            open_gap_penalty = &_intron_open_gap_penalty;
//            extend_gap_penalty = &_intron_extend_gap_penalty;
//            vP = vProfile["intron"] + ((nucleotideCodeSubstitutionMatrix.get_dna_acid_map()[_dna_d[i] ])*(this->_segLen)); // Correct part of the vProfile
//        }else if( EXON == this->checkElements(i)  ){
//            open_gap_penalty = & _exon_open_gap_penalty;
//            extend_gap_penalty = & _exon_extend_gap_penalty;
//            vP = vProfile["exon"] + ((nucleotideCodeSubstitutionMatrix.get_dna_acid_map()[_dna_d[i] ])*(this->_segLen)); // Correct part of the vProfile
//        }else if( START == this->checkElements(i) || STOP == this->checkElements(i) ){
//            open_gap_penalty = & _start_stop_codon_open_gap_penalty;
//            extend_gap_penalty = & _start_stop_codon_extend_gap_penalty;
//            vP = vProfile["start_stop"] + ((nucleotideCodeSubstitutionMatrix.get_dna_acid_map()[_dna_d[i] ])*(this->_segLen)); // Correct part of the vProfile
//        }else if( SPLICEDONOR == this->checkElements(i) || SPLICEACCEPTOR == this->checkElements(i) ){
//            open_gap_penalty = & _splice_sites_open_gap_penalty;
//            extend_gap_penalty = & _splice_sites_extend_gap_penalty;
//            vP = vProfile["splice_sites"] + ((nucleotideCodeSubstitutionMatrix.get_dna_acid_map()[_dna_d[i] ])*(this->_segLen)); // Correct part of the vProfile
//        }
//
//
//        // INDEL begin vector
//        vGapO = _mm256_set1_epi32(*open_gap_penalty); // here I use epi32, because for the score matrix I used int32_t
//        // INDEL extension vector
//        vGapE = _mm256_set1_epi32(*extend_gap_penalty);
//
//        if( i==0 ){
//            for(j=0; j<_segLen; ++j){
//                _mm256_storeu_si256(pvE+i, vGapE);
//            }
//        }
//
//        vH = _mm256_loadu_si256(pvHStore+(this->_segLen-1));// the last one
//        vH = _mm256_alignr_epi8(vH, _mm256_permute2x128_si256(vH, vH, _MM_SHUFFLE(0, 0, 2, 0)), 16 - 4);
//
//        pv = pvHLoad;
//        pvHLoad = pvHStore;
//        pvHStore = pv;
//
//        // inner loop to process the query sequence
//        for (j = 0; j < this->_segLen; ++j) {
//            vH = _mm256_add_epi32(vH, _mm256_loadu_si256(vP+j)); // add the scoring profile to vH
//            vHsave = _mm256_add_epi32(vH, vZero);
//            // Get max from vH, vE
//            vE = _mm256_loadu_si256(pvE+j);
//            vH = _mm256_max_epi32(vH, vE);
//
//            vcmp = _mm256_cmpeq_epi32(vH, vHsave);
//            _mm256_storeu_si256(trackMatch + j, vcmp);
//
//            // Save vH values.
//            _mm256_storeu_si256(pvHStore + j, vH);
//
//            // Update vE value. for next database base pair  here should be get the exact score as the standard smith-waterman algorithm
//            vH = _mm256_add_epi32(vH, vGapO);
//            vE = _mm256_add_epi32(vE, vGapE);
//            vE = _mm256_max_epi32(vE, vH);
//            _mm256_storeu_si256(pvE+j, vE);
//
////            vE = _mm256_add_epi32( vH, _mm256_add_epi32(_mm256_mul_epi32(_mm256_mul_epi32(vcmp, vMinuteOne),vGapO), _mm256_mul_epi32(_mm256_add_epi32(vcmp, vOne),vGapE)) );
//            _mm256_storeu_si256(pvE+j, vE);
//            // Load the next vH.
//            vH = _mm256_loadu_si256(pvHLoad + j);
//        }
//
//        // put the score into the score matrix begin
//        for( j=0; j<this->_segLen; ++j ){
//            a = ( int32_t*)(pvHStore+j);
//            b = ( int32_t*)(trackMatch+j);
//            for( z=0; z<8; ++z ){
//                if( (j + z*this->_segLen) < _length_of_q ){
//                    _similarity_matrix[j+z*this->_segLen+1][i+1] = a[z];
//                    _track_matrix[j+z*this->_segLen+1][i+1] = b[z];
//                }
//            }
//        }
//
//        lastOne = 1; // insertion
//        for( j=0; j<_length_of_q; ++j ){
//            if( lastOne == 1 ){
//                deletionScore = _similarity_matrix[j][i+1] + (int32_t)(*extend_gap_penalty);
//                if( deletionScore > _similarity_matrix[j+1][i+1] ){
//                    _similarity_matrix[j+1][i+1] = deletionScore;
//                    _track_matrix[j+1][i+1] = 1;
//                }else if ( deletionScore == _similarity_matrix[j+1][i+1] && _track_matrix[j+1][i+1] ==0 ){
//                    _track_matrix[j+1][i+1] = 1;
//                }else{
//                    lastOne=0;//MATCH OR MISMATCH
//                }
//            }else{
//                deletionScore = _similarity_matrix[j][i+1]+(int32_t)(*open_gap_penalty);
//                if( deletionScore > _similarity_matrix[j+1][i+1] ){
//                    _similarity_matrix[j+1][i+1] = deletionScore;
//                    _track_matrix[j+1][i+1] = 1;
//                } else if ( deletionScore == _similarity_matrix[j+1][i+1] && _track_matrix[j+1][i+1] ==0 ){
//                    _track_matrix[j+1][i+1] = 1;
//                } else{
//                    lastOne=0;//MATCH OR MISMATCH
//                }
//            }
//        // put the score into the score matrix end
///*
//            // test for new start and stop codon score method start
//            if (i == _startCodonPosition + 2 && j > 2) {
//                std::stringstream targetthree;
//                targetthree << _dna_d[i - 3];
//                targetthree << _dna_d[i - 2];
//                targetthree << _dna_d[i - 1];
//                std::string startThree = targetthree.str();
//                if (nucleotideCodeSubstitutionMatrix.getPossibleStartCodons().find(startThree) !=
//                    nucleotideCodeSubstitutionMatrix.getPossibleStartCodons().end()) {
//                    _similarity_matrix[j - 2][i - 2] =
//                            _similarity_matrix[j - 3][i - 3] + nucleotideCodeSubstitutionMatrix.get_start_stop_codon_subsitition_matrix()[0][0];
//                    _track_matrix[j - 2][i - 2] = -1;
//                    _similarity_matrix[j - 1][i - 1] =
//                            _similarity_matrix[j - 2][i - 2] + nucleotideCodeSubstitutionMatrix.get_start_stop_codon_subsitition_matrix()[0][0];
//                    _track_matrix[j - 1][i - 1] = -1;
//                    _similarity_matrix[j][i] =
//                            _similarity_matrix[j - 1][i - 1] + nucleotideCodeSubstitutionMatrix.get_start_stop_codon_subsitition_matrix()[0][0];
//                    _track_matrix[j][i] = -1;
//                }
//            }
//
//            if (i == _stopCodonPosition + 2 && j > 2) {
//                std::stringstream targetthree;
//                targetthree << _dna_d[i - 3];
//                targetthree << _dna_d[i - 2];
//                targetthree << _dna_d[i - 1];
//                std::string stopThree = targetthree.str();
//                // the the query has been the end, and the snp is not higher than insertion, do not change
//                if (nucleotideCodeSubstitutionMatrix.getPossibleStopCodons().find(stopThree) !=
//                    nucleotideCodeSubstitutionMatrix.getPossibleStopCodons().end()
//                    && !(j == _length_of_q &&
//                         (_similarity_matrix[j - 3][i - 3] + 3 * nucleotideCodeSubstitutionMatrix.get_start_stop_codon_subsitition_matrix()[0][0]) <
//                         _similarity_matrix[j][i])) {
//                    _similarity_matrix[j - 2][i - 2] =
//                            _similarity_matrix[j - 3][i - 3] + nucleotideCodeSubstitutionMatrix.get_start_stop_codon_subsitition_matrix()[0][0];
//                    _track_matrix[j - 2][i - 2] = -1;
//                    _similarity_matrix[j - 1][i - 1] =
//                            _similarity_matrix[j - 2][i - 2] + nucleotideCodeSubstitutionMatrix.get_start_stop_codon_subsitition_matrix()[0][0];
//                    _track_matrix[j - 1][i - 1] = -1;
//                    _similarity_matrix[j][i] =
//                            _similarity_matrix[j - 1][i - 1] + nucleotideCodeSubstitutionMatrix.get_start_stop_codon_subsitition_matrix()[0][0];
//                    _track_matrix[j][i] = -1;
//                }
//            }        // test for new start and stop codon score method end*/
//        }
//
//        //updata vH begin
//        for( j=0; j<this->_segLen; j++ ){
//            a = (int32_t*)(pvHStore+j);
//            for( z=0; z<8; ++z ){
//                if( (j + z*this->_segLen) < _length_of_q ){
//                    a[z] = _similarity_matrix[j+z*this->_segLen+1][i+1];
//                }
//            }
//        }//updata vH end
//    }
//    free(pvE);
//    free(pvHStore);
//    free(pvHLoad);
//    free(trackMatch);
//
//    free(vProfile["exon"]);
//    free(vProfile["intron"]);
//    free(vProfile["splice_sites"]);
//    free(vProfile["start_stop"]);
//}
//




// Calculating similarity matrix
void alignNeedlemanForTranscript_simd_avx2int32::calculate_similarity( std::map<std::string, __m256i*> & vProfile,
        NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix){

    const __m256i vZero = _mm256_set1_epi32(0);
    const __m256i vMinuteOne = _mm256_set1_epi32(-1);
    const __m256i vOne = _mm256_set1_epi32(1);

    __m256i* pvHStore = (__m256i*) calloc(this->_segLen, sizeof(__m256i));
    __m256i* pvHLoad = (__m256i*) calloc(this->_segLen, sizeof(__m256i));
    __m256i* trackMatch = (__m256i*) calloc(this->_segLen, sizeof(__m256i));
    __m256i* trackDel = (__m256i*) calloc(this->_segLen, sizeof(__m256i));
    __m256i* pv;
    __m256i* vP;
    __m256i vcmp, vGapO, vGapE;

    int32_t i, j, z, t;
    int32_t * open_gap_penalty;
    int32_t * extend_gap_penalty;

    //for outer loop
    __m256i vE, vH, vHsave;
    int32_t deletionScore;
    int32_t* a;
    int32_t* b;
    int32_t* c;
    int lastOne;
    for (i = 0; i<_length_of_d; ++i) {//important
        if( INTRON == this->checkElements(i+1) ) {
            open_gap_penalty = &_intron_open_gap_penalty;
            extend_gap_penalty = &_intron_extend_gap_penalty;
            vP = vProfile["intron"] + ((nucleotideCodeSubstitutionMatrix.get_dna_acid_map()[_dna_d[i] ])*(this->_segLen)); // Correct part of the vProfile
        }else if( EXON == this->checkElements(i+1)  ){
            open_gap_penalty = & _exon_open_gap_penalty;
            extend_gap_penalty = & _exon_extend_gap_penalty;
            vP = vProfile["exon"] + ((nucleotideCodeSubstitutionMatrix.get_dna_acid_map()[_dna_d[i] ])*(this->_segLen)); // Correct part of the vProfile
        }else if( START == this->checkElements(i+1) || STOP == this->checkElements(i+1) ){
            open_gap_penalty = & _start_stop_codon_open_gap_penalty;
            extend_gap_penalty = & _start_stop_codon_extend_gap_penalty;
            vP = vProfile["start_stop"] + ((nucleotideCodeSubstitutionMatrix.get_dna_acid_map()[_dna_d[i] ])*(this->_segLen)); // Correct part of the vProfile
        }else if( SPLICEDONOR == this->checkElements(i+1) || SPLICEACCEPTOR == this->checkElements(i+1) ){
            open_gap_penalty = & _splice_sites_open_gap_penalty;
            extend_gap_penalty = & _splice_sites_extend_gap_penalty;
            vP = vProfile["splice_sites"] + ((nucleotideCodeSubstitutionMatrix.get_dna_acid_map()[_dna_d[i] ])*(this->_segLen)); // Correct part of the vProfile
        }

        // INDEL begin vector
        vGapO = _mm256_set1_epi32(*open_gap_penalty); // here I use epi32, because for the score matrix I used int32_t
        // INDEL extension vector
        vGapE = _mm256_set1_epi32(*extend_gap_penalty);

        vH = _mm256_loadu_si256(pvHStore+(this->_segLen-1));// the last one
        vH = _mm256_alignr_epi8(vH, _mm256_permute2x128_si256(vH, vH, _MM_SHUFFLE(0, 0, 2, 0)), 16 - 4);

        pv = pvHLoad;
        pvHLoad = pvHStore;
        pvHStore = pv;

        // inner loop to process the query sequence
        for (j = 0; j < this->_segLen; ++j) {
            vcmp = _mm256_loadu_si256(trackDel + j);
            __m256i t1;
            if( 0 == i ){
                vE = _mm256_add_epi32(vGapE, _mm256_loadu_si256(pvHLoad + j));
            } else {
                t1 = _mm256_add_epi32(_mm256_mullo_epi32 (_mm256_add_epi32(vcmp, vOne), vGapO),
                        _mm256_mullo_epi32(_mm256_mullo_epi32 (vcmp, vMinuteOne), vGapE));
                vE = _mm256_add_epi32(t1, _mm256_loadu_si256(pvHLoad + j));
            }

            vH = _mm256_add_epi32(vH, _mm256_loadu_si256(vP+j)); // add the scoring profile to vH
            vHsave = _mm256_add_epi32(vH, vZero);
            // Get max from vH, vE

            vH = _mm256_max_epi32(vH, vE);
            vcmp = _mm256_cmpeq_epi32(vH, vHsave);
            // Save vcmp values.
            _mm256_storeu_si256(trackMatch + j, vcmp);
            vcmp = _mm256_cmpeq_epi32(vH, vE);
            _mm256_storeu_si256(trackDel + j, vcmp);

            // Save vH values.
            _mm256_storeu_si256(pvHStore + j, vH);

            // Load the next vH.
            vH = _mm256_loadu_si256(pvHLoad + j);
        }

        // put the score into the score matrix begin
        for( j=0; j<this->_segLen; ++j ){
            a = ( int32_t*)(pvHStore+j);
            b = ( int32_t*)(trackMatch+j);
            c = ( int32_t*)(trackDel+j);
            for( z=0; z<8; ++z ){
                if( (j + z*this->_segLen) < _length_of_q ){
                    _similarity_matrix[j+z*this->_segLen+1][i+1] = a[z];
                    _track_match[j+z*this->_segLen+1][i+1] = b[z];
                    _track_del[j+z*this->_segLen+1][i+1] = c[z];
                }
            }
        }

        lastOne = 1; // insertion
        for( j=0; j<_length_of_q; ++j ){
            if( lastOne == 1 ){
                deletionScore = _similarity_matrix[j][i+1] + (int32_t)(*extend_gap_penalty);
                if( deletionScore > _similarity_matrix[j+1][i+1] ){
                    _similarity_matrix[j+1][i+1] = deletionScore;
                    _track_match[j+1][i+1] = 0;
                    _track_del[j+1][i+1] = 0;
                    _track_ins[j+1][i+1] = -1;

                    t=0;
                    while( j>=t && _track_ins[j-t][i+1]==-1 && _track_match[j-t][i+1]==-1 ){
                        _track_match[j-t][i+1] = 0;
                        ++t;
                    }
                    t=0;
                    while( j>=t && _track_ins[j-t][i+1]==-1 && _track_del[j-t][i+1]==-1 ){
                        _track_del[j-t][i+1] = 0;
                        ++t;
                    }

                } else if ( deletionScore == _similarity_matrix[j+1][i+1] ){
                    _track_ins[j+1][i+1] = -1;
                } else {
                    _track_ins[j+1][i+1] = 0;
                    lastOne=0;//MATCH OR MISMATCH
                }
            }else{
                deletionScore = _similarity_matrix[j][i+1]+(int32_t)(*open_gap_penalty);
                if( deletionScore > _similarity_matrix[j+1][i+1] ){
                    _similarity_matrix[j+1][i+1] = deletionScore;
                    _track_match[j+1][i+1] = 0;
                    _track_del[j+1][i+1] = 0;
                    _track_ins[j+1][i+1] = -1;
                    lastOne = 1;
                } else if ( deletionScore == _similarity_matrix[j+1][i+1] ){
                    _track_ins[j+1][i+1] = -1;
                    lastOne = 1;
                } else {
                    _track_ins[j+1][i+1] = 0;
                }
            }
            // put the score into the score matrix end
        }

        //update with track matrices _track_del begin
        for( j=0; j<_length_of_q; ++j ){
            if( _track_del[j+1][i+1]==-1 && _track_match[j+1][i+1]==0 && _track_ins[j+1][i+1]==0 ){
                t=0;
                while( i>=t && _track_del[j+1][i-t]==-1 && _track_match[j+1][i-t]==-1 ){
                    _track_match[j+1][i-t] = 0;
                    ++t;
                }
                t=0;
                while( i>=t && _track_del[j+1][i-t]==-1 && _track_ins[j+1][i-t]==-1 ){
                    _track_ins[j+1][i-t] = 0;
                    ++t;
                }
            }
        }
        //update with track matrices _track_del end


        //update with track matrices _track_match begin
        for( j=0; j<_length_of_q; ++j ){
            if( _track_match[j+1][i+1]==-1 && _track_del[j+1][i+1]==0 && _track_ins[j+1][i+1]==0 ){
                t=0;
                while( i>=t && j>=t && _track_del[j-t][i-t]==-1 && _track_match[j-t][i-t]==-1 ){
                    _track_del[j-t][i-t] = 0;
                    ++t;
                }
                t=0;
                while( i>=t && j>=t && _track_match[j-t][i-t]==-1 && _track_ins[j-t][i-t]==-1 ){
                    _track_ins[j-t][i-t] = 0;
                    ++t;
                }
            }
        }
        //update with track matrices _track_match end

        //updata vH begin
        for( j=0; j<this->_segLen; j++ ){
            a = (int32_t*)(pvHStore+j);
            b = ( int32_t*)(trackMatch+j);
            c = ( int32_t*)(trackDel+j);
            for( z=0; z<8; ++z ){
                if( (j + z*this->_segLen) < _length_of_q ){
                    a[z] = _similarity_matrix[j+z*this->_segLen+1][i+1];
                    b[z] = _track_match[j+z*this->_segLen+1][i+1];
                    c[z] = _track_del[j+z*this->_segLen+1][i+1];
                }
            }
        }//updata vH end
    }
    free(pvHStore);
    free(pvHLoad);
    free(trackMatch);
    free(trackDel);
    free(vProfile["exon"]);
    free(vProfile["intron"]);
    free(vProfile["splice_sites"]);
    free(vProfile["start_stop"]);
}

void alignNeedlemanForTranscript_simd_avx2int32::get_optimal_alignment( NucleotideCodeSubstitutionMatrix& _nucleotideCodeSubstitutionMatrix ) {

    _alignment_q = "";
    _alignment_d = "";
    _infor = "";
    std::stack<char> SQ, SD;
    std::stack<char> SI;

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
        SI.push('-');
    }
    while (i > 0 || j > 0) {
        char element;
        switch (checkElements(j)){
            case INTRON: element='I'; break;
            case EXON: element='E'; break;
            case START: element='S'; break;
            case STOP: element='T'; break;
            case SPLICEDONOR: element='P'; break;
            case SPLICEACCEPTOR: element='P'; break;
        }
        if (i == 0) {
            SQ.push('-');
            SD.push(_dna_d[j-1]);
            SI.push(element);
            --j;
        } else if (j == 0) {
            SQ.push(_dna_q[i-1]);
            SD.push('-');
            SI.push('-');
            --i;
        }else{
            if ( _track_match[i][j]==-1  ) {
                SQ.push(_dna_q[i - 1]);
                SD.push(_dna_d[j - 1]);
                SI.push(element);
                i -= 1;
                j -= 1;
            } else if (_track_ins[i][j] == -1 ) {//insertion
                SQ.push(_dna_q[i - 1]);
                SD.push('-');
                SI.push('-');
                i -= 1;
            } else {        // Going to S(i, j-1) //deletion
                SQ.push('-');
                SD.push(_dna_d[j - 1]);
                SI.push(element);
                j -= 1;
            }
        }
    }
    while (!SQ.empty()) {
        _alignment_q += SQ.top();
        _alignment_d += SD.top();
        _infor += SI.top();
        SQ.pop();
        SD.pop();
        SI.pop();
    }
}

void alignNeedlemanForTranscript_simd_avx2int32::print_results() {
    std::cout << _alignment_q << std::endl;
    std::cout << _alignment_d << std::endl;
    std::cout << _infor << std::endl;
}

std::string alignNeedlemanForTranscript_simd_avx2int32::getAlignment_q(){
    return _alignment_q;
}
std::string alignNeedlemanForTranscript_simd_avx2int32::getAlignment_d(){
    return _alignment_d;
}
