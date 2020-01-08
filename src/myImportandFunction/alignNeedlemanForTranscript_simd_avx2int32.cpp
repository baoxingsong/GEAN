/*
 * =====================================================================================
 *
 *       Filename:  alignNeedlemanForTranscript_simd.cpp
 *
 *    Description:
 *
 *        Version:  1.1
 *        Created:  30/December/2018 22:11:39
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Baoxing Song (songbx.me), songbaoxing168@163.com
 *   Organization:  The department of life science, Qiannan Normal University for Nationalities
 *
 * =====================================================================================
 */

/*************************************************************************


I borrowed some code from SSW Library Mengyao Zhao with License: MIT
https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library


************************************************************************/

#include "alignNeedlemanForTranscript_simd_avx2int32.h"





int checkElements( const int & startCodonPosition, const int & stopCodonPosition, const std::vector<SpliceSitePosition>& spliceSitePositions, const int & position ){
    if( position < startCodonPosition || position > stopCodonPosition+2  ){
        return 1;//INTRON;
    }
    for(std::vector<SpliceSitePosition>::size_type i = 0; i != spliceSitePositions.size(); i++ ){
        SpliceSitePosition spliceSite = spliceSitePositions[i];
        if(position > (spliceSite.getDonorSpliceSitePosition()+1) && position < (spliceSite.getAcceptorSpliceSitePosition()-1) ){
            return 1;//INTRON;
        }else if( position== spliceSite.getAcceptorSpliceSitePosition() || position == spliceSite.getAcceptorSpliceSitePosition()-1){
            return 2;//SPLICEACCEPTOR;
        }else if( position== spliceSite.getDonorSpliceSitePosition() || position == spliceSite.getDonorSpliceSitePosition()+1  ){
            return 2;//SPLICEDONOR;
        }
    }
    if (position == startCodonPosition || position == startCodonPosition+1 || position == startCodonPosition+2 ){
        return 3;//START;
    }
    if (position == stopCodonPosition || position == stopCodonPosition+1 || position == stopCodonPosition+2 ){
        return 3;//STOP;
    }
    return 0;//EXON;
}

char element_int_to_char( const int & element ){
    switch(element) {
        case  0: return 'E';
        case  1: return 'I';
        case  2: return 'S';
        case  3: return 'P';
    }
}

void alignNeedlemanForTranscript_simd_avx2int32_old(std::string& dna_d,
       std::string& dna_q, int & startCodonPosition, int & stopCodonPosition,
       std::vector<SpliceSitePosition>& spliceSitePositions,
       std::map<std::string, std::string>& parameters,
       NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix,
       std::string & alignment_q, std::string & alignment_d, std::string & infor){
    int32_t _exon_match_score = stoi(get_parameters("alignmentExonMatchP", parameters));
    int32_t _exon_mismatch_score = stoi(get_parameters("alignmentExonMismatchP", parameters));

    int32_t _intron_match_score = stoi(get_parameters("alignmentIntronMatchP", parameters));
    int32_t _intron_mismatch_score  = stoi(get_parameters("alignmentIntronMismatchP", parameters));

    int32_t _splice_sites_match_score = stoi(get_parameters("alignmentSpliceSitesMatchP", parameters));
    int32_t _splice_sites_mismatch_score = stoi(get_parameters("alignmentSpliceSitesMismatchP", parameters));

    int32_t _start_stop_codon_match_score = stoi(get_parameters("alignmentStartStopCodonMatchP", parameters));
    int32_t _start_stop_codon_mismatch_score = stoi(get_parameters("alignmentStartStopCodonMismatchP", parameters));

    int32_t * _open_gap_penalty = new int32_t[4];
    int32_t * _extend_gap_penalty = new int32_t[4];

    _open_gap_penalty[0]=stoi(get_parameters("alignmentExonOpenGapP", parameters));
    _extend_gap_penalty[0]=stoi(get_parameters("alignmentExonExtendGapP", parameters));

    _open_gap_penalty[1]=stoi(get_parameters("alignmentIntronOpenGapP", parameters));
    _extend_gap_penalty[1]=stoi(get_parameters("alignmentIntronExtendGapP", parameters));

    _open_gap_penalty[2]=stoi(get_parameters("alignmentSpliceSitesOpenGapP", parameters));
    _extend_gap_penalty[2]=stoi(get_parameters("alignmentSpliceSitesExtendGapP", parameters));

    _open_gap_penalty[3]=stoi(get_parameters("alignmentStartStopCodonOpenGapP", parameters));
    _extend_gap_penalty[3]=stoi(get_parameters("alignmentStartStopCodonExtendGapP", parameters));

    size_t length_of_q = dna_q.length();
    size_t length_of_d = dna_d.length();

    int8_t bias = 0;
    size_t n = 5;
    size_t segLen = ( length_of_q+ 7)/8;
    int32_t i, j, ai, k;
/*
    //check potential overflow problem begin
    int32_t maxLength = length_of_d;
    if( maxLength <  length_of_q){
        maxLength = length_of_q;
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
    if( (maxLength * maxScore) > pow(2,31)){
        std::cerr << "the sequence is too long, could not be optimized with this pipeline" << std::endl;
        return;
    }
    int32_t minimumPenalty = _exon_mismatch_score;

    if( minimumPenalty > _intron_mismatch_score ){
        minimumPenalty = _intron_mismatch_score;
    }
    if( minimumPenalty > _splice_sites_mismatch_score ){
        minimumPenalty = _splice_sites_mismatch_score;
    }
    if( minimumPenalty > _start_stop_codon_mismatch_score ){
        minimumPenalty = _start_stop_codon_mismatch_score;
    }

    for( i=0; i<4; ++i ){
        if( minimumPenalty > _open_gap_penalty[i] ){
            minimumPenalty =  _open_gap_penalty[i];
        }
    }
    if( (maxLength * minimumPenalty) < -pow(2,31)){
        std::cerr << "the sequence is too long, could not be optimized with this pipeline" << std::endl;
        return;
    } //check potential overflow problem end

    for( i=0; i<4; ++i ){
        if( _open_gap_penalty[i] > _extend_gap_penalty[i] ){
            std::cerr << "it is not reasonable to have open gap penalty larger than extend gap penalty" << std::endl;
        }
    }//warning message about score strategy end
*/
    int32_t * similarity_matrix = new int32_t[length_of_q + 1];
    for ( i = 0; i<= length_of_q; ++i) {
        similarity_matrix[i] = 0;
    }
    std::vector<bool> * trace_match = new std::vector<bool> [length_of_q + 1];
    std::vector<bool> * trace_del = new std::vector<bool> [length_of_q + 1];
    std::vector<bool> * trace_ins = new std::vector<bool> [length_of_q + 1];
    for ( i = 0; i < (length_of_q + 1); ++i) {
        trace_match[i].resize(length_of_d+1);
        trace_del[i].resize(length_of_d+1);
        trace_ins[i].resize(length_of_d+1);
    }

    for ( i = 0; i <= 1; i++) {
        for ( j = 0; j <= length_of_d; j++) {
            trace_match[i][j] = true;
            trace_del[i][j] = true;
            trace_ins[i][j] = true;
        }
    }
    for ( j = 0; j <= 1; j++){
        for ( i = 0; i <= length_of_q; i++) {
            trace_match[i][j] = true;
            trace_del[i][j] = true;
            trace_ins[i][j] = true;
        }
    }
    int8_t * ref_elements = new int8_t[length_of_d];
    for (i = 0; i<length_of_d; ++i) {//important
        ref_elements[i]=checkElements(startCodonPosition, stopCodonPosition, spliceSitePositions, i + 1); // because the gff records is start from 1
    }
    //std::cout << "alignNeedlemanForTranscript_simd_avx2int32 line 151" << std::endl;

   // std::cout << "alignNeedlemanForTranscript_simd_avx2int32 line 186" << std::endl;


    // vProfile begain
    __m256i* exon_profile = (__m256i*)malloc(n * segLen * sizeof(__m256i));
    __m256i* intron_profile = (__m256i*)malloc(n * segLen * sizeof(__m256i));
    __m256i* splice_sites_profile = (__m256i*)malloc(n * segLen * sizeof(__m256i));
    __m256i* start_stop_profile = (__m256i*)malloc(n * segLen * sizeof(__m256i));


    for( ai=0; ai<n; ++ai ){
        //vProfile[a] = new __m256i[segLen];
        int32_t* exon_s = (int32_t*)(exon_profile+ai*segLen);
        int32_t* intron_s = (int32_t*)(intron_profile+ai*segLen);
        int32_t* splice_sites_s = (int32_t*)(splice_sites_profile+ai*segLen);
        int32_t* start_stop_s = (int32_t*)(start_stop_profile+ai*segLen);
        for( i=0; i<segLen; ++i ){
            j = i;
            for( k=0; k<8; ++k ){
                *exon_s++ = j>= length_of_q ? bias : nucleotideCodeSubstitutionMatrix.get_exon_subsitition_matrix()[ai][nucleotideCodeSubstitutionMatrix.get_dna_acid_map(dna_q[j])];
                *intron_s++ = j>= length_of_q ? bias : nucleotideCodeSubstitutionMatrix.get_intron_subsitition_matrix()[ai][nucleotideCodeSubstitutionMatrix.get_dna_acid_map(dna_q[j])];
                *splice_sites_s++ = j>= length_of_q ? bias : nucleotideCodeSubstitutionMatrix.get_splice_sites_subsitition_matrix()[ai][nucleotideCodeSubstitutionMatrix.get_dna_acid_map(dna_q[j])];
                *start_stop_s++ = j>= length_of_q ? bias : nucleotideCodeSubstitutionMatrix.get_start_stop_codon_subsitition_matrix()[ai][nucleotideCodeSubstitutionMatrix.get_dna_acid_map(dna_q[j])];
                j += segLen;
            }
        }
    }
    __m256i** vProfile =  new __m256i*[4];
    vProfile[0] = exon_profile;
    vProfile[1] = intron_profile;
    vProfile[2] = splice_sites_profile;
    vProfile[3] = start_stop_profile;
    // vProfile end

    //std::cout << "alignNeedlemanForTranscript_simd_avx2int32 line 220  vProfile end" << std::endl;


    //calculate_similarity begin

    const __m256i vZero = _mm256_set1_epi32(0);
    const __m256i minimum = _mm256_set1_epi32(-1147483648);
    __m256i* pvHStore = (__m256i*) calloc(segLen, sizeof(__m256i));
    __m256i* pvHLoad = (__m256i*) calloc(segLen, sizeof(__m256i));
    __m256i* pvE = (__m256i*) calloc(segLen, sizeof(__m256i));
    __m256i* pvEC = (__m256i*) calloc(segLen, sizeof(__m256i));
    __m256i* pvMC = (__m256i*) calloc(segLen, sizeof(__m256i));
    __m256i* pv;
    __m256i* vP;
    __m256i e;
    __m256i vGapO, vGapE, next_vGapO, next_vGapE, vTemp, vTemp2;
    int z;
    int32_t cmp;
    int32_t * a;
    int32_t * b;
    int32_t * c;
    int32_t * d;
    //for outer loop
    __m256i vH;
    for (i = 0; i<length_of_d; ++i) {//important
        //std::cout << "alignNeedlemanForTranscript_simd_avx2int32 line 251" << std::endl;
        //std::cout << "alignNeedlemanForTranscript_simd_avx2int32 line 255" << std::endl;
        vP = vProfile[ref_elements[i]] + ((nucleotideCodeSubstitutionMatrix.get_dna_acid_map(dna_d[i]))*(segLen)); // Correct part of the vProfile
        //std::cout << "alignNeedlemanForTranscript_simd_avx2int32 line 257" << std::endl;
        // INDEL begin vector
        vGapO = _mm256_set1_epi32(_open_gap_penalty[ref_elements[i]]); // here I use epi32, because for the score matrix I used int32_t
        // INDEL extension vector
        vGapE = _mm256_set1_epi32(_extend_gap_penalty[ref_elements[i]]);
        if( i == length_of_d-1){
            next_vGapO=vGapO;
            next_vGapE=vGapE;
        }else{
            next_vGapO=_mm256_set1_epi32(_open_gap_penalty[ref_elements[i+1]]);
            next_vGapE= _mm256_set1_epi32(_extend_gap_penalty[ref_elements[i+1]]);
        }
        //std::cout << "alignNeedlemanForTranscript_simd_avx2int32 line 269" << std::endl;
        vH = _mm256_loadu_si256(pvHStore+(segLen-1));// the last one
        //vH = pvHStore[_segLen - 1];
        //std::cout << "alignNeedlemanForTranscript_simd_avx2int32 line 271" << std::endl;
        vH = _mm256_alignr_epi8(vH, _mm256_permute2x128_si256(vH, vH, _MM_SHUFFLE(0, 0, 2, 0)), 16 - 4);

        //vMaxColumn = vZero; /* vMaxColumn is used to record the max values of column i. */
        pv = pvHLoad;
        pvHLoad = pvHStore;
        pvHStore = pv;
        //std::cout << "alignNeedlemanForTranscript_simd_avx2int32 line 280" << std::endl;
        // inner loop to process the query sequence
//        std::cout << "alignNeedlemanForTranscript_simd_avx2int32 line 272" << std::endl;
        for (j = 0; j < segLen; ++j) {
            if( 0 == i ){
                e = vGapE;
            }else{
                e = _mm256_loadu_si256(pvE + j);
            }
            vH = _mm256_add_epi32(vH, _mm256_loadu_si256(vP+j)); // add the scoring profile to vH
            _mm256_storeu_si256(pvEC+j, e);
            _mm256_storeu_si256(pvMC+j, vH);
            vH = _mm256_max_epi32(vH, e);
            //vMaxColumn = _mm256_max_epi32(vMaxColumn, vH);

            /* Save vH values. */
            _mm256_storeu_si256(pvHStore + j, vH);
            /* Update vE value. */
            vH = _mm256_add_epi32(vH, next_vGapO); // this is for next column
            e = _mm256_add_epi32(e, next_vGapE);
            e = _mm256_max_epi32(e, vH);
            _mm256_storeu_si256(pvE + j, e);

            /* Load the next vH. */
            vH = _mm256_loadu_si256(pvHLoad + j);
        }

        for( j=0; j<segLen; ++j ) {
            a = (int32_t *) (pvHStore + j);
            b = (int32_t *) (pvMC + j);
            c = (int32_t *) (pvEC + j);
            for (z = 0; z < 8; ++z) {
                if ((j + z * segLen) < length_of_q) {
                    similarity_matrix[j+z*segLen+1] = a[z];
                    trace_match[j + z * segLen + 1][i + 1] = (a[z] == b[z]);
                    trace_del[j + z * segLen + 1][i + 1] = (a[z] == c[z]);
                }
            }
        }

        int lastOne = 1, t; // insertion
        int32_t deletionScore;
        for( j=0; j<length_of_q; ++j ){
            if( lastOne == 1 ){
                deletionScore = similarity_matrix[j] + _extend_gap_penalty[ref_elements[i]];
                if( deletionScore > similarity_matrix[j+1] ){
                    similarity_matrix[j+1] = deletionScore;
                    trace_match[j+1][i+1] = false;
                    trace_del[j+1][i+1] = false;
                    trace_ins[j+1][i+1] = true;

                    t=0;
                    while( j>=t && trace_ins[j-t][i+1] && trace_match[j-t][i+1] ){
                        trace_match[j-t][i+1] = false;
                        ++t;
                    }
                    t=0;
                    while( j>=t && trace_ins[j-t][i+1] && trace_del[j-t][i+1] ){
                        trace_del[j-t][i+1] = false;
                        ++t;
                    }
                } else if ( deletionScore == similarity_matrix[j+1] ){
                    trace_ins[j+1][i+1] = true;
                } else {
                    trace_ins[j+1][i+1] = false;
                    lastOne=0;//MATCH OR MISMATCH
                }
            }else{
                deletionScore = similarity_matrix[j] + _open_gap_penalty[ref_elements[i]];
                if( deletionScore > similarity_matrix[j+1] ){
                    similarity_matrix[j+1] = deletionScore;
                    trace_match[j+1][i+1] = false;
                    trace_del[j+1][i+1] = false;
                    trace_ins[j+1][i+1] = true;
                    lastOne = 1;
                } else if ( deletionScore == similarity_matrix[j+1] ){
                    trace_ins[j+1][i+1] = true;
                    lastOne = 1;
                } else {
                    trace_ins[j+1][i+1] = false;
                }
            }
            // put the score into the score matrix end
        }

        //update with track matrices _track_del begin
        for( j=0; j<length_of_q; ++j ){
            if( trace_del[j+1][i+1] && trace_match[j+1][i+1]==false && !trace_ins[j+1][i+1] ){
                t=0;
                while( i>=t && trace_del[j+1][i-t] && trace_match[j+1][i-t] ){
                    trace_match[j+1][i-t] = false;
                    ++t;
                }
                t=0;
                while( i>=t && trace_del[j+1][i-t] && trace_ins[j+1][i-t] ){
                    trace_ins[j+1][i-t] = false;
                    ++t;
                }
            }
        }
        //update with track matrices _track_del end

        //update with track matrices _track_match begin
        for( j=0; j<length_of_q; ++j ){
            if( trace_match[j+1][i+1] && !trace_del[j+1][i+1] && !trace_ins[j+1][i+1] ){
                t=0;
                while( i>=t && j>=t && trace_del[j-t][i-t] && trace_match[j-t][i-t] ){
                    trace_del[j-t][i-t] = false;
                    ++t;
                }
                t=0;
                while( i>=t && j>=t && trace_match[j-t][i-t] && trace_ins[j-t][i-t] ){
                    trace_ins[j-t][i-t] = false;
                    ++t;
                }
            }
        }
        //update with track matrices _track_match end

        //updata vH begin
        for( j=0; j<segLen; j++ ){
            a = (int32_t*)(pvHStore+j);
            for( z=0; z<8; ++z ){
                if( (j + z*segLen) < length_of_q ){
                    a[z] = similarity_matrix[j+z*segLen+1];
                }
            }
        }//updata vH end
    }
//    std::cout << "alignNeedlemanForTranscript_simd_avx2int32 line 399" << std::endl;
    for( j=0; j<segLen; ++j ){
        a = ( int32_t*)(pvHStore+j);
        for( z=0; z<8; ++z ){
            if( (j + z*segLen) < length_of_q ){
                similarity_matrix[j+z*segLen+1] = a[z];
            }
        }
    }
    //std::cout << "alignNeedlemanForTranscript_simd_avx2int32 line 335" << std::endl;
    free(pvHStore);
    //std::cout << "alignNeedlemanForTranscript_simd_avx2int32 line 421" << std::endl;
    free(pvHLoad);
    //std::cout << "alignNeedlemanForTranscript_simd_avx2int32 line 423" << std::endl;
    free(pvE);
    //std::cout << "alignNeedlemanForTranscript_simd_avx2int32 line 425" << std::endl;
    free(pvEC);
    free(pvMC);
    //std::cout << "alignNeedlemanForTranscript_simd_avx2int32 line 426" << std::endl;
    free(vProfile[0]);
    free(vProfile[1]);
    free(vProfile[2]);
    free(vProfile[3]);
    //std::cout << "alignNeedlemanForTranscript_simd_avx2int32 line 429" << std::endl;
    delete [] vProfile;
    //calculate_similarity end


    //std::cout << "alignNeedlemanForTranscript_simd_avx2int32 line 426  calculate_similarity end" << std::endl;

    // get_optimal_alignment begin

    std::stack<char> SQ, SD, SI;

    i = length_of_q;
    j = length_of_d;
    int32_t highestScore = similarity_matrix[i];
    for( k=length_of_q; k >0; --k ){
        if( similarity_matrix[k] > highestScore ){
            i = k;
            highestScore= similarity_matrix[k];
        }
    }
    for( k=length_of_q; k>i; --k ){
        SQ.push(dna_q[k-1]);
        SD.push('-');
        SI.push('-');
    }
    while (i > 0 || j > 0) {
        char element = element_int_to_char(ref_elements[j - 1]);
        if (i == 0) {
            SQ.push('-');
            SD.push(dna_d[j - 1]);
            SI.push(element);
            --j;
        } else if (j == 0) {
            SQ.push(dna_q[i - 1]);
            SD.push('-');
            SI.push('-');
            --i;
        } else {
            if (trace_match[i][j]) {
                SQ.push(dna_q[i - 1]);
                SD.push(dna_d[j - 1]);
                SI.push(element);
                --i;
                --j;
            } /*else if (_trace_del[i][j]) {//insertion
                SQ.push('-');
                SD.push(_dna_d[j - 1]);
                SI.push(element);
                --j;
            } else {        // Going to S(i, j-1) //deletion
                SQ.push(_dna_q[i - 1]);
                SD.push('-');
                SI.push('-');
                --i;
            }*/

            else if (trace_ins[i][j] ) {//insertion
                SQ.push(dna_q[i - 1]);
                SD.push('-');
                SI.push('-');
                --i;
            } else {        // Going to S(i, j-1) //deletion
                SQ.push('-');
                SD.push(dna_d[j - 1]);
                SI.push(element);
                --j;
            }

        }
    }
    while (!SQ.empty()) {
        alignment_q += SQ.top();
        alignment_d += SD.top();
        infor += SI.top();
        SQ.pop();
        SD.pop();
        SI.pop();
    }
    // get_optimal_alignment end



    delete[] _open_gap_penalty;
    delete[] _extend_gap_penalty;
    delete[] similarity_matrix;
    delete[] trace_match;
    delete[] trace_del;
    delete[] trace_ins;
    delete[] ref_elements;
}























































































void alignNeedlemanForTranscript_simd_avx2int32(std::string& dna_d,
                                                    std::string& dna_q, int & startCodonPosition, int & stopCodonPosition,
                                                    std::vector<SpliceSitePosition>& spliceSitePositions,
                                                    std::map<std::string, std::string>& parameters,
                                                    NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix,
                                                    std::string & alignment_q, std::string & alignment_d, std::string & infor){
    int32_t _exon_match_score = stoi(get_parameters("alignmentExonMatchP", parameters));
    int32_t _exon_mismatch_score = stoi(get_parameters("alignmentExonMismatchP", parameters));

    int32_t _intron_match_score = stoi(get_parameters("alignmentIntronMatchP", parameters));
    int32_t _intron_mismatch_score  = stoi(get_parameters("alignmentIntronMismatchP", parameters));

    int32_t _splice_sites_match_score = stoi(get_parameters("alignmentSpliceSitesMatchP", parameters));
    int32_t _splice_sites_mismatch_score = stoi(get_parameters("alignmentSpliceSitesMismatchP", parameters));

    int32_t _start_stop_codon_match_score = stoi(get_parameters("alignmentStartStopCodonMatchP", parameters));
    int32_t _start_stop_codon_mismatch_score = stoi(get_parameters("alignmentStartStopCodonMismatchP", parameters));

    int32_t * _open_gap_penalty = new int32_t[4];
    int32_t * _extend_gap_penalty = new int32_t[4];

    _open_gap_penalty[0]=stoi(get_parameters("alignmentExonOpenGapP", parameters));
    _extend_gap_penalty[0]=stoi(get_parameters("alignmentExonExtendGapP", parameters));

    _open_gap_penalty[1]=stoi(get_parameters("alignmentIntronOpenGapP", parameters));
    _extend_gap_penalty[1]=stoi(get_parameters("alignmentIntronExtendGapP", parameters));

    _open_gap_penalty[2]=stoi(get_parameters("alignmentSpliceSitesOpenGapP", parameters));
    _extend_gap_penalty[2]=stoi(get_parameters("alignmentSpliceSitesExtendGapP", parameters));

    _open_gap_penalty[3]=stoi(get_parameters("alignmentStartStopCodonOpenGapP", parameters));
    _extend_gap_penalty[3]=stoi(get_parameters("alignmentStartStopCodonExtendGapP", parameters));

    size_t length_of_q = dna_q.length();
    size_t length_of_d = dna_d.length();

    int8_t bias = 0;
    size_t n = 5;
    size_t segLen = ( length_of_q+ 7)/8;
    int32_t i, j, ai, k;

    int32_t * similarity_matrix = new int32_t[length_of_q + 1];
    for ( i = 0; i<= length_of_q; ++i) {
        similarity_matrix[i] = 0;
    }
    std::vector<bool> * trace_match = new std::vector<bool> [length_of_q + 1];
    std::vector<bool> * trace_del = new std::vector<bool> [length_of_q + 1];
    std::vector<bool> * trace_ins = new std::vector<bool> [length_of_q + 1];
    for ( i = 0; i < (length_of_q + 1); ++i) {
        trace_match[i].resize(length_of_d+1);
        trace_del[i].resize(length_of_d+1);
        trace_ins[i].resize(length_of_d+1);
    }

    for ( i = 0; i <= 1; i++) {
        for ( j = 0; j <= length_of_d; j++) {
            trace_match[i][j] = true;
            trace_del[i][j] = true;
            trace_ins[i][j] = true;
        }
    }
    for ( j = 0; j <= 1; j++){
        for ( i = 0; i <= length_of_q; i++) {
            trace_match[i][j] = true;
            trace_del[i][j] = true;
            trace_ins[i][j] = true;
        }
    }
    int8_t * ref_elements = new int8_t[length_of_d];
    for (i = 0; i<length_of_d; ++i) {//important
        ref_elements[i]=checkElements(startCodonPosition, stopCodonPosition, spliceSitePositions, i + 1); // because the gff records is start from 1
    }

    // vProfile begain
    __m256i* exon_profile = (__m256i*)malloc(n * segLen * sizeof(__m256i));
    __m256i* intron_profile = (__m256i*)malloc(n * segLen * sizeof(__m256i));
    __m256i* splice_sites_profile = (__m256i*)malloc(n * segLen * sizeof(__m256i));
    __m256i* start_stop_profile = (__m256i*)malloc(n * segLen * sizeof(__m256i));

    for( ai=0; ai<n; ++ai ){
        //vProfile[a] = new __m256i[segLen];
        int32_t* exon_s = (int32_t*)(exon_profile+ai*segLen);
        int32_t* intron_s = (int32_t*)(intron_profile+ai*segLen);
        int32_t* splice_sites_s = (int32_t*)(splice_sites_profile+ai*segLen);
        int32_t* start_stop_s = (int32_t*)(start_stop_profile+ai*segLen);
        for( i=0; i<segLen; ++i ){
            j = i;
            for( k=0; k<8; ++k ){
                *exon_s++ = j>= length_of_q ? bias : nucleotideCodeSubstitutionMatrix.get_exon_subsitition_matrix()[ai][nucleotideCodeSubstitutionMatrix.get_dna_acid_map(dna_q[j])];
                *intron_s++ = j>= length_of_q ? bias : nucleotideCodeSubstitutionMatrix.get_intron_subsitition_matrix()[ai][nucleotideCodeSubstitutionMatrix.get_dna_acid_map(dna_q[j])];
                *splice_sites_s++ = j>= length_of_q ? bias : nucleotideCodeSubstitutionMatrix.get_splice_sites_subsitition_matrix()[ai][nucleotideCodeSubstitutionMatrix.get_dna_acid_map(dna_q[j])];
                *start_stop_s++ = j>= length_of_q ? bias : nucleotideCodeSubstitutionMatrix.get_start_stop_codon_subsitition_matrix()[ai][nucleotideCodeSubstitutionMatrix.get_dna_acid_map(dna_q[j])];
                j += segLen;
            }
        }
    }
    __m256i** vProfile =  new __m256i*[4];
    vProfile[0] = exon_profile;
    vProfile[1] = intron_profile;
    vProfile[2] = splice_sites_profile;
    vProfile[3] = start_stop_profile;
    // vProfile end

    //calculate_similarity begin
    const __m256i vZero = _mm256_set1_epi32(0);
    const __m256i minimum = _mm256_set1_epi32(-1147483648);
    __m256i* pvHStore = (__m256i*) calloc(segLen, sizeof(__m256i));
    __m256i* pvHLoad = (__m256i*) calloc(segLen, sizeof(__m256i));
    __m256i* pvE = (__m256i*) calloc(segLen, sizeof(__m256i));
    __m256i* pvEC = (__m256i*) calloc(segLen, sizeof(__m256i));
    __m256i* pvMC = (__m256i*) calloc(segLen, sizeof(__m256i));
    __m256i* pv;
    __m256i* vP;
    __m256i e;
    __m256i vGapO, vGapE, next_vGapO, next_vGapE, vTemp, vTemp2;
    int z;
    int32_t cmp;
    int32_t * a;
    int32_t * b;
    int32_t * c;
    int32_t * d;
    //for outer loop
    __m256i vH;
    for (i = 0; i<length_of_d; ++i) {//important
        vP = vProfile[ref_elements[i]] + ((nucleotideCodeSubstitutionMatrix.get_dna_acid_map(dna_d[i]))*(segLen)); // Correct part of the vProfile
        // INDEL begin vector
        vGapO = _mm256_set1_epi32(_open_gap_penalty[ref_elements[i]]); // here I use epi32, because for the score matrix I used int32_t
        // INDEL extension vector
        vGapE = _mm256_set1_epi32(_extend_gap_penalty[ref_elements[i]]);
        if( i == length_of_d-1){
            next_vGapO=vGapO;
            next_vGapE=vGapE;
        }else{
            next_vGapO=_mm256_set1_epi32(_open_gap_penalty[ref_elements[i+1]]);
            next_vGapE= _mm256_set1_epi32(_extend_gap_penalty[ref_elements[i+1]]);
        }
        vH = _mm256_loadu_si256(pvHStore+(segLen-1));// the last one

        vH = _mm256_alignr_epi8(vH, _mm256_permute2x128_si256(vH, vH, _MM_SHUFFLE(0, 0, 2, 0)), 16 - 4);

        pv = pvHLoad;
        pvHLoad = pvHStore;
        pvHStore = pv;

        for (j = 0; j < segLen; ++j) {
            if( 0 == i ){
                e = vGapE;
            }else{
                e = _mm256_loadu_si256(pvE + j);
            }
            vH = _mm256_add_epi32(vH, _mm256_loadu_si256(vP+j)); // add the scoring profile to vH
            _mm256_storeu_si256(pvEC+j, e);
            _mm256_storeu_si256(pvMC+j, vH);
            vH = _mm256_max_epi32(vH, e);
            //vMaxColumn = _mm256_max_epi32(vMaxColumn, vH);

            /* Save vH values. */
            _mm256_storeu_si256(pvHStore + j, vH);
            /* Update vE value. */
            vH = _mm256_add_epi32(vH, next_vGapO); // this is for next column
            e = _mm256_add_epi32(e, next_vGapE);
            e = _mm256_max_epi32(e, vH);
            _mm256_storeu_si256(pvE + j, e);

            /* Load the next vH. */
            vH = _mm256_loadu_si256(pvHLoad + j);
        }

        for( j=0; j<segLen; ++j ) {
            a = (int32_t *) (pvHStore + j);
            b = (int32_t *) (pvMC + j);
            c = (int32_t *) (pvEC + j);
            for (z = 0; z < 8; ++z) {
                if ((j + z * segLen) < length_of_q) {
                    similarity_matrix[j+z*segLen+1] = a[z];
                    trace_match[j + z * segLen + 1][i + 1] = (a[z] == b[z]);
                    trace_del[j + z * segLen + 1][i + 1] = (a[z] == c[z]);
                }
            }
        }

        int lastOne = 1, t; // insertion
        int32_t deletionScore;
        for( j=0; j<length_of_q; ++j ){
            if( lastOne == 1 ){
                deletionScore = similarity_matrix[j] + _extend_gap_penalty[ref_elements[i]];
                if( deletionScore > similarity_matrix[j+1] ){
                    similarity_matrix[j+1] = deletionScore;
                    trace_match[j+1][i+1] = false;
                    trace_del[j+1][i+1] = false;
                    trace_ins[j+1][i+1] = true;

                    t=0;
                    while( j>=t && trace_ins[j-t][i+1] && trace_match[j-t][i+1] ){
                        trace_match[j-t][i+1] = false;
                        ++t;
                    }
                    t=0;
                    while( j>=t && trace_ins[j-t][i+1] && trace_del[j-t][i+1] ){
                        trace_del[j-t][i+1] = false;
                        ++t;
                    }
                } else if ( deletionScore == similarity_matrix[j+1] ){
                    trace_ins[j+1][i+1] = true;
                } else {
                    trace_ins[j+1][i+1] = false;
                    lastOne=0;//MATCH OR MISMATCH
                }
            }else{
                deletionScore = similarity_matrix[j] + _open_gap_penalty[ref_elements[i]];
                if( deletionScore > similarity_matrix[j+1] ){
                    similarity_matrix[j+1] = deletionScore;
                    trace_match[j+1][i+1] = false;
                    trace_del[j+1][i+1] = false;
                    trace_ins[j+1][i+1] = true;
                    lastOne = 1;
                } else if ( deletionScore == similarity_matrix[j+1] ){
                    trace_ins[j+1][i+1] = true;
                    lastOne = 1;
                } else {
                    trace_ins[j+1][i+1] = false;
                }
            }// put the score into the score matrix end
        }

        //update with track matrices _track_del begin
        for( j=0; j<length_of_q; ++j ){
            if( trace_del[j+1][i+1] && trace_match[j+1][i+1]==false && !trace_ins[j+1][i+1] ){
                t=0;
                while( i>=t && trace_del[j+1][i-t] && trace_match[j+1][i-t] ){
                    trace_match[j+1][i-t] = false;
                    ++t;
                }
                t=0;
                while( i>=t && trace_del[j+1][i-t] && trace_ins[j+1][i-t] ){
                    trace_ins[j+1][i-t] = false;
                    ++t;
                }
            }
        }
        //update with track matrices _track_del end

        //update with track matrices _track_match begin
        for( j=0; j<length_of_q; ++j ){
            if( trace_match[j+1][i+1] && !trace_del[j+1][i+1] && !trace_ins[j+1][i+1] ){
                t=0;
                while( i>=t && j>=t && trace_del[j-t][i-t] && trace_match[j-t][i-t] ){
                    trace_del[j-t][i-t] = false;
                    ++t;
                }
                t=0;
                while( i>=t && j>=t && trace_match[j-t][i-t] && trace_ins[j-t][i-t] ){
                    trace_ins[j-t][i-t] = false;
                    ++t;
                }
            }
        }//update with track matrices _track_match end

        //updata vH begin
        for( j=0; j<segLen; j++ ){
            a = (int32_t*)(pvHStore+j);
            for( z=0; z<8; ++z ){
                if( (j + z*segLen) < length_of_q ){
                    a[z] = similarity_matrix[j+z*segLen+1];
                }
            }
        }//updata vH end
    }

    for( j=0; j<segLen; ++j ){
        a = ( int32_t*)(pvHStore+j);
        for( z=0; z<8; ++z ){
            if( (j + z*segLen) < length_of_q ){
                similarity_matrix[j+z*segLen+1] = a[z];
            }
        }
    }
    free(pvHStore);
    free(pvHLoad);
    free(pvE);
    free(pvEC);
    free(pvMC);
    free(vProfile[0]);
    free(vProfile[1]);
    free(vProfile[2]);
    free(vProfile[3]);
    delete [] vProfile;
    //calculate_similarity end

    // get_optimal_alignment begin
    std::stack<char> SQ, SD, SI;
    i = length_of_q;
    j = length_of_d;
    int32_t highestScore = similarity_matrix[i];
    for( k=length_of_q; k >0; --k ){
        if( similarity_matrix[k] > highestScore ){
            i = k;
            highestScore= similarity_matrix[k];
        }
    }
    for( k=length_of_q; k>i; --k ){
        SQ.push(dna_q[k-1]);
        SD.push('-');
        SI.push('-');
    }
    while (i > 0 || j > 0) {
        char element = element_int_to_char(ref_elements[j - 1]);
        if (i == 0) {
            SQ.push('-');
            SD.push(dna_d[j - 1]);
            SI.push(element);
            --j;
        } else if (j == 0) {
            SQ.push(dna_q[i - 1]);
            SD.push('-');
            SI.push('-');
            --i;
        } else {
            if (trace_match[i][j]) {
                SQ.push(dna_q[i - 1]);
                SD.push(dna_d[j - 1]);
                SI.push(element);
                --i;
                --j;
            }else if (trace_ins[i][j] ) {//insertion
                SQ.push(dna_q[i - 1]);
                SD.push('-');
                SI.push('-');
                --i;
            } else {        // Going to S(i, j-1) //deletion
                SQ.push('-');
                SD.push(dna_d[j - 1]);
                SI.push(element);
                --j;
            }
        }
    }
    while (!SQ.empty()) {
        alignment_q += SQ.top();
        alignment_d += SD.top();
        infor += SI.top();
        SQ.pop();
        SD.pop();
        SI.pop();
    }
    // get_optimal_alignment end


    delete[] _open_gap_penalty;
    delete[] _extend_gap_penalty;
    delete[] similarity_matrix;
    delete[] trace_match;
    delete[] trace_del;
    delete[] trace_ins;
    delete[] ref_elements;
}