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

alignNeedlemanForTranscript_simd_avx2int32::alignNeedlemanForTranscript_simd_avx2int32(std::string& dna_d,
                                                                                       std::string& dna_q, int & startCodonPosition, int & stopCodonPosition,
                                                                                       std::vector<SpliceSitePosition>& splitSitePositions,
                                                                                       std::map<std::string, std::string>& parameters, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix){
    //std::cout << "NeedlemanWunschForTranscript begin" << std::endl;
    this->_dna_q = dna_q;
    this->_dna_d = dna_d;
    this->_exon_match_score = stoi(get_parameters("alignmentExonMatchP", parameters));
    this->_exon_mismatch_score = stoi(get_parameters("alignmentExonMismatchP", parameters));

    this->_intron_match_score = stoi(get_parameters("alignmentIntronMatchP", parameters));
    this->_intron_mismatch_score  = stoi(get_parameters("alignmentIntronMismatchP", parameters));

    this->_splice_sites_match_score = stoi(get_parameters("alignmentSpliceSitesMatchP", parameters));
    this->_splice_sites_mismatch_score = stoi(get_parameters("alignmentSpliceSitesMismatchP", parameters));

    this->_start_stop_codon_match_score = stoi(get_parameters("alignmentStartStopCodonMatchP", parameters));
    this->_start_stop_codon_mismatch_score = stoi(get_parameters("alignmentStartStopCodonMismatchP", parameters));

    _open_gap_penalty = new int32_t[4];
    _extend_gap_penalty = new int32_t[4];
    _open_gap_penalty[0]=stoi(get_parameters("alignmentExonOpenGapP", parameters));
    _extend_gap_penalty[0]=stoi(get_parameters("alignmentExonExtendGapP", parameters));

    _open_gap_penalty[1]=stoi(get_parameters("alignmentIntronOpenGapP", parameters));
    _extend_gap_penalty[1]=stoi(get_parameters("alignmentIntronExtendGapP", parameters));

    _open_gap_penalty[2]=stoi(get_parameters("alignmentSpliceSitesOpenGapP", parameters));
    _extend_gap_penalty[2]=stoi(get_parameters("alignmentSpliceSitesExtendGapP", parameters));

    _open_gap_penalty[3]=stoi(get_parameters("alignmentStartStopCodonOpenGapP", parameters));
    _extend_gap_penalty[3]=stoi(get_parameters("alignmentStartStopCodonExtendGapP", parameters));

    this->_startCodonPosition = startCodonPosition;
    this->_stopCodonPosition = stopCodonPosition;
    this->_spliceSitePositions = splitSitePositions;

    this->_length_of_q = this->_dna_q.length();
    this->_length_of_d = this->_dna_d.length();

    this->bias = 0;
    this->_n = 5;
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
    int32_t i, j;
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

    this->_similarity_matrix = new int32_t[this->_length_of_q + 1];
    for ( i = 0; i<= _length_of_q; ++i) {
        _similarity_matrix[i] = 0;
    }
    this->_trace_match = new std::vector<bool> [this->_length_of_q + 1];
    this->_trace_del = new std::vector<bool> [this->_length_of_q + 1];
    this->_trace_ins = new std::vector<bool> [this->_length_of_q + 1];
    for ( i = 0; i < (this->_length_of_q + 1); ++i) {
        _trace_match[i].resize(_length_of_d+1);
        _trace_del[i].resize(_length_of_d+1);
        _trace_ins[i].resize(_length_of_d+1);
    }

    for ( i = 0; i <= 1; i++) {
        for ( j = 0; j <= _length_of_d; j++) {
            _trace_match[i][j] = true;
            _trace_del[i][j] = true;
            _trace_ins[i][j] = true;
        }
    }
    for ( j = 0; j <= 1; j++){
        for ( i = 0; i <= _length_of_q; i++) {
            _trace_match[i][j] = true;
            _trace_del[i][j] = true;
            _trace_ins[i][j] = true;
        }
    }
    _ref_elements = new int8_t[this->_length_of_d];
    for (i = 0; i<_length_of_d; ++i) {//important
        _ref_elements[i]=this->checkElements(i + 1); // because the gff records is start from 1
    }
    //std::cout << "alignNeedlemanForTranscript_simd_avx2int32 line 151" << std::endl;
    __m256i** vProfile = query_profile_avx2_byte(nucleotideCodeSubstitutionMatrix);
    //std::cout << "alignNeedlemanForTranscript_simd_avx2int32 line 153" << std::endl;
    this->calculate_similarity(vProfile, nucleotideCodeSubstitutionMatrix);
    //std::cout << "alignNeedlemanForTranscript_simd_avx2int32 line 155" << std::endl;
    this->get_optimal_alignment(nucleotideCodeSubstitutionMatrix);
    //std::cout << "alignNeedlemanForTranscript_simd_avx2int32 line 157" << std::endl;
   // this->print_results();
}

alignNeedlemanForTranscript_simd_avx2int32::~alignNeedlemanForTranscript_simd_avx2int32() {
    delete[] this->_open_gap_penalty;
    delete[] this->_extend_gap_penalty;
    delete[] this->_similarity_matrix;
    delete[] this->_trace_match;
    delete[] this->_trace_del;
    delete[] this->_trace_ins;
    delete[] this->_ref_elements;
}

int alignNeedlemanForTranscript_simd_avx2int32::checkElements( int position ){
    if( position < this->_startCodonPosition || position > this->_stopCodonPosition+2  ){
        return 1;//INTRON;
    }
    for(std::vector<SpliceSitePosition>::size_type i = 0; i != this->_spliceSitePositions.size(); i++ ){
        SpliceSitePosition spliceSite = this->_spliceSitePositions[i];
        if(position > (spliceSite.getDonorSpliceSitePosition()+1) && position < (spliceSite.getAcceptorSpliceSitePosition()-1) ){
            return 1;//INTRON;
        }else if( position== spliceSite.getAcceptorSpliceSitePosition() || position == spliceSite.getAcceptorSpliceSitePosition()-1){
            return 2;//SPLICEACCEPTOR;
        }else if( position== spliceSite.getDonorSpliceSitePosition() || position == spliceSite.getDonorSpliceSitePosition()+1  ){
            return 2;//SPLICEDONOR;
        }
    }
    if (position == this->_startCodonPosition || position == this->_startCodonPosition+1 || position == this->_startCodonPosition+2 ){
        return 3;//START;
    }
    if (position == this->_stopCodonPosition || position == this->_stopCodonPosition+1 || position == this->_stopCodonPosition+2 ){
        return 3;//STOP;
    }
    return 0;//EXON;
}

__m256i** alignNeedlemanForTranscript_simd_avx2int32::query_profile_avx2_byte(
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
                *exon_s++ = j>= _length_of_q ? bias : nucleotideCodeSubstitutionMatrix.get_exon_subsitition_matrix()[a][nucleotideCodeSubstitutionMatrix.get_dna_acid_map(_dna_q[j])];
                *intron_s++ = j>= _length_of_q ? bias : nucleotideCodeSubstitutionMatrix.get_intron_subsitition_matrix()[a][nucleotideCodeSubstitutionMatrix.get_dna_acid_map(_dna_q[j])];
                *splice_sites_s++ = j>= _length_of_q ? bias : nucleotideCodeSubstitutionMatrix.get_splice_sites_subsitition_matrix()[a][nucleotideCodeSubstitutionMatrix.get_dna_acid_map(_dna_q[j])];
                *start_stop_s++ = j>= _length_of_q ? bias : nucleotideCodeSubstitutionMatrix.get_start_stop_codon_subsitition_matrix()[a][nucleotideCodeSubstitutionMatrix.get_dna_acid_map(_dna_q[j])];
                j += this->_segLen;
            }
        }
    }
    __m256i** vProfile =  new __m256i*[4];
    vProfile[0] = exon_profile;
    vProfile[1] = intron_profile;
    vProfile[2] = splice_sites_profile;
    vProfile[3] = start_stop_profile;
    return vProfile;
}

// Calculating similarity matrix
void alignNeedlemanForTranscript_simd_avx2int32::calculate_similarity_not_very_good( __m256i** vProfile, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix){
    //std::cout << "alignNeedlemanForTranscript_simd_avx2int32 line 229" << std::endl;
    const __m256i vZero = _mm256_set1_epi32(0);
    const __m256i minimum = _mm256_set1_epi32(-1147483648);
    //std::cout << "alignNeedlemanForTranscript_simd_avx2int32 line 231" << std::endl;
    __m256i* pvHStore = (__m256i*) calloc(this->_segLen, sizeof(__m256i));
    __m256i* pvHLoad = (__m256i*) calloc(this->_segLen, sizeof(__m256i));
    __m256i* pvE = (__m256i*) calloc(this->_segLen, sizeof(__m256i));
    __m256i* pvEC = (__m256i*) calloc(this->_segLen, sizeof(__m256i));
    __m256i* pvFC = (__m256i*) calloc(this->_segLen, sizeof(__m256i));
    __m256i* pvMC = (__m256i*) calloc(this->_segLen, sizeof(__m256i));
    __m256i* pv;
    __m256i* vP;
    __m256i vF;
    __m256i e;
    __m256i vGapO, vGapE, next_vGapO, next_vGapE, vTemp, vTemp2;
    //std::cout << "alignNeedlemanForTranscript_simd_avx2int32 line 242" << std::endl;
    int z;
    int32_t i, j, k, cmp;
    int32_t * a;
    int32_t * b;
    int32_t * c;
    int32_t * d;
    //for outer loop
    __m256i vH;
    for (i = 0; i<_length_of_d; ++i) {//important
        //std::cout << "alignNeedlemanForTranscript_simd_avx2int32 line 251" << std::endl;
        vF = minimum; /* Initialize F value to 0.
							   Any errors to vH values will be corrected in the Lazy_F loop.
							 */
        a = ( int32_t*)(&vF);
        a[0] = _extend_gap_penalty[_ref_elements[i]];
        //std::cout << "alignNeedlemanForTranscript_simd_avx2int32 line 255" << std::endl;
        vP = vProfile[_ref_elements[i]] + ((nucleotideCodeSubstitutionMatrix.get_dna_acid_map(_dna_d[i]))*(this->_segLen)); // Correct part of the vProfile
        //std::cout << "alignNeedlemanForTranscript_simd_avx2int32 line 257" << std::endl;
        // INDEL begin vector
        vGapO = _mm256_set1_epi32(_open_gap_penalty[_ref_elements[i]]); // here I use epi32, because for the score matrix I used int32_t
        // INDEL extension vector
        vGapE = _mm256_set1_epi32(_extend_gap_penalty[_ref_elements[i]]);
        if( i == _length_of_d-1){
            next_vGapO=vGapO;
            next_vGapE=vGapE;
        }else{
            next_vGapO=_mm256_set1_epi32(_open_gap_penalty[_ref_elements[i+1]]);
            next_vGapE= _mm256_set1_epi32(_extend_gap_penalty[_ref_elements[i+1]]);
        }
        //std::cout << "alignNeedlemanForTranscript_simd_avx2int32 line 269" << std::endl;
        vH = _mm256_loadu_si256(pvHStore+(this->_segLen-1));// the last one
        //vH = pvHStore[_segLen - 1];
        //std::cout << "alignNeedlemanForTranscript_simd_avx2int32 line 271" << std::endl;
        vH = _mm256_alignr_epi8(vH, _mm256_permute2x128_si256(vH, vH, _MM_SHUFFLE(0, 0, 2, 0)), 16 - 4);

        //vMaxColumn = vZero; /* vMaxColumn is used to record the max values of column i. */
        pv = pvHLoad;
        pvHLoad = pvHStore;
        pvHStore = pv;
        //std::cout << "alignNeedlemanForTranscript_simd_avx2int32 line 280" << std::endl;
        // inner loop to process the query sequence
        for (j = 0; j < this->_segLen; ++j) {
            if( 0 == i ){
                e = vGapE;
            }else{
                e = _mm256_loadu_si256(pvE + j);
            }
            vH = _mm256_add_epi32(vH, _mm256_loadu_si256(vP+j)); // add the scoring profile to vH
                _mm256_storeu_si256(pvEC+j, e);
                _mm256_storeu_si256(pvMC+j, vH);
                _mm256_storeu_si256(pvFC+j, vF);
            vH = _mm256_max_epi32(vH, e);
            vH = _mm256_max_epi32(vH, vF);
            //vMaxColumn = _mm256_max_epi32(vMaxColumn, vH);

            /* Save vH values. */
            _mm256_storeu_si256(pvHStore + j, vH);
            /* Update vE value. */
            vH = _mm256_add_epi32(vH, next_vGapO); // this is for next column
            e = _mm256_add_epi32(e, next_vGapE);
            e = _mm256_max_epi32(e, vH);
            _mm256_storeu_si256(pvE + j, e);

            /* Update vF value. */
            vH = _mm256_loadu_si256(pvHStore + j);
            vH = _mm256_add_epi32(vH, vGapO);
            vF = _mm256_add_epi32(vF, vGapE);
            vF = _mm256_max_epi32(vF, vH);

            /* Load the next vH. */
            vH = _mm256_loadu_si256(pvHLoad + j);
        }
        /* Lazy_F loop: has been revised to disallow adjecent insertion and then deletion, so don't update E(i, j), learn from SWPS3 */
        /* reset pointers to the start of the saved data */
        // it turns out the original method could not give the F matrix
        // if you modify that, it is not faster than the standard way

        j=0;
        vF = _mm256_alignr_epi8(vF, _mm256_permute2x128_si256(vF, vF, _MM_SHUFFLE(0, 0, 2, 0)), 16 - 4); //shift 4 byte (int32)
        vH = _mm256_loadu_si256 (pvHStore + j);
        a = ( int32_t*)(&vF);
        b = ( int32_t*)(&vH);
        int matchScore = b[0];
        a[0] = _extend_gap_penalty[_ref_elements[i]] > matchScore ? _extend_gap_penalty[_ref_elements[i]] : matchScore;
        _mm256_storeu_si256(pvFC+j, vF);

        vTemp2 = _mm256_add_epi32 (vH, vGapO);
        vTemp = _mm256_sub_epi32 (vF, vTemp2);
//        vTemp = _mm256_cmpeq_epi32 (vTemp, vZero);
//        cmp  = _mm256_movemask_epi8 (vTemp);
        a = (int32_t *) (&vTemp);
        //std::cout << "cmp " << cmp << std::endl;
        while ( a[0]>0 || a[1]>0 || a[2]>0 || a[3]>0 || a[4]>0 || a[5]>0 || a[6]>0 || a[7]>0){
                //std::cout << "line 321 cmp " << cmp << std::endl;
            vH = _mm256_max_epi32 (vH, vF);
            _mm256_storeu_si256 (pvHStore + j, vH);

            vF = _mm256_add_epi32 (vF, vGapE);
            ++j;
            if (j >= this->_segLen){
                j = 0;
                vF = _mm256_alignr_epi8(vF, _mm256_permute2x128_si256(vF, vF, _MM_SHUFFLE(0, 0, 2, 0)), 16 - 4); //shift 4 byte (int32)
                b = ( int32_t*)(&vF);
                b[0] = _extend_gap_penalty[_ref_elements[i]] > matchScore ? _extend_gap_penalty[_ref_elements[i]] : matchScore;
            }
                _mm256_storeu_si256(pvFC+j, vF);
            vH = _mm256_loadu_si256 (pvHStore + j);
            vTemp2 = _mm256_add_epi32 (vH, vGapO);
            vTemp = _mm256_sub_epi32 (vF, vTemp2);

            vF = _mm256_add_epi32(vF, vGapE);
            vF = _mm256_max_epi32(vF, vTemp2);

            a = (int32_t *) (&vTemp);
        }

        //std::cout << "alignNeedlemanForTranscript_simd_avx2int32 line 320" << std::endl;
        for( j=0; j<this->_segLen; ++j ){
            a = ( int32_t*)(pvHStore+j);
            b = ( int32_t*)(pvMC+j);
            c = ( int32_t*)(pvEC+j);
            d = ( int32_t*)(pvFC+j);
            for( z=0; z<8; ++z ){
                if( (j + z*this->_segLen) < _length_of_q ){
                    //std::cout << j+z*this->_segLen+1 << " " << i+1 << " " << a[z] << " " << b[z] << " " << c[z] << std::endl;
                    _trace_match[j+z*this->_segLen+1][i+1] = ( a[z] == b[z]);
                    _trace_del[j+z*this->_segLen+1][i+1] = ( a[z] == c[z] );
                    _trace_ins[j+z*this->_segLen+1][i+1] = ( a[z] == d[z] );
                }
            }
        }
        int t;

        //update with track matrices _trace_ins begin
        for( j=0; j<_length_of_q; ++j ){
            if( _trace_ins[j+1][i+1] && !_trace_del[j+1][i+1] && !_trace_match[j+1][i+1] ){
                t=0;
                while( i>=t && j>=t && _trace_del[j-t][i-t] && _trace_ins[j-t][i-t] ){
                    _trace_del[j-t][i-t] = false;
                    ++t;
                }
                t=0;
                while( i>=t && j>=t && _trace_match[j-t][i-t] && _trace_ins[j-t][i-t] ){
                    _trace_match[j-t][i-t] = false;
                    ++t;
                }
            }
        }//update with track matrices _trace_ins end

        //update with track matrices _track_del begin
        for( j=0; j<_length_of_q; ++j ){
            if( _trace_del[j+1][i+1] && _trace_match[j+1][i+1]==false && !_trace_ins[j+1][i+1] ){
                t=0;
                while( i>=t && _trace_del[j+1][i-t] && _trace_match[j+1][i-t] ){
                    _trace_match[j+1][i-t] = false;
                    ++t;
                }
                t=0;
                while( i>=t && _trace_del[j+1][i-t] && _trace_ins[j+1][i-t] ){
                    _trace_ins[j+1][i-t] = false;
                    ++t;
                }
            }
        }
        //update with track matrices _track_del end


        //update with track matrices _track_match begin
        for( j=0; j<_length_of_q; ++j ){
            if( _trace_match[j+1][i+1] && !_trace_del[j+1][i+1] && !_trace_ins[j+1][i+1] ){
                t=0;
                while( i>=t && j>=t && _trace_del[j-t][i-t] && _trace_match[j-t][i-t] ){
                    _trace_del[j-t][i-t] = false;
                    ++t;
                }
                t=0;
                while( i>=t && j>=t && _trace_match[j-t][i-t] && _trace_ins[j-t][i-t] ){
                    _trace_ins[j-t][i-t] = false;
                    ++t;
                }
            }
        }
        //update with track matrices _track_match end
    }
    //std::cout << "alignNeedlemanForTranscript_simd_avx2int32 line 341" << std::endl;
    for( j=0; j<this->_segLen; ++j ){
        a = ( int32_t*)(pvHStore+j);
        for( z=0; z<8; ++z ){
            if( (j + z*this->_segLen) < _length_of_q ){
                _similarity_matrix[j+z*this->_segLen+1] = a[z];
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
    free(pvFC);
    //std::cout << "alignNeedlemanForTranscript_simd_avx2int32 line 426" << std::endl;
    free(vProfile[0]);
    free(vProfile[1]);
    free(vProfile[2]);
    free(vProfile[3]);
    //std::cout << "alignNeedlemanForTranscript_simd_avx2int32 line 429" << std::endl;
    delete [] vProfile;
    //std::cout << "alignNeedlemanForTranscript_simd_avx2int32 line 431" << std::endl;
}


// Calculating similarity matrix
void alignNeedlemanForTranscript_simd_avx2int32::calculate_similarity( __m256i** vProfile, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix){
    const __m256i vZero = _mm256_set1_epi32(0);
    const __m256i minimum = _mm256_set1_epi32(-1147483648);
    __m256i* pvHStore = (__m256i*) calloc(this->_segLen, sizeof(__m256i));
    __m256i* pvHLoad = (__m256i*) calloc(this->_segLen, sizeof(__m256i));
    __m256i* pvE = (__m256i*) calloc(this->_segLen, sizeof(__m256i));
    __m256i* pvEC = (__m256i*) calloc(this->_segLen, sizeof(__m256i));
    __m256i* pvMC = (__m256i*) calloc(this->_segLen, sizeof(__m256i));
    __m256i* pv;
    __m256i* vP;
    __m256i e;
    __m256i vGapO, vGapE, next_vGapO, next_vGapE, vTemp, vTemp2;
    int z;
    int32_t i, j, k, cmp;
    int32_t * a;
    int32_t * b;
    int32_t * c;
    int32_t * d;
    //for outer loop
    __m256i vH;
    for (i = 0; i<_length_of_d; ++i) {//important
        //std::cout << "alignNeedlemanForTranscript_simd_avx2int32 line 251" << std::endl;
        //std::cout << "alignNeedlemanForTranscript_simd_avx2int32 line 255" << std::endl;
        vP = vProfile[_ref_elements[i]] + ((nucleotideCodeSubstitutionMatrix.get_dna_acid_map(_dna_d[i]))*(this->_segLen)); // Correct part of the vProfile
        //std::cout << "alignNeedlemanForTranscript_simd_avx2int32 line 257" << std::endl;
        // INDEL begin vector
        vGapO = _mm256_set1_epi32(_open_gap_penalty[_ref_elements[i]]); // here I use epi32, because for the score matrix I used int32_t
        // INDEL extension vector
        vGapE = _mm256_set1_epi32(_extend_gap_penalty[_ref_elements[i]]);
        if( i == _length_of_d-1){
            next_vGapO=vGapO;
            next_vGapE=vGapE;
        }else{
            next_vGapO=_mm256_set1_epi32(_open_gap_penalty[_ref_elements[i+1]]);
            next_vGapE= _mm256_set1_epi32(_extend_gap_penalty[_ref_elements[i+1]]);
        }
        //std::cout << "alignNeedlemanForTranscript_simd_avx2int32 line 269" << std::endl;
        vH = _mm256_loadu_si256(pvHStore+(this->_segLen-1));// the last one
        //vH = pvHStore[_segLen - 1];
        //std::cout << "alignNeedlemanForTranscript_simd_avx2int32 line 271" << std::endl;
        vH = _mm256_alignr_epi8(vH, _mm256_permute2x128_si256(vH, vH, _MM_SHUFFLE(0, 0, 2, 0)), 16 - 4);

        //vMaxColumn = vZero; /* vMaxColumn is used to record the max values of column i. */
        pv = pvHLoad;
        pvHLoad = pvHStore;
        pvHStore = pv;
        //std::cout << "alignNeedlemanForTranscript_simd_avx2int32 line 280" << std::endl;
        // inner loop to process the query sequence
        for (j = 0; j < this->_segLen; ++j) {
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

        for( j=0; j<this->_segLen; ++j ) {
            a = (int32_t *) (pvHStore + j);
            b = (int32_t *) (pvMC + j);
            c = (int32_t *) (pvEC + j);
            for (z = 0; z < 8; ++z) {
                if ((j + z * this->_segLen) < _length_of_q) {
                    _similarity_matrix[j+z*this->_segLen+1] = a[z];
                    _trace_match[j + z * this->_segLen + 1][i + 1] = (a[z] == b[z]);
                    _trace_del[j + z * this->_segLen + 1][i + 1] = (a[z] == c[z]);
                }
            }
        }

        int lastOne = 1, t; // insertion
        int32_t deletionScore;
        for( j=0; j<_length_of_q; ++j ){
            if( lastOne == 1 ){
                deletionScore = _similarity_matrix[j] + _extend_gap_penalty[_ref_elements[i]];
                if( deletionScore > _similarity_matrix[j+1] ){
                    _similarity_matrix[j+1] = deletionScore;
                    _trace_match[j+1][i+1] = false;
                    _trace_del[j+1][i+1] = false;
                    _trace_ins[j+1][i+1] = true;

                    t=0;
                    while( j>=t && _trace_ins[j-t][i+1] && _trace_match[j-t][i+1] ){
                        _trace_match[j-t][i+1] = false;
                        ++t;
                    }
                    t=0;
                    while( j>=t && _trace_ins[j-t][i+1] && _trace_del[j-t][i+1] ){
                        _trace_del[j-t][i+1] = false;
                        ++t;
                    }
                } else if ( deletionScore == _similarity_matrix[j+1] ){
                    _trace_ins[j+1][i+1] = true;
                } else {
                    _trace_ins[j+1][i+1] = false;
                    lastOne=0;//MATCH OR MISMATCH
                }
            }else{
                deletionScore = _similarity_matrix[j] + _open_gap_penalty[_ref_elements[i]];
                if( deletionScore > _similarity_matrix[j+1] ){
                    _similarity_matrix[j+1] = deletionScore;
                    _trace_match[j+1][i+1] = false;
                    _trace_del[j+1][i+1] = false;
                    _trace_ins[j+1][i+1] = true;
                    lastOne = 1;
                } else if ( deletionScore == _similarity_matrix[j+1] ){
                    _trace_ins[j+1][i+1] = true;
                    lastOne = 1;
                } else {
                    _trace_ins[j+1][i+1] = false;
                }
            }
            // put the score into the score matrix end
        }

        //update with track matrices _track_del begin
        for( j=0; j<_length_of_q; ++j ){
            if( _trace_del[j+1][i+1] && _trace_match[j+1][i+1]==false && !_trace_ins[j+1][i+1] ){
                t=0;
                while( i>=t && _trace_del[j+1][i-t] && _trace_match[j+1][i-t] ){
                    _trace_match[j+1][i-t] = false;
                    ++t;
                }
                t=0;
                while( i>=t && _trace_del[j+1][i-t] && _trace_ins[j+1][i-t] ){
                    _trace_ins[j+1][i-t] = false;
                    ++t;
                }
            }
        }
        //update with track matrices _track_del end

        //update with track matrices _track_match begin
        for( j=0; j<_length_of_q; ++j ){
            if( _trace_match[j+1][i+1] && !_trace_del[j+1][i+1] && !_trace_ins[j+1][i+1] ){
                t=0;
                while( i>=t && j>=t && _trace_del[j-t][i-t] && _trace_match[j-t][i-t] ){
                    _trace_del[j-t][i-t] = false;
                    ++t;
                }
                t=0;
                while( i>=t && j>=t && _trace_match[j-t][i-t] && _trace_ins[j-t][i-t] ){
                    _trace_ins[j-t][i-t] = false;
                    ++t;
                }
            }
        }
        //update with track matrices _track_match end

        //updata vH begin
        for( j=0; j<this->_segLen; j++ ){
            a = (int32_t*)(pvHStore+j);
            for( z=0; z<8; ++z ){
                if( (j + z*this->_segLen) < _length_of_q ){
                    a[z] = _similarity_matrix[j+z*this->_segLen+1];
                }
            }
        }//updata vH end
    }
    //std::cout << "alignNeedlemanForTranscript_simd_avx2int32 line 341" << std::endl;
    for( j=0; j<this->_segLen; ++j ){
        a = ( int32_t*)(pvHStore+j);
        for( z=0; z<8; ++z ){
            if( (j + z*this->_segLen) < _length_of_q ){
                _similarity_matrix[j+z*this->_segLen+1] = a[z];
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
    std::cout << "alignNeedlemanForTranscript_simd_avx2int32 line 431" << std::endl;
}

char element_int_to_char( int element ){
    switch(element) {
        case  0: return 'E';
        case  1: return 'I';
        case  2: return 'S';
        case  3: return 'P';
    }
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
    int32_t highestScore = this->_similarity_matrix[i];
    for( k=_length_of_q; k >0; --k ){
        if( _similarity_matrix[k] > highestScore ){
            i = k;
            highestScore= _similarity_matrix[k];
        }
    }
    for( k=_length_of_q; k>i; --k ){
        SQ.push(_dna_q[k-1]);
        SD.push('-');
        SI.push('-');
    }
    while (i > 0 || j > 0) {
        char element = element_int_to_char(this->_ref_elements[j - 1]);
        if (i == 0) {
            SQ.push('-');
            SD.push(_dna_d[j - 1]);
            SI.push(element);
            --j;
        } else if (j == 0) {
            SQ.push(_dna_q[i - 1]);
            SD.push('-');
            SI.push('-');
            --i;
        } else {
            if (_trace_match[i][j]) {
                SQ.push(_dna_q[i - 1]);
                SD.push(_dna_d[j - 1]);
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

            else if (_trace_ins[i][j] ) {//insertion
                SQ.push(_dna_q[i - 1]);
                SD.push('-');
                SI.push('-');
                --i;
            } else {        // Going to S(i, j-1) //deletion
                SQ.push('-');
                SD.push(_dna_d[j - 1]);
                SI.push(element);
                --j;
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
