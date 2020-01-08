//
// Created by song on 8/5/18.
//

#include "alignSlidingWindow.h"





__m256i* query_profile_avx2_int32( const int & _segLen, const size_t & _length_of_q, const std::string & _dna_q,
                                   NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix, const int & _n){
    __m256i* vProfile = (__m256i*)malloc(_n * _segLen * sizeof(__m256i));
    int32_t a, i, j, k;
    for( a =0; a<_n; ++a ){
        int32_t* s = (int32_t*)(vProfile+a*_segLen);
        for( i=0; i<_segLen; ++i ){
            j = i;
            for( k=0; k<8; ++k ){
                *s++ = j>= _length_of_q ? 0 : nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[a][nucleotideCodeSubstitutionMatrix.get_dna_acid_map(_dna_q[j])];
                j += _segLen;
            }
        }
    }
    return vProfile;
}

__m256i* query_profile_avx2_int16( const int & _segLen, const size_t & _length_of_q, const std::string & _dna_q,
                                   NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix, const int & _n){
    __m256i* vProfile = (__m256i*)malloc(_n * _segLen * sizeof(__m256i));
    int16_t a, i, j, k;
    for( a =0; a<_n; ++a ){
        int16_t* s = (int16_t*)(vProfile+a*_segLen);
        for( i=0; i<_segLen; ++i ){
            j = i;
            for( k=0; k<16; ++k ){
                *s++ = j>= _length_of_q ? 0 : nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[a][nucleotideCodeSubstitutionMatrix.get_dna_acid_map(_dna_q[j])];
                j += _segLen;
            }
        }
    }
    return vProfile;
}
__m256i* query_profile_avx2_int8( const int & _segLen, const size_t & _length_of_q, const std::string & _dna_q,
                                  NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix, const int & _n){
    __m256i* vProfile = (__m256i*)malloc(_n * _segLen * sizeof(__m256i));
    int8_t a, i, j, k;
    for( a =0; a<_n; ++a ){
        int8_t* s = (int8_t*)(vProfile+a*_segLen);
        for( i=0; i<_segLen; ++i ){
            j = i;
            for( k=0; k<32; ++k ){
                *s++ = j>= _length_of_q ? 0 : nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[a][nucleotideCodeSubstitutionMatrix.get_dna_acid_map(_dna_q[j])];
                j += _segLen;
            }
        }
    }
    return vProfile;
}

void alignment_avx_int32(const std::string& _dna_q, const std::string& _dna_d, std::stack<char> & SQ, std::stack<char> & SD,
                         const int8_t & _open_gap_penalty, const int8_t & _extend_gap_penalty,
                         NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix){

    size_t _length_of_q = _dna_q.length();
    size_t _length_of_d= _dna_d.length();
    size_t _segLen = ( _length_of_q+ 7)/8;

    __m256i* vProfile = query_profile_avx2_int32(_segLen, _length_of_q, _dna_q, nucleotideCodeSubstitutionMatrix, 5 );

    __m256i* pvHStore = (__m256i*) calloc(_segLen, sizeof(__m256i));
    __m256i* pvHLoad = (__m256i*) calloc(_segLen, sizeof(__m256i));
    __m256i* pvE = (__m256i*) calloc(_segLen, sizeof(__m256i)); //vE is the score from left side
    __m256i* pv;
    __m256i* vP;

    int32_t i, j, z;
    int32_t deletionScore;
    int32_t ** _similarity_matrix = new int32_t*[_length_of_q+1];
    for (i = 0; i < (_length_of_q+1); ++i) {
        _similarity_matrix[i] = new int32_t[_length_of_d+1];
    }
    i=0;
    for ( j = 0; j <= _length_of_d; ++j) {
        _similarity_matrix[i][j] = 0;
    }
    j=0;
    for (i = 0; i < (_length_of_q+1); ++i) {
        _similarity_matrix[i][j] = 0;
    }

    // INDEL begin vector
    const __m256i vGapO = _mm256_set1_epi32(_open_gap_penalty);
    // here I use epi32, because for the score matrix I used int32_t
    // INDEL extension vector
    const __m256i vGapE = _mm256_set1_epi32(_extend_gap_penalty);

    for(i=0; i<_segLen; ++i){
        _mm256_storeu_si256(pvE+i, vGapE);
    }
    //for outer loop
    __m256i vE, vH;
    int32_t* a;

    for (i = 0; (i != _length_of_d); ++i) {//important
        vH = _mm256_loadu_si256(pvHStore+(_segLen-1));// the last one
        vH = _mm256_alignr_epi8(vH, _mm256_permute2x128_si256(vH, vH, _MM_SHUFFLE(0, 0, 2, 0)), 16 - 4);

        // Shift the value in vH left by 4 byte. one score matrix cell
        // _mm256_slli_si256 not not work here,
        // because it shifts 8-bit [byte] elements within a 128-bit lane of the source vector

        vP = vProfile + ((nucleotideCodeSubstitutionMatrix.get_dna_acid_map(_dna_d[i]))*(_segLen)); // Correct part of the vProfile

        pv = pvHLoad;
        pvHLoad = pvHStore;
        pvHStore = pv;

        // inner loop to process the query sequence
        for (j = 0; j < _segLen; ++j) {
            vH = _mm256_add_epi32(vH, _mm256_loadu_si256(vP+j)); // add the scoring profile to vH
            // Get max from vH, vE
            vE = _mm256_loadu_si256(pvE+j);
            vH = _mm256_max_epi32(vH, vE);

            // Save vH values.
            _mm256_storeu_si256(pvHStore + j, vH);

            // Update vE value. for next database base pair, here should be get the exact score as general
            // smith-waterman algorithm (while it does not)
            vH = _mm256_add_epi32(vH, vGapO);
            vE = _mm256_add_epi32(vE, vGapE);
            vE = _mm256_max_epi32(vE, vH);
            _mm256_storeu_si256(pvE+j, vE);

            // Load the next vH.
            vH = _mm256_loadu_si256(pvHLoad + j);
        }

        // put the score into the score matrix begin
        for( j=0; j<_segLen; ++j ){
            a = ( int32_t*)(pvHStore+j);
            for( z=0; z<8; ++z ){
                if( (j + z*_segLen) < _length_of_q ){
                    _similarity_matrix[j+z*_segLen+1][i+1] = a[z];
                }
            }
        }

        int lastOne = 1; // INDEL
        for( j=0; j<_length_of_q; ++j ){
            if( lastOne == 1 ){
                deletionScore = _similarity_matrix[j][i+1]+_extend_gap_penalty;
                if( deletionScore >= _similarity_matrix[j+1][i+1]  ){
                    _similarity_matrix[j+1][i+1] = deletionScore;
                }else{
                    lastOne=0;//MATCH OR MISMATCH
                }
            }else{
                deletionScore = _similarity_matrix[j][i+1]+_open_gap_penalty;
                if( deletionScore >= _similarity_matrix[j+1][i+1]  ){
                    _similarity_matrix[j+1][i+1] = deletionScore;
                    lastOne=1;
                }else{
                    lastOne=0;//MATCH OR MISMATCH
                }
            }
        }
        // put the score into the score matrix end

        //updata vH begin
        for( j=0; j<_segLen; j++ ){
            a = (int32_t*)(pvHStore+j);
            for( z=0; z<8; ++z ){
                if( (j + z*_segLen) < _length_of_q ){
                    a[z] = _similarity_matrix[j+z*_segLen+1][i+1];
                }
            }
        }//updata vH end
    }
    free(pvE);
    free(pvHStore);
    free(vProfile);
    free(pvHLoad);

    //track step begin
    // this function works only for this purpose, since here I trimmed the sequences
    int index_i;
    int index_j;
    int S;

    size_t k;
    i = _length_of_q;
    j = _length_of_d;

    int32_t highestScore = _similarity_matrix[i][j];
    for( k=_length_of_q; k >0; --k ){
        if( _similarity_matrix[k][j] > highestScore ){
            i = k;
            highestScore= _similarity_matrix[k][j];
        }
    }

    for( k=_length_of_d; k >0; --k ){
        if( _similarity_matrix[_length_of_q][k] > highestScore ){
            j = k;
            i = _length_of_q;
            highestScore= _similarity_matrix[_length_of_q][k];
        }
    }
    // here, I do not track from the last elements of the matrix, so this is a kind of trim
    // This is not really global sequence alignment

    while (i != 0 || j != 0) {
        if (i == 0) {
            SQ.push('-');
            SD.push(_dna_d[j-1]);
            --j;
        } else if (j == 0) {
            SQ.push(_dna_q[i-1]);
            SD.push('-');
            --i;
        } else {
            index_i = nucleotideCodeSubstitutionMatrix.get_dna_acid_map(_dna_q[i-1]);
            index_j = nucleotideCodeSubstitutionMatrix.get_dna_acid_map(_dna_d[j-1]);
            S = nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[index_i][index_j];
            //  (_dna_q[ii-1] == _dna_d[jj-1]) ? match_score : -mismatch_score;
            if ( _similarity_matrix[i][j] == (_similarity_matrix[i-1][j-1] + S) ) {
                SQ.push(_dna_q[i-1]);
                SD.push(_dna_d[j-1]);
                --i; --j;
            } else if (_similarity_matrix[i-1][j] >= _similarity_matrix[i][j-1]) {
                SQ.push(_dna_q[i-1]);
                SD.push('-');
                --i;
            } else {
                SQ.push('-');
                SD.push(_dna_d[j-1]);
                --j;
            }
        }
    }//track end
    for ( i = 0; i < (_length_of_q+1); i++ ) {
        delete [] _similarity_matrix[i];
    }
    delete [] _similarity_matrix;
}

void alignment_avx_int16(const std::string& _dna_q, const std::string& _dna_d, std::stack<char> & SQ, std::stack<char> & SD,
                         const int8_t & _open_gap_penalty, const int8_t & _extend_gap_penalty,
                         NucleotideCodeSubstitutionMatrix & nucleotideCodeSubstitutionMatrix){

    size_t _length_of_q = _dna_q.length();
    size_t _length_of_d= _dna_d.length();
    size_t _segLen = ( _length_of_q+ 15)/16;

    __m256i* vProfile = query_profile_avx2_int16(_segLen, _length_of_q, _dna_q, nucleotideCodeSubstitutionMatrix, 5 );

    __m256i* pvHStore = (__m256i*) calloc(_segLen, sizeof(__m256i));
    __m256i* pvHLoad = (__m256i*) calloc(_segLen, sizeof(__m256i));
    __m256i* pvE = (__m256i*) calloc(_segLen, sizeof(__m256i)); //vE is the score from left side
    __m256i* pv;
    __m256i* vP;

    int16_t i, j, z;
    int16_t deletionScore;
    int16_t ** _similarity_matrix = new int16_t*[_length_of_q+1];
    for (i = 0; i < (_length_of_q+1); ++i) {
        _similarity_matrix[i] = new int16_t[_length_of_d+1];
    }
    i=0;
    for ( j = 0; j <= _length_of_d; ++j) {
        _similarity_matrix[i][j] = 0;
    }
    j=0;
    for (i = 0; i < (_length_of_q+1); ++i) {
        _similarity_matrix[i][j] = 0;
    }

    // INDEL begin vector
    const __m256i vGapO = _mm256_set1_epi16(_open_gap_penalty);
    // here I use epi16, because for the score matrix I used int16_t
    // INDEL extension vector
    const __m256i vGapE = _mm256_set1_epi16(_extend_gap_penalty);

    for(i=0; i<_segLen; ++i){
        _mm256_storeu_si256(pvE+i, vGapE);
    }
    //for outer loop
    __m256i vE, vH;
    int16_t* a;

    for (i = 0; (i != _length_of_d); ++i) {//important
        vH = _mm256_loadu_si256(pvHStore+(_segLen-1));// the last one
        vH = _mm256_alignr_epi8(vH, _mm256_permute2x128_si256(vH, vH, _MM_SHUFFLE(0, 0, 2, 0)), 16 - 2);

        // Shift the value in vH left by 4 byte. one score matrix cell
        // _mm256_slli_si256 not not work here,
        // because it shifts 8-bit [byte] elements within a 128-bit lane of the source vector

        vP = vProfile + ((nucleotideCodeSubstitutionMatrix.get_dna_acid_map(_dna_d[i]))*(_segLen)); // Correct part of the vProfile

        pv = pvHLoad;
        pvHLoad = pvHStore;
        pvHStore = pv;

        // inner loop to process the query sequence
        for (j = 0; j < _segLen; ++j) {
            vH = _mm256_add_epi16(vH, _mm256_loadu_si256(vP+j)); // add the scoring profile to vH
            // Get max from vH, vE
            vE = _mm256_loadu_si256(pvE+j);
            vH = _mm256_max_epi16(vH, vE);

            // Save vH values.
            _mm256_storeu_si256(pvHStore + j, vH);

            // Update vE value. for next database base pair, here should be get the exact score as general
            // smith-waterman algorithm (while it does not)
            vH = _mm256_add_epi16(vH, vGapO);
            vE = _mm256_add_epi16(vE, vGapE);
            vE = _mm256_max_epi16(vE, vH);
            _mm256_storeu_si256(pvE+j, vE);

            // Load the next vH.
            vH = _mm256_loadu_si256(pvHLoad + j);
        }

        // put the score into the score matrix begin
        for( j=0; j<_segLen; ++j ){
            a = ( int16_t*)(pvHStore+j);
            for( z=0; z<16; ++z ){
                if( (j + z*_segLen) < _length_of_q ){
                    _similarity_matrix[j+z*_segLen+1][i+1] = a[z];
                }
            }
        }

        int lastOne = 1; // INDEL
        for( j=0; j<_length_of_q; ++j ){
            if( lastOne == 1 ){
                deletionScore = _similarity_matrix[j][i+1]+_extend_gap_penalty;
                if( deletionScore >= _similarity_matrix[j+1][i+1]  ){
                    _similarity_matrix[j+1][i+1] = deletionScore;
                }else{
                    lastOne=0;//MATCH OR MISMATCH
                }
            }else{
                deletionScore = _similarity_matrix[j][i+1]+_open_gap_penalty;
                if( deletionScore >= _similarity_matrix[j+1][i+1]  ){
                    _similarity_matrix[j+1][i+1] = deletionScore;
                    lastOne=1;
                }else{
                    lastOne=0;//MATCH OR MISMATCH
                }
            }
        }
        // put the score into the score matrix end

        //updata vH begin
        for( j=0; j<_segLen; j++ ){
            a = (int16_t*)(pvHStore+j);
            for( z=0; z<16; ++z ){
                if( (j + z*_segLen) < _length_of_q ){
                    a[z] = _similarity_matrix[j+z*_segLen+1][i+1];
                }
            }
        }//updata vH end
    }
    free(pvE);
    free(pvHStore);
    free(vProfile);
    free(pvHLoad);

    //track step begin
    // this function works only for this purpose, since here I trimmed the sequences
    int index_i;
    int index_j;
    int S;

    size_t k;
    i = _length_of_q;
    j = _length_of_d;

    int16_t highestScore = _similarity_matrix[i][j];
    for( k=_length_of_q; k >0; --k ){
        if( _similarity_matrix[k][j] > highestScore ){
            i = k;
            highestScore= _similarity_matrix[k][j];
        }
    }

    for( k=_length_of_d; k >0; --k ){
        if( _similarity_matrix[_length_of_q][k] > highestScore ){
            j = k;
            i = _length_of_q;
            highestScore= _similarity_matrix[_length_of_q][k];
        }
    }
    // here, I do not track from the last elements of the matrix, so this is a kind of trim
    // This is not really global sequence alignment
    //std::cout << "i " << i << " j " << j << " highestScore " << highestScore << std::endl;

    while (i != 0 || j != 0) {
        if (i == 0) {
            SQ.push('-');
            SD.push(_dna_d[j-1]);
            --j;
        } else if (j == 0) {
            SQ.push(_dna_q[i-1]);
            SD.push('-');
            --i;
        } else {
            index_i = nucleotideCodeSubstitutionMatrix.get_dna_acid_map(_dna_q[i-1]);
            index_j = nucleotideCodeSubstitutionMatrix.get_dna_acid_map(_dna_d[j-1]);
            S = nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[index_i][index_j];
            //  (_dna_q[ii-1] == _dna_d[jj-1]) ? match_score : -mismatch_score;
            if ( _similarity_matrix[i][j] == (_similarity_matrix[i-1][j-1] + S) ) {
                SQ.push(_dna_q[i-1]);
                SD.push(_dna_d[j-1]);
                --i; --j;
            } else if (_similarity_matrix[i-1][j] >= _similarity_matrix[i][j-1]) {
                SQ.push(_dna_q[i-1]);
                SD.push('-');
                --i;
            } else {
                SQ.push('-');
                SD.push(_dna_d[j-1]);
                --j;
            }
        }
    }//track end
    for ( i = 0; i < (_length_of_q+1); i++ ) {
        delete [] _similarity_matrix[i];
    }
    delete [] _similarity_matrix;
}

void alignment_avx_int8(const std::string& _dna_q, const std::string& _dna_d, std::stack<char> & SQ, std::stack<char> & SD,
                         const int8_t & _open_gap_penalty, const int8_t & _extend_gap_penalty,
                        NucleotideCodeSubstitutionMatrix & nucleotideCodeSubstitutionMatrix){

    size_t _length_of_q = _dna_q.length();
    size_t _length_of_d= _dna_d.length();
    size_t _segLen = ( _length_of_q+31)/32;
    __m256i* vProfile = query_profile_avx2_int8(_segLen, _length_of_q, _dna_q, nucleotideCodeSubstitutionMatrix, 5 );
    __m256i* pvHStore = (__m256i*) calloc(_segLen, sizeof(__m256i));
    __m256i* pvHLoad = (__m256i*) calloc(_segLen, sizeof(__m256i));
    __m256i* pvE = (__m256i*) calloc(_segLen, sizeof(__m256i)); //vE is the score from left side
    __m256i* pv;
    __m256i* vP;

    int8_t i, j, z;

    int8_t ** _similarity_matrix = new int8_t*[_length_of_q+1];
    for (i = 0; i < (_length_of_q+1); ++i) {
        _similarity_matrix[i] = new int8_t[_length_of_d+1];
    }
    i=0;
    for ( j = 0; j <= _length_of_d; ++j) {
        _similarity_matrix[i][j] = 0;
    }
    j=0;
    for (i = 0; i < (_length_of_q+1); ++i) {
        _similarity_matrix[i][j] = 0;
    }

    // INDEL begin vector
    const __m256i vGapO = _mm256_set1_epi8(_open_gap_penalty);
    // here I use epi8, because for the score matrix I used int8_t
    // INDEL extension vector
    const __m256i vGapE = _mm256_set1_epi8(_extend_gap_penalty);

    for(i=0; i<_segLen; ++i){
        _mm256_storeu_si256(pvE+i, vGapE);
    }
    //for outer loop
    __m256i vE, vH;
    int8_t* a;
    int8_t deletionScore;
    for (i = 0; (i != _length_of_d); ++i) {//important
        vH = _mm256_loadu_si256(pvHStore+(_segLen-1));// the last one
        vH = _mm256_alignr_epi8(vH, _mm256_permute2x128_si256(vH, vH, _MM_SHUFFLE(0, 0, 2, 0)), 16 - 1);

        // Shift the value in vH left by 4 byte. one score matrix cell
        // _mm256_slli_si256 not not work here,
        // because it shifts 8-bit [byte] elements within a 128-bit lane of the source vector

        vP = vProfile + ((nucleotideCodeSubstitutionMatrix.get_dna_acid_map(_dna_d[i]))*(_segLen)); // Correct part of the vProfile

        pv = pvHLoad;
        pvHLoad = pvHStore;
        pvHStore = pv;
        // inner loop to process the query sequence
        for (j = 0; j < _segLen; ++j) {
            vH = _mm256_add_epi8(vH, _mm256_loadu_si256(vP+j)); // add the scoring profile to vH
            // Get max from vH, vE
            vE = _mm256_loadu_si256(pvE+j);
            vH = _mm256_max_epi8(vH, vE);
            // Save vH values.
            _mm256_storeu_si256(pvHStore + j, vH);

            // Update vE value. for next database base pair, here should be get the exact score as general
            // smith-waterman algorithm (while it does not)
            vH = _mm256_add_epi8(vH, vGapO);
            vE = _mm256_add_epi8(vE, vGapE);
            vE = _mm256_max_epi8(vE, vH);
            _mm256_storeu_si256(pvE+j, vE);

            // Load the next vH.
            vH = _mm256_loadu_si256(pvHLoad + j);
        }
        // put the score into the score matrix begin
        for( j=0; j<_segLen; ++j ){
            a = ( int8_t*)(pvHStore+j);
            for( z=0; z<32; ++z ){
                if( (j + z*_segLen) < _length_of_q ){
                    _similarity_matrix[j+z*_segLen+1][i+1] = a[z];
                }
            }
        }

        int lastOne = 1; // INDEL
        for( j=0; j<_length_of_q; ++j ){
            if( lastOne == 1 ){
                deletionScore = _similarity_matrix[j][i+1]+_extend_gap_penalty;
                if( deletionScore >= _similarity_matrix[j+1][i+1]  ){
                    _similarity_matrix[j+1][i+1] = deletionScore;
                }else{
                    lastOne=0;//MATCH OR MISMATCH
                }
            }else{
                deletionScore = _similarity_matrix[j][i+1]+_open_gap_penalty;
                if( deletionScore >= _similarity_matrix[j+1][i+1]  ){
                    _similarity_matrix[j+1][i+1] = deletionScore;
                    lastOne=1;
                }else{
                    lastOne=0;//MATCH OR MISMATCH
                }
            }
        }
        // put the score into the score matrix end

        //updata vH begin
        for( j=0; j<_segLen; j++ ){
            a = (int8_t*)(pvHStore+j);
            for( z=0; z<32; ++z ){
                if( (j + z*_segLen) < _length_of_q ){
                    a[z] = _similarity_matrix[j+z*_segLen+1][i+1];
                }
            }
        }//updata vH end
    }
    free(pvE);
    free(pvHStore);
    free(vProfile);
    free(pvHLoad);

    //track step begin
    // this function works only for this purpose, since here I trimmed the sequences
    int index_i;
    int index_j;
    int S;

    size_t k;
    i = _length_of_q;
    j = _length_of_d;

    int8_t highestScore = _similarity_matrix[i][j];
    for( k=_length_of_q; k >0; --k ){
        if( _similarity_matrix[k][j] > highestScore ){
            i = k;
            highestScore= _similarity_matrix[k][j];
        }
    }

    for( k=_length_of_d; k >0; --k ){
        if( _similarity_matrix[_length_of_q][k] > highestScore ){
            j = k;
            i = _length_of_q;
            highestScore= _similarity_matrix[_length_of_q][k];
        }
    }
    // here, I do not track from the last elements of the matrix, so this is a kind of trim
    // This is not really global sequence alignment

    while (i != 0 || j != 0) {
        if (i == 0) {
            SQ.push('-');
            SD.push(_dna_d[j-1]);
            --j;
        } else if (j == 0) {
            SQ.push(_dna_q[i-1]);
            SD.push('-');
            --i;
        } else {
            index_i = nucleotideCodeSubstitutionMatrix.get_dna_acid_map(_dna_q[i-1]);
            index_j = nucleotideCodeSubstitutionMatrix.get_dna_acid_map(_dna_d[j-1]);
            S = nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[index_i][index_j];
            //  (_dna_q[ii-1] == _dna_d[jj-1]) ? match_score : -mismatch_score;
            if ( _similarity_matrix[i][j] == (_similarity_matrix[i-1][j-1] + S) ) {
                SQ.push(_dna_q[i-1]);
                SD.push(_dna_d[j-1]);
                --i; --j;
            } else if (_similarity_matrix[i-1][j] >= _similarity_matrix[i][j-1]) {
                SQ.push(_dna_q[i-1]);
                SD.push('-');
                --i;
            } else {
                SQ.push('-');
                SD.push(_dna_d[j-1]);
                --j;
            }
        }
    }//track end
    for ( i = 0; i < (_length_of_q+1); i++ ) {
        delete [] _similarity_matrix[i];
    }
    delete [] _similarity_matrix;
}

//void alignSlidingWindow_old_version( const std::string& dna_q, const std::string& dna_d,
void alignSlidingWindow( const std::string& dna_q, const std::string& dna_d,
        std::string & _alignment_q, std::string & _alignment_d, const int & slidingWindowSize,
        std::map<std::string, std::string>& parameters, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix ){
////    std::cout << slidingWindowSize << " " << dna_q.length() << " " << dna_d.length() << std::endl;
//    int8_t substitute_matrix[5][5] = {
//            //  A,     T,     C,     G,     N
//            {    1,    -1,    -1,    -1,    0} ,   //A
//            {   -1,     1,    -1,    -1,    0} ,   //T
//            {   -1,    -1,     1,    -1,    0} ,   //C
//            {   -1,    -1,    -1,     1,    0},     //G
//            {    0,     0,     0,     0,    0}     //N
//    };
//    int8_t ** _substitute_matrix = new int8_t*[5];
;
//    for( i=0; i<5;++i ){
//        _substitute_matrix[i] = new int8_t[5];
//        for( j=0; j<5;++j ){
//            _substitute_matrix[i][j] = substitute_matrix[i][j];
//        }
//    }
//
//    std::map<char , int8_t> _dna_acid_map;
//    _dna_acid_map['A']=0;
//    _dna_acid_map['a']=0;
//    _dna_acid_map['T']=1;
//    _dna_acid_map['t']=1;
//    _dna_acid_map['C']=2;
//    _dna_acid_map['c']=2;
//    _dna_acid_map['G']=3;
//    _dna_acid_map['g']=3;
    int8_t _open_gap_penalty=stoi(get_parameters("alignmentOpenGapP", parameters));
    int8_t _extend_gap_penalty=stoi(get_parameters("alignmentExtendGapP", parameters));
    size_t _length_of_q=dna_q.size();
    size_t _length_of_d=dna_d.size();

    //2^15 = 32768
    //of the maximum length of the windowSize of is about 32000/2 = 16000
    size_t databaseStart=1;
    size_t databaseEnd = 0;
    size_t queryStart=1;
    size_t queryEnd = 0;
    while( databaseStart<_length_of_d && queryStart< _length_of_q){
        databaseEnd=databaseStart+slidingWindowSize;
        queryEnd=queryStart+slidingWindowSize;
        if( databaseEnd>_length_of_d ){
            databaseEnd=_length_of_d;
        }
        if( queryEnd>_length_of_q ){
            queryEnd=_length_of_q;
        }
        std::string qSeq = getSubsequence(dna_q, queryStart, queryEnd );
        std::string dSeq = getSubsequence(dna_d, databaseStart, databaseEnd );
//        std::cout << qSeq << std::endl << dSeq << std::endl;
        std::stack<char> SQ;
        std::stack<char> SD;
        if( slidingWindowSize>1073741824 ){
            std::cout << "the windows size is too large" << std::endl;
            exit(1);
        }else if( slidingWindowSize>16384 ){
            alignment_avx_int32(qSeq, dSeq, SQ, SD, _open_gap_penalty, _extend_gap_penalty, nucleotideCodeSubstitutionMatrix);
        } else if( slidingWindowSize>65 ){
            alignment_avx_int16(qSeq, dSeq, SQ, SD, _open_gap_penalty, _extend_gap_penalty, nucleotideCodeSubstitutionMatrix);
        }else{
            alignment_avx_int8(qSeq, dSeq, SQ, SD, _open_gap_penalty, _extend_gap_penalty, nucleotideCodeSubstitutionMatrix); // this function is not correct do not run it
        }
        while (!SQ.empty()) {
            _alignment_q += SQ.top();
            _alignment_d += SD.top();
            if( SQ.top() != '-' ){
                ++queryStart;
            }
            if( SD.top() != '-' ){
                ++databaseStart;
            }
            SQ.pop();
            SD.pop();
        }
//        std::cout << _alignment_q << std::endl << _alignment_d << std::endl;
    }
    while( databaseStart<=_length_of_d ){
        _alignment_q += '-';
        _alignment_d += dna_d[databaseStart-1];
        ++databaseStart;
    }
    while( queryStart<=_length_of_q ){
        _alignment_q += dna_q[queryStart-1];
        _alignment_d += '-';
        ++queryStart;
    }
//    std::cout << "sequence alignment done" << std::endl;
}



//void alignSlidingWindow_old_version( const std::string& dna_q, const std::string& dna_d,
void alignSlidingWindow( const std::string& dna_q, const std::string& dna_d,
                         std::string & _alignment_q, std::string & _alignment_d, const int & slidingWindowSize,
                         const int & startShiftDistance, const int & endShiftDistance , std::map<std::string, std::string>& parameters,
                         NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix ){
////    std::cout << slidingWindowSize << " " << dna_q.length() << " " << dna_d.length() << std::endl;
//    int8_t substitute_matrix[5][5] = {
//            //  A,     T,     C,     G,     N
//            {    1,    -1,    -1,    -1,    0} ,   //A
//            {   -1,     1,    -1,    -1,    0} ,   //T
//            {   -1,    -1,     1,    -1,    0} ,   //C
//            {   -1,    -1,    -1,     1,    0},     //G
//            {    0,     0,     0,     0,    0}     //N
//    };
//    int8_t ** _substitute_matrix = new int8_t*[5];
//    for( i=0; i<5;++i ){
//        _substitute_matrix[i] = new int8_t[5];
//        for( j=0; j<5;++j ){
//            _substitute_matrix[i][j] = substitute_matrix[i][j];
//        }
//    }
//
//    std::map<char , int8_t> _dna_acid_map;
//    _dna_acid_map['A']=0;
//    _dna_acid_map['a']=0;
//    _dna_acid_map['T']=1;
//    _dna_acid_map['t']=1;
//    _dna_acid_map['C']=2;
//    _dna_acid_map['c']=2;
//    _dna_acid_map['G']=3;
//    _dna_acid_map['g']=3;
    int shiftdistance = startShiftDistance;
    if( shiftdistance < endShiftDistance ){
        shiftdistance = endShiftDistance;
    }

    int8_t _open_gap_penalty=stoi(get_parameters("alignmentOpenGapP", parameters));
    int8_t _extend_gap_penalty=stoi(get_parameters("alignmentExtendGapP", parameters));
    size_t _length_of_q=dna_q.size();
    size_t _length_of_d=dna_d.size();

    //2^15 = 32768
    //of the maximum length of the windowSize of is about 32000/2 = 16000
    size_t databaseStart=1;
    size_t databaseEnd = 0;
    size_t queryStart=1;
    size_t queryEnd = 0;
    int thisSlidingWindowSize = slidingWindowSize + shiftdistance; //give the first one a large window size
    while( databaseStart<_length_of_d && queryStart< _length_of_q){
        databaseEnd=databaseStart+thisSlidingWindowSize;
        queryEnd=queryStart+thisSlidingWindowSize;
        if( databaseEnd>_length_of_d ){
            databaseEnd=_length_of_d;
        }
        if( queryEnd>_length_of_q ){
            queryEnd=_length_of_q;
        }
        std::string qSeq = getSubsequence(dna_q, queryStart, queryEnd );
        std::string dSeq = getSubsequence(dna_d, databaseStart, databaseEnd );
//        std::cout << qSeq << std::endl << dSeq << std::endl;
        std::stack<char> SQ;
        std::stack<char> SD;
        if( thisSlidingWindowSize>1073741824 ){
            std::cout << "the windows size is too large" << std::endl;
            exit(1);
        }else if( thisSlidingWindowSize>16384 ){
            alignment_avx_int32(qSeq, dSeq, SQ, SD, _open_gap_penalty, _extend_gap_penalty, nucleotideCodeSubstitutionMatrix);
        } else if( thisSlidingWindowSize>65 ){
            alignment_avx_int16(qSeq, dSeq, SQ, SD, _open_gap_penalty, _extend_gap_penalty, nucleotideCodeSubstitutionMatrix);
        }else{
            alignment_avx_int8(qSeq, dSeq, SQ, SD, _open_gap_penalty, _extend_gap_penalty, nucleotideCodeSubstitutionMatrix); // this function is not correct do not run it
        }
        while (!SQ.empty()) {
            _alignment_q += SQ.top();
            _alignment_d += SD.top();
            if( SQ.top() != '-' ){
                ++queryStart;
            }
            if( SD.top() != '-' ){
                ++databaseStart;
            }
            SQ.pop();
            SD.pop();
        }
        thisSlidingWindowSize=slidingWindowSize;
//        std::cout << _alignment_q << std::endl << _alignment_d << std::endl;
    }
    while( databaseStart<=_length_of_d ){
        _alignment_q += '-';
        _alignment_d += dna_d[databaseStart-1];
        ++databaseStart;
    }
    while( queryStart<=_length_of_q ){
        _alignment_q += dna_q[queryStart-1];
        _alignment_d += '-';
        ++queryStart;
    }
//    std::cout << "sequence alignment done" << std::endl;
}



static const int32_t SCORE_OUT_BANDED_ALIGNMENT_REGION = -536870912;


void bandAlign( const std::string& dna_q, const std::string& dna_d,
                         std::string & _alignment_q, std::string & _alignment_d, const int & slidingWindowSize,
                         const int & startShiftDistance, const int & endShiftDistance , std::map<std::string, std::string>& parameters,
                         NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix ){
    int sw=slidingWindowSize;
    sw = slidingWindowSize>startShiftDistance?slidingWindowSize:startShiftDistance;
    sw = slidingWindowSize>endShiftDistance?slidingWindowSize:endShiftDistance;

    bandAlign( dna_q, dna_d, _alignment_q, _alignment_d, sw, parameters, nucleotideCodeSubstitutionMatrix );
}
void bandAlign( const std::string& dna_q, const std::string& dna_d,
                         std::string & _alignment_q, std::string & _alignment_d, const int & slidingWindowSize,
                         std::map<std::string, std::string>& parameters,
                         NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix ){

    int8_t _open_gap_penalty=stoi(get_parameters("alignmentOpenGapP", parameters));
    int8_t _extend_gap_penalty=stoi(get_parameters("alignmentExtendGapP", parameters));

    int length1 = dna_d.length();
    int length2 = dna_q.length();
    int32_t i=0, j=0, mscore=0;

    uint8_t d=0;

    Matrix T(length1+1, length2+1);

    int32_t * t;
    int32_t * M1 = new int32_t [length2 + 1]; //M1 and M2 is for the previous column and the current column
    int32_t * M2 = new int32_t [length2 + 1];

    int32_t * F = new int32_t [length2 + 1];

    std::fill_n(M1, length2+1, SCORE_OUT_BANDED_ALIGNMENT_REGION);
    std::fill_n(M2, length2+1, SCORE_OUT_BANDED_ALIGNMENT_REGION);
    std::fill_n(F, length2+1, SCORE_OUT_BANDED_ALIGNMENT_REGION);
    M1[0]=0;
    M2[0]=0;

    int32_t e;
    int32_t  start, end, h, f;
    for ( i=1; i<=length1; ++i ){
        start = i - slidingWindowSize;
        if ( start < 1 ){
            start = 1;
        }
        end = i + slidingWindowSize;
        if( end > length2 ){
            end = length2;
        }
        e = SCORE_OUT_BANDED_ALIGNMENT_REGION;
        F[end] = SCORE_OUT_BANDED_ALIGNMENT_REGION;
        for (j = start; j <= end; ++j) {
            mscore = M1[j - 1] + nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[nucleotideCodeSubstitutionMatrix.get_dna_acid_map(dna_d[i-1])][nucleotideCodeSubstitutionMatrix.get_dna_acid_map(dna_q[j-1])];

            d = mscore > F[j] ? 0 : 1;
            mscore = mscore > F[j] ? mscore : F[j];

            d = mscore > e ? d : 2;
            mscore = mscore > e ? mscore : e;

            M2[j] = mscore;

            h = mscore + _open_gap_penalty;
            f = F[j] + _extend_gap_penalty;
            d |= f >= h ? 1 << 3 : 0;
            f = f >= h ? f : h;
            F[j] = f;

            e += _extend_gap_penalty;
            d |= e >= h ? 1 << 4 : 0;
            e = e >= h ? e : h;

            T.set(i, j, d);
        }

        t = M1; // this pointer switch, should be very fast
        M1 = M2;
        M2 = t;
    }
    std::vector<uint32_t> cigar;

    uint32_t op = 0;
    uint32_t length = 1;
    // trace back begin
    // For speed up and RAM saving purpose, this is just am approximation tracing back implementation
    int32_t ii = length1;
    int32_t jj = length2;
    int tmp, state=0;
    while (ii>0 && jj>0) {
        tmp = T.get(ii, jj);
        if (state == 0) state = tmp & 7; // if requesting the H state, find state one maximizes it.
        else if (!(tmp >> (state + 2) & 1)) state = 0; // if requesting other states, _state_ stays the same if it is a continuation; otherwise, set to H
        if (state == 0) state = tmp & 7;
        if (state == 0){
            op=0;
            --ii;
            --jj;
        }else if (state == 1 || state == 3){
            op =2;
            --ii;
        }else{
            op =1;
            --jj;
        }

        if( cigar.size() == 0 || op != (cigar[cigar.size() - 1]&0xf) ){
            cigar.push_back(length << 4 | op);
        }else{
            cigar[cigar.size()-1] += length<<4;
        }
    }
    while( ii>0 ){
        op =2;
        --ii;
        if( cigar.size() == 0 || op != (cigar[cigar.size() - 1]&0xf) ){
            cigar.push_back(length << 4 | op);
        }else{
            cigar[cigar.size()-1] += length<<4;
        }
    }
    while( jj>0 ){
        op =1;
        --jj;
        if( cigar.size() == 0 || op != (cigar[cigar.size() - 1]&0xf) ){
            cigar.push_back(length << 4 | op);
        }else{
            cigar[cigar.size()-1] += length<<4;
        }
    }

    delete[] M1;
    delete[] M2;
    delete[] F;

    int32_t refPosition = 0;
    int32_t queryPostion = 0;

    for( i=cigar.size()-1; i >= 0; --i ){
        uint32_t cigarLength = cigar[i]>>4;
        uint32_t cigarType = cigar[i]&0xf;
        // "MID"
        if( 0 == cigarType ){
            _alignment_d += dna_d.substr( refPosition, cigarLength );
            _alignment_q += dna_q.substr( queryPostion, cigarLength );
            refPosition += cigarLength;
            queryPostion += cigarLength;
        } else if ( 1 == cigarType ){
            _alignment_d += std::string(cigarLength, '-');
            _alignment_q += dna_q.substr( queryPostion, cigarLength );
            queryPostion += cigarLength;
        }else if (2 == cigarType){
            _alignment_d += dna_d.substr( refPosition, cigarLength );
            _alignment_q += std::string(cigarLength, '-');
            refPosition += cigarLength;
        }else{
            std::cout << "some thing unknown happened with the cigar vector parser" << std::endl;
        }
    }
}

/*
 * currently, I still have some difficulty to use the minimap2 algorithm
 * **/
void alignSlidingWindow_new_version( const std::string& dna_q, const std::string& dna_d,
//void alignSlidingWindow( const std::string& dna_q, const std::string& dna_d,
                                     std::string & _alignment_q, std::string & _alignment_d, const int & slidingWindowSize,
                                     std::map<std::string, std::string>& parameters, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix ){
    //for hengAlignment function begin
    //mm_tbuf_t *tbuf = mm_tbuf_init();
    //ksw_extz_t ez;
    //memset(&ez, 0, sizeof(ksw_extz_t));
    //for hengAlignment end

    int8_t _open_gap_penalty=stoi(get_parameters("alignmentOpenGapP", parameters));
    int8_t _extend_gap_penalty=stoi(get_parameters("alignmentExtendGapP", parameters));
    size_t _length_of_q=dna_q.size();
    size_t _length_of_d=dna_d.size();

    size_t databaseStart=1;
    size_t databaseEnd = 0;
    size_t queryStart=1;
    size_t queryEnd = 0;
    while( databaseStart<_length_of_d && queryStart< _length_of_q){
        databaseEnd=databaseStart+slidingWindowSize;
        queryEnd=queryStart+slidingWindowSize;
        if( databaseEnd>_length_of_d ){
            databaseEnd=_length_of_d;
        }
        if( queryEnd>_length_of_q ){
            queryEnd=_length_of_q;
        }
        std::string qSeq = getSubsequence(dna_q, queryStart, queryEnd );
        std::string dSeq = getSubsequence(dna_d, databaseStart, databaseEnd );
        std::stack<char> SQ;
        std::stack<char> SD;
        //std::cout << "databaseStart " << databaseStart << " queryStart " << queryStart << std::endl;
        //hengAlign(qSeq, dSeq, SQ, SD, _open_gap_penalty, _extend_gap_penalty, nucleotideCodeSubstitutionMatrix, &(tbuf->km), &ez);
        hengAlign(qSeq, dSeq, SQ, SD, _open_gap_penalty, _extend_gap_penalty, nucleotideCodeSubstitutionMatrix);
        //std::cout << "databaseStart " << databaseStart << " queryStart " << queryStart << std::endl;
        while (!SQ.empty()) {
            _alignment_q += SQ.top();
            _alignment_d += SD.top();
            if( SQ.top() != '-' ){
                ++queryStart;
            }
            if( SD.top() != '-' ){
                ++databaseStart;
            }
            SQ.pop();
            SD.pop();
        }
    }
    while( databaseStart<=_length_of_d ){
        _alignment_q += '-';
        _alignment_d += dna_d[databaseStart-1];
        ++databaseStart;
    }
    while( queryStart<=_length_of_q ){
        _alignment_q += dna_q[queryStart-1];
        _alignment_d += '-';
        ++queryStart;
    }
    //mm_tbuf_destroy(tbuf);
}
