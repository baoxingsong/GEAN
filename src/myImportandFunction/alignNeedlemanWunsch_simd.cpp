//
// Created by Baoxing song on 25.04.18.
//
//
// this code is the needleman-wunsch implenmation with SSE/AVX technology
// for the first line and first column, the values are always 0, so it is not really needleman wunsch algorithm

//todo has not been well tested yet
#include "alignNeedlemanWunsch_simd.h"


NeedlemanWunsch_simd_Slow::NeedlemanWunsch_simd_Slow(const std::string & dna_q, const std::string & dna_d, const int8_t & match_score,
                                           const int8_t & mis_match_score, const int8_t & open_gap_penalty, const int8_t & extend_gap_penalty,  NucleotideCodeSubstitutionMatrix & _nucleotideCodeSubstitutionMatrix) {

    this->_dna_q = dna_q;
    this->_dna_d = dna_d;
    this->_open_gap_penalty = open_gap_penalty;
    this->_extend_gap_penalty = extend_gap_penalty;
    this->_length_of_q = this->_dna_q.length();
    this->_length_of_d = this->_dna_d.length();
    this->bias = 0;
    this->_n = 5;

    this->_segLen = ( _length_of_q+ 7)/8;
    //score matrix and track matrix begin

    this->_similarity_matrix = new int32_t*[this->_length_of_q+1];
    int i, j;
    for (i = 0; i < (this->_length_of_q+1); ++i) {
        this->_similarity_matrix[i] = new int32_t[this->_length_of_d+1];
    }
    i=0;
    for ( j = 0; j <= _length_of_d; ++j) {
        this->_similarity_matrix[i][j] = 0;
    }
    j=0;
    for (i = 0; i < (this->_length_of_q+1); ++i) {
        this->_similarity_matrix[i][j] = 0;
    }
    this->_track_match = new int32_t*[this->_length_of_q + 1];
    this->_track_del = new int32_t*[this->_length_of_q + 1];
    this->_track_ins = new int32_t*[this->_length_of_q + 1];
    for ( i = 0; i < (this->_length_of_q + 1); ++i) {
        _track_match[i] = new int32_t[this->_length_of_d + 1];
        _track_del[i] = new int32_t[this->_length_of_d + 1];
        _track_ins[i] = new int32_t[this->_length_of_d + 1];
    }
    for ( i = 0; i <= 1; i++) {
        for ( j = 0; j <= _length_of_d; ++j) {
            _track_match[i][j] = -1;
            _track_del[i][j] = -1;
            _track_ins[i][j] = -1;
        }
    }
    for ( j = 0; j <= 1; j++){
        for ( i = 0; i <= _length_of_q; ++i) {
            _track_match[i][j] = -1;
            _track_del[i][j] = -1;
            _track_ins[i][j] = -1;
        }
    }

    this->_substitute_matrix = new int32_t*[_n];
    for( i=0; i<_n; ++i){
        this->_substitute_matrix[i] = new int32_t[_n];
        for( j=0; j<_n; ++j){
            if( _nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[i][j] > 0 ){
                this->_substitute_matrix[i][j] = _nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[i][j] * match_score;
            }else{
                this->_substitute_matrix[i][j] = -(_nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[i][j] * mis_match_score);
            }
        }
    }

    // here let's assume we have sequence with length 10M, the score of exon match is 36, then the total score could be 3.6e+08
    // the value range of int32 is about -2e+09 -2e+09, so it should be enough
    //check potential overflow problem begin
    uint64_t maxLength = _length_of_d;
    if( maxLength <  _length_of_q){
        maxLength = _length_of_q;
    }
    int32_t maxScore = match_score;
    if( (maxLength * maxScore) > pow(2,31)){
        std::cerr << "the sequence is too long, could not be optimized with this pipeline" << std::endl;
        return;
    }
    int32_t minScore = open_gap_penalty;
    if( (maxLength * minScore) < -pow(2,31)){
        std::cerr << "the sequence is too long, could not be optimized with this pipeline" << std::endl;
        return;
    } //check potential overflow problem end
    __m256i* vProfile = this->query_profile_avx2_byte(_nucleotideCodeSubstitutionMatrix);
    this->ssw_avx2(vProfile, _nucleotideCodeSubstitutionMatrix);
    this->get_optimal_alignment();
}

__m256i* NeedlemanWunsch_simd_Slow::query_profile_avx2_byte(  NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix){ // dna_q is the query, (the non-reference sequence, with gene structure unknown)
    __m256i* vProfile = (__m256i*)malloc(this->_n * this->_segLen * sizeof(__m256i));
    int32_t a, i, j, k;
    for( a =0; a<_n; ++a ){
        //vProfile[a] = new __m256i[segLen];
        int32_t* s = (int32_t*)(vProfile+a*this->_segLen);
        for( i=0; i<this->_segLen; ++i ){
            j = i;
            for( k=0; k<8; ++k ){
                *s++ = j>= _length_of_q ? bias : this->_substitute_matrix[a][nucleotideCodeSubstitutionMatrix.get_dna_acid_map(_dna_q[j])];
                j += this->_segLen;
            }
        }
    }
    return vProfile;
}

void NeedlemanWunsch_simd_Slow::ssw_avx2(__m256i* vProfile, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix ){
    const __m256i vZero = _mm256_set1_epi32(0);
    const __m256i vMinuteOne = _mm256_set1_epi32(-1);
    const __m256i vOne = _mm256_set1_epi32(1);

    __m256i* pvHStore = (__m256i*) calloc(this->_segLen, sizeof(__m256i));
    __m256i* pvHLoad = (__m256i*) calloc(this->_segLen, sizeof(__m256i));
    __m256i* trackMatch = (__m256i*) calloc(this->_segLen, sizeof(__m256i));
    __m256i* trackDel = (__m256i*) calloc(this->_segLen, sizeof(__m256i));
    __m256i* pv;
    __m256i* vP;
    __m256i vcmp, vGapO= _mm256_set1_epi32(_open_gap_penalty), vGapE=_mm256_set1_epi32(_extend_gap_penalty);

    int32_t i, j, z, t;
    //for outer loop
    __m256i vE, vH, vHsave;
    int32_t deletionScore;
    int32_t* a;
    int32_t* b;
    int32_t* c;
    int lastOne;
    for (i = 0; i<_length_of_d; ++i) {//important

        vP = vProfile + ((nucleotideCodeSubstitutionMatrix.get_dna_acid_map(_dna_d[i] ))*(this->_segLen));
        vH = _mm256_loadu_si256(pvHStore+(this->_segLen-1));// the last one
        vH = _mm256_alignr_epi8(vH, _mm256_permute2x128_si256(vH, vH, _MM_SHUFFLE(0, 0, 2, 0)), 16 - 4);

        pv = pvHLoad;
        pvHLoad = pvHStore;
        pvHStore = pv;
        // inner loop to process the query sequence
        for (j = 0; j < this->_segLen; ++j) {
            vcmp = _mm256_loadu_si256(trackDel + j);
            if( 0 == i ){
                vE = _mm256_add_epi32(vGapE, _mm256_loadu_si256(pvHLoad + j));
            } else {
                __m256i t1 = _mm256_add_epi32(_mm256_mullo_epi32 (_mm256_add_epi32(vcmp, vOne), vGapO),
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
                deletionScore = _similarity_matrix[j][i+1] + (int32_t)(_extend_gap_penalty);
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
                deletionScore = _similarity_matrix[j][i+1]+(int32_t)(_open_gap_penalty);
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
    free(vProfile);
}


// Trace back step.
void NeedlemanWunsch_simd_Slow::get_optimal_alignment( ) {
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
            if ( _track_match[i][j]==-1  ) {
                SQ.push(_dna_q[i - 1]);
                SD.push(_dna_d[j - 1]);
                i -= 1;
                j -= 1;
            } else if (_track_ins[i][j] == -1 ) {//insertion
                SQ.push(_dna_q[i - 1]);
                SD.push('-');
                i -= 1;
            } else {        // Going to S(i, j-1) //deletion
                SQ.push('-');
                SD.push(_dna_d[j - 1]);
                j -= 1;
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

void NeedlemanWunsch_simd_Slow::print_results() {
    std::cout << _alignment_q << std::endl;
    std::cout << _alignment_d << std::endl;

}

NeedlemanWunsch_simd_Slow::~NeedlemanWunsch_simd_Slow(){
    int i;
    for ( i = 0; i < (_length_of_q + 1); i++) {
        delete[] this->_similarity_matrix[i];
        delete[] this->_track_match[i];
        delete[] this->_track_del[i];
        delete[] this->_track_ins[i];
    }
    delete[] this->_similarity_matrix;
    delete[] this->_track_match;
    delete[] this->_track_del;
    delete[] this->_track_ins;
    for( i = 0; i < this->_n; i++ ){
        delete [] this->_substitute_matrix[i];
    }
    delete [] this->_substitute_matrix;

}

int32_t ** NeedlemanWunsch_simd_Slow::getSimilarity_matrix(){
    return this->_similarity_matrix;
}

std::string NeedlemanWunsch_simd_Slow::getAlignment_q(){
    return _alignment_q;
}
std::string NeedlemanWunsch_simd_Slow::getAlignment_d(){
    return _alignment_d;
}




























NeedlemanWunsch_simd_Fast::NeedlemanWunsch_simd_Fast(const std::string & dna_q, const std::string & dna_d, const int8_t & match_score,
                                           const int8_t & mis_match_score, const int8_t & open_gap_penalty, const int8_t & extend_gap_penalty, const NucleotideCodeSubstitutionMatrix & _nucleotideCodeSubstitutionMatrix) {

    this->_dna_q = dna_q;
    this->_dna_d = dna_d;
    this->_open_gap_penalty = open_gap_penalty*(int8_t)36;
    this->_extend_gap_penalty = extend_gap_penalty*(int8_t)36;
    this->_length_of_q = this->_dna_q.length();
    this->_length_of_d = this->_dna_d.length();
    this->bias = 0;
    this->_n = 5;

    this->_segLen = ( _length_of_q+ 7)/8;
    //score matrix and track matrix begin

    this->_similarity_matrix = new int32_t*[this->_length_of_q+1];
    int i, j;
    for (i = 0; i < (this->_length_of_q+1); ++i) {
        this->_similarity_matrix[i] = new int32_t[this->_length_of_d+1];
    }
    i=0;
    for ( j = 0; j <= _length_of_d; ++j) {
        this->_similarity_matrix[i][j] = 0;
    }
    j=0;
    for (i = 0; i < (this->_length_of_q+1); ++i) {
        this->_similarity_matrix[i][j] = 0;
    }
    this->_substitute_matrix = new int32_t*[_n];
    for( i=0; i<_n; ++i){
        this->_substitute_matrix[i] = new int32_t[_n];
        for( j=0; j<_n; ++j){
            if( _nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[i][j] > 0 ){
                this->_substitute_matrix[i][j] = _nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[i][j] * match_score;
            }else{
                this->_substitute_matrix[i][j] = -(_nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[i][j] * mis_match_score);
            }
        }
    }

    // here let's assume we have sequence with length 10M, the score of exon match is 36, then the total score could be 3.6e+08
    // the value range of int32 is about -2e+09 -2e+09, so it should be enough
    //check potential overflow problem begin
    uint64_t maxLength = _length_of_d;
    if( maxLength <  _length_of_q){
        maxLength = _length_of_q;
    }
    int32_t maxScore = 36 * match_score;
    if( (maxLength * maxScore) > pow(2,31)){
        std::cerr << "the sequence is too long, could not be optimized with this pipeline" << std::endl;
        return;
    }
    int32_t minScore = -36 * open_gap_penalty;
    if( (maxLength * minScore) < -pow(2,31)){
        std::cerr << "the sequence is too long, could not be optimized with this pipeline" << std::endl;
        return;
    } //check potential overflow problem end

}

__m256i* NeedlemanWunsch_simd_Fast::query_profile_avx2_byte( NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix){ // dna_q is the query, (the non-reference sequence, with gene structure unknown)
    __m256i* vProfile = (__m256i*)malloc(this->_n * this->_segLen * sizeof(__m256i));
    int32_t a, i, j, k;
    for( a =0; a<_n; ++a ){
        //vProfile[a] = new __m256i[segLen];
        int32_t* s = (int32_t*)(vProfile+a*this->_segLen);
        for( i=0; i<this->_segLen; ++i ){
            j = i;
            for( k=0; k<8; ++k ){
                *s++ = j>= _length_of_q ? bias : this->_substitute_matrix[a][nucleotideCodeSubstitutionMatrix.get_dna_acid_map(_dna_q[j])];
                j += this->_segLen;
            }
        }
    }
    return vProfile;
}

void NeedlemanWunsch_simd_Fast::ssw_avx2(__m256i* vProfile, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix ){
    //* Define 0 vector.
    const __m256i vZero = _mm256_set1_epi32(0);

    __m256i* pvHStore = (__m256i*) calloc(this->_segLen, sizeof(__m256i));
    __m256i* pvHLoad = (__m256i*) calloc(this->_segLen, sizeof(__m256i));
    __m256i* pvE = (__m256i*) calloc(this->_segLen, sizeof(__m256i)); //vE is the score from left side
    __m256i* pv;
    __m256i* vP;

    int32_t i, j, z;

    // INDEL begin vector
    const __m256i vGapO = _mm256_set1_epi32(this->_open_gap_penalty); // here I use epi32, because for the score matrix I used int32_t
    // INDEL extension vector
    const __m256i vGapE = _mm256_set1_epi32(this->_extend_gap_penalty);

    for(i=0; i<_segLen; ++i){
        _mm256_storeu_si256(pvE+i, vGapE);
    }
    //for outer loop
    __m256i vE, vH;
    double deletionScore;
    const int32_t* a;

    for (i = 0; (i != _length_of_d); ++i) {//important
        vH = _mm256_loadu_si256(pvHStore+(this->_segLen-1));// the last one

        vH = _mm256_alignr_epi8(vH, _mm256_permute2x128_si256(vH, vH, _MM_SHUFFLE(0, 0, 2, 0)), 16 - 4);
        // Shift the value in vH left by 4 byte. one score matrix cell
        // _mm256_slli_si256 not not work here, because it shifts 8-bit [byte] elements within a 128-bit lane of the source vector

        vP = vProfile + ((nucleotideCodeSubstitutionMatrix.get_dna_acid_map(_dna_d[i] ))*(this->_segLen)); // Correct part of the vProfile

        pv = pvHLoad;
        pvHLoad = pvHStore;
        pvHStore = pv;

        // inner loop to process the query sequence
        for (j = 0; j < this->_segLen; ++j) {
            vH = _mm256_add_epi32(vH, _mm256_loadu_si256(vP+j)); // add the scoring profile to vH
            // Get max from vH, vE
            vE = _mm256_loadu_si256(pvE+j);
            vH = _mm256_max_epi32(vH, vE);

            // Save vH values.
            _mm256_storeu_si256(pvHStore + j, vH);

            // Update vE value. for next database base pair  here should be get the exact score as general smith-waterman algorithm
            vH = _mm256_add_epi32(vH, vGapO);
            vE = _mm256_add_epi32(vE, vGapE);
            vE = _mm256_max_epi32(vE, vH);
            _mm256_storeu_si256(pvE+j, vE);

            // Load the next vH.
            vH = _mm256_loadu_si256(pvHLoad + j);
        }

        // put the score into the score matrix begin
        for( j=0; j<this->_segLen; ++j ){
            a = ( int32_t*)(pvHStore+j);

            for( z=0; z<8; ++z ){
                if( (j + z*this->_segLen) < _length_of_q ){
                    _similarity_matrix[j+z*this->_segLen+1][i+1] = a[z];
                }
            }
        }

        int lastOne = 1; // INDEL
        for( j=0; j<_length_of_q; ++j ){
            if( lastOne == 1 ){
                deletionScore = _similarity_matrix[j][i+1]+_extend_gap_penalty;
                if( deletionScore > _similarity_matrix[j+1][i+1]  ){
                    _similarity_matrix[j+1][i+1] = deletionScore;
                }else{
                    lastOne=0;//MATCH OR MISMATCH
                }
            }else{
                deletionScore = _similarity_matrix[j][i+1]+_open_gap_penalty;
                if( deletionScore > _similarity_matrix[j+1][i+1]  ){
                    _similarity_matrix[j+1][i+1] = deletionScore;
                }else{
                    lastOne=0;//MATCH OR MISMATCH
                }
            }
        }
        // put the score into the score matrix end

        //updata vH begin
        for( j=0; j<this->_segLen; j++ ){
            auto* a = (int32_t*)(pvHStore+j);
            for( z=0; z<8; ++z ){
                if( (j + z*this->_segLen) < _length_of_q ){
                    a[z] = _similarity_matrix[j+z*this->_segLen+1][i+1];
                }
            }
        }//updata vH end
    }
    free(pvE);
    free(pvHStore);
    free(vProfile);
    free(pvHLoad);
}



void NeedlemanWunsch_simd_Fast::get_optimal_alignment( NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix ){

    _alignment_q = "";
    _alignment_d = "";

    std::stack<char> SA, SB;
    size_t k;
    size_t i = _length_of_q, j = _length_of_d;
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
        SA.push(_dna_q[k-1]);
        SB.push('-');
    }

    for( k=_length_of_d; k>j; --k ){
        SA.push('-');
        SB.push(_dna_d[k-1]);
    }
    while (i != 0 || j != 0) {
        if (i == 0) {
            SA.push('-');
            SB.push(_dna_d[j-1]);
            --j;
        } else if (j == 0) {
            SA.push(_dna_q[i-1]);
            SB.push('-');
            --i;
        } else {
            int index_i = nucleotideCodeSubstitutionMatrix.get_dna_acid_map(_dna_q[i-1]);
            int index_j = nucleotideCodeSubstitutionMatrix.get_dna_acid_map(_dna_d[j-1]);
            int S = this->_substitute_matrix[index_i][index_j];//  (_dna_q[ii-1] == _dna_d[jj-1]) ? match_score : -mismatch_score;
            if ( this->_similarity_matrix[i][j] == (this->_similarity_matrix[i-1][j-1] + S) ) {
                SA.push(_dna_q[i-1]);
                SB.push(_dna_d[j-1]);
                --i; --j;
            } else if (_similarity_matrix[i-1][j] > _similarity_matrix[i][j-1]) {
                SA.push(_dna_q[i-1]);
                SB.push('-');
                --i;
            } else {
                SA.push('-');
                SB.push(_dna_d[j-1]);
                --j;
            }
        }
    }
    while (!SA.empty()) {
        _alignment_q += SA.top();
        _alignment_d += SB.top();
        SA.pop();
        SB.pop();
    }
}

void NeedlemanWunsch_simd_Fast::print_results() {

    int32_t _similarity=0, _identity=0;
    size_t i;
    for ( i = 0; i < _alignment_q.length(); i++) {
        std::cout << _alignment_q[i];
    }
    std::cout << "\n";


    for (i = 0; i < _alignment_d.length(); i++){
        std::cout << _alignment_d[i];
    }
    std::cout << "\n";
    std::cout << std::setfill(' ') << std::setw(50) << "\n";

    double percentage;
    std::cout << "Score: " << _similarity_matrix[_length_of_q][_length_of_d] << "\n";

    std::cout << "Length: " << _alignment_q.length() << " (with gaps)\n";

    percentage = ((double)_identity / (double)_alignment_q.length()) * 100;
    std::cout << "Identity: " << _identity << "/" <<_alignment_q.length() << " ( %" << percentage << " ) " << "\n";

    percentage = ((double)_similarity / (double)_alignment_q.length()) * 100;
    std::cout << "Similarity: " << _similarity << "/" << _alignment_q.length() <<" ( %" << percentage << " ) " << "\n";

}

NeedlemanWunsch_simd_Fast::~NeedlemanWunsch_simd_Fast(){
    int32_t i;
    for ( i = 0; i < (_length_of_q+1); i++ ) {
        delete [] this->_similarity_matrix[i];
    }
    delete [] this->_similarity_matrix;
    for( i = 0; i < this->_n; i++ ){
        delete [] this->_substitute_matrix[i];
    }
    delete [] this->_substitute_matrix;

}

int32_t ** NeedlemanWunsch_simd_Fast::getSimilarity_matrix(){
    return this->_similarity_matrix;
}

std::string NeedlemanWunsch_simd_Fast::get_alignment_q(){
    return this->_alignment_q;
}
std::string NeedlemanWunsch_simd_Fast::get_alignment_d(){
    return this->_alignment_d;
}
