//
// Created by Baoxing song on 2018-12-26.
//

/*
 * I ever read the source code of minimap2 to write the source code of this file
 * Some code fragments is borrowed from minimap2
 * */

#include "hengAlign.h"

#define KSW_NEG_INF -0x40000000

#define kcalloc(km, count, size) calloc((count), (size))
#define krealloc(km, ptr, size) realloc((ptr), (size))
#define kmalloc(km, size) malloc((size))
#define kfree(km, ptr) free((ptr))

typedef struct {
    //uint32_t max:31, zdropped:1;
    //int max_q, max_t;      // max extension coordinate
    int mqe, mqe_t;        // max score when reaching the end of query
    int mte, mte_q;        // max score when reaching the end of target
    //int score;             // max score reaching both ends; may be KSW_NEG_INF
    int m_cigar, n_cigar;
    //int reach_end;
    uint32_t *cigar;
} ksw_extz_t;


static inline void ksw_reset_extz(ksw_extz_t *ez){
    /*ez->max_q = ez->max_t = */ez->mqe_t = ez->mte_q = -1;
    /*ez->max = 0, ez->score = */ez->mqe = ez->mte = KSW_NEG_INF;
    ez->n_cigar = 0;//, /*ez->zdropped = 0, ez->reach_end = 0*/;
}

typedef struct header_t {
    size_t size;
    struct header_t *ptr;
} header_t;


struct mm_tbuf_s {
    void *km;
    int rep_len, frag_gap;
};
typedef struct {
    header_t base, *loop_head, *core_head; /* base is a zero-sized block always kept in the loop */
} kmem_t;

void *km_init(void)
{
    return calloc(1, sizeof(kmem_t));
}
int mm_dbg_flag = 0;
// memory buffer for thread-local storage during mapping
typedef struct mm_tbuf_s mm_tbuf_t;
mm_tbuf_t *mm_tbuf_init(void){
    mm_tbuf_t *b;
    b = (mm_tbuf_t*)calloc(1, sizeof(mm_tbuf_t));
    if (!(mm_dbg_flag & 1)) b->km = km_init();
    return b;
}

static inline uint32_t *ksw_push_cigar(void *km, int *n_cigar, int *m_cigar, uint32_t *cigar, uint32_t op, int len){
    if (*n_cigar == 0 || op != (cigar[(*n_cigar) - 1]&0xf)) {
        if (*n_cigar == *m_cigar) {
            *m_cigar = *m_cigar? (*m_cigar)<<1 : 4;
            cigar = (uint32_t*)krealloc(km, cigar, (*m_cigar) << 2);
        }
        cigar[(*n_cigar)++] = len<<4 | op;
    } else cigar[(*n_cigar)-1] += len<<4;
    return cigar;
}

static inline void ksw_backtrack(void *km, const uint8_t *p, const int *off, const int *off_end, int n_col, int i0, int j0,
                                 int *m_cigar_, int *n_cigar_, uint32_t **cigar_){ // p[] - lower 3 bits: which type gets the max; bit
    int n_cigar = 0, m_cigar = *m_cigar_, i = i0, j = j0, r, state = 0;
    uint32_t *cigar = *cigar_, tmp;
    while (i >= 0 && j >= 0) { // at the beginning of the loop, _state_ tells us which state to check
        int force_state = -1;
        r = i + j;
        if (i < off[r]) force_state = 2;
        if (off_end && i > off_end[r]) force_state = 1;
        tmp = force_state < 0? p[(size_t)r * n_col + i - off[r]] : 0;
//        std::cout << "line 112 tmp " << tmp << " force_state " << force_state<< std::endl;
        if (state == 0) state = tmp & 7; // if requesting the H state, find state one maximizes it.
        else if (!(tmp >> (state + 2) & 1)) state = 0; // if requesting other states, _state_ stays the same if it is a continuation; otherwise, set to H
        if (state == 0) state = tmp & 7; // TODO: probably this line can be merged into the "else if" line right above; not 100% sure
        if (force_state >= 0) state = force_state;
        if (state == 0) cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, 0, 1), --i, --j; // match
        else if (state == 1 || (state == 3 )) cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, 2, 1), --i; // deletion
        else cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, 1, 1), --j; // insertion
    }
    if (i >= 0) cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, 2, i + 1); // first deletion
    if (j >= 0) cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, 1, j + 1); // first insertion
//    for (i = 0; i < n_cigar>>1; ++i) {// reverse CIGAR
//        tmp = cigar[i], cigar[i] = cigar[n_cigar - 1 - i], cigar[n_cigar - 1 - i] = tmp;
//        std::cout << "reverse" << std::endl;
//    }
    *m_cigar_ = m_cigar, *n_cigar_ = n_cigar, *cigar_ = cigar;
}


void hengAlign(const std::string& _dna_q, const std::string& _dna_d, std::stack<char> & SQ, std::stack<char> & SD,
                         const int8_t & _open_gap_penalty, const int8_t & _extend_gap_penalty,
                         NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix) {
    mm_tbuf_t *tbuf = mm_tbuf_init();
    void *km = tbuf->km;
    int qlen = _dna_q.length();
    int dlen = _dna_d.length();
    int8_t q = -(_open_gap_penalty - _extend_gap_penalty);
    int8_t e = - _extend_gap_penalty;

//    r=i+j
//    t=i

#define __dp_code_block1 \
	z = _mm_add_epi8(_mm_load_si128(&s[t]), qe2_);  \
	xt1 = _mm_load_si128(&x[t]);                     /* xt1 <- x[r-1][t..t+15] */ \
	tmp = _mm_srli_si128(xt1, 15);                   /* tmp <- x[r-1][t+15] */ \
	xt1 = _mm_or_si128(_mm_slli_si128(xt1, 1), x1_); /* xt1 <- x[r-1][t-1..t+14] */ \
	x1_ = tmp; \
	vt1 = _mm_load_si128(&v[t]);                     /* vt1 <- v[r-1][t..t+15] */ \
	tmp = _mm_srli_si128(vt1, 15);                   /* tmp <- v[r-1][t+15] */ \
	vt1 = _mm_or_si128(_mm_slli_si128(vt1, 1), v1_); /* vt1 <- v[r-1][t-1..t+14] */ \
	v1_ = tmp; \
	a = _mm_add_epi8(xt1, vt1);                      /* a <- x[r-1][t-1..t+14] + v[r-1][t-1..t+14] */ \
	ut = _mm_load_si128(&u[t]);                      /* ut <- u[t..t+15] */ \
	b = _mm_add_epi8(_mm_load_si128(&y[t]), ut);     /* b <- y[r-1][t..t+15] + u[r-1][t..t+15] */


#define __dp_code_block2 \
	z = _mm_max_epu8(z, b);                          /* z = max(z, b); this works because both are non-negative */ \
	z = _mm_min_epu8(z, max_sc_); \
	_mm_store_si128(&u[t], _mm_sub_epi8(z, vt1));    /* u[r][t..t+15] <- z - v[r-1][t-1..t+14] */ \
	_mm_store_si128(&v[t], _mm_sub_epi8(z, ut));     /* v[r][t..t+15] <- z - u[r-1][t..t+15] */ \
	z = _mm_sub_epi8(z, q_); \
	a = _mm_sub_epi8(a, z); \
	b = _mm_sub_epi8(b, z);

    int r, t, qe = q + e, n_col_, *off = 0, *off_end = 0, tlen_, qlen_, last_st, last_en/*, wl, wr, max_sc, min_sc*/;
    int32_t *H = 0;//, H0 = 0, last_H0_t = 0;
    uint8_t *qr, *sf, *mem, *mem2 = 0;
    __m128i q_, qe2_, zero_, flag1_, flag2_, flag8_, flag16_, sc_mch_, sc_mis_, sc_N_, m1_, max_sc_;
    __m128i *u, *v, *x, *y, *s, *p = 0;

//    std::cout << "line 173" << std::endl;
    ksw_extz_t ez1;
    memset(&ez1, 0, sizeof(ksw_extz_t));
    ksw_extz_t * ez = &ez1;
//    std::cout << "line 177" << std::endl;
    ksw_reset_extz(ez);
//    std::cout << "line 179" << std::endl;
    int m = 5;
    zero_   = _mm_set1_epi8(0);
    q_      = _mm_set1_epi8(q); //open gap penalty
    qe2_    = _mm_set1_epi8((q + e) * 2);
    flag1_  = _mm_set1_epi8(1);
    flag2_  = _mm_set1_epi8(2);
    flag8_  = _mm_set1_epi8(0x08);
    flag16_ = _mm_set1_epi8(0x10);
    sc_mch_ = _mm_set1_epi8(nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[0][0]);
    sc_mis_ = _mm_set1_epi8(nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[0][1]);
    sc_N_   = nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[m-1][m-1]==0 ?
              _mm_set1_epi8(-_extend_gap_penalty) : _mm_set1_epi8(nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[m-1][m-1]);
    m1_     = _mm_set1_epi8(m - 1); // wildcard it is N of the sequence
    max_sc_ = _mm_set1_epi8(nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[0][0] + (q + e) * 2);

    int w = dlen > qlen? dlen : qlen;

    //max_sc = nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[0][0];
    //min_sc = nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[0][1];

    tlen_ = (dlen + 15) / 16;
    n_col_ = qlen < dlen? qlen : dlen;
    n_col_ = ((n_col_ < w + 1? n_col_ : w + 1) + 15) / 16 + 1;
    qlen_ = (qlen + 15) / 16;

    mem = (uint8_t*)kcalloc(km, tlen_ * 6 + qlen_ + 1, 16);
    u = (__m128i*)(((size_t)mem + 15) >> 4 << 4); // 16-byte aligned
    v = u + tlen_;
    x = v + tlen_;
    y = x + tlen_;
    s = y + tlen_;
    sf = (uint8_t*)(s + tlen_);
    qr = sf + tlen_ * 16;
    // the baove code is for memory location and put the pointer at the correct position
    //std::cout << "tlen_ " << tlen_ << std::endl;

    H = (int32_t*)kmalloc(km, tlen_ * 16 * 4);
    for (t = 0; t < tlen_ * 16; ++t){
        H[t] = KSW_NEG_INF;
    }

    mem2 = (uint8_t*)kmalloc(km, ((size_t)(qlen + dlen - 1) * n_col_ + 1) * 16);
    p = (__m128i*)(((size_t)mem2 + 15) >> 4 << 4);
    off = (int*)kmalloc(km, (qlen + dlen - 1) * sizeof(int) * 2);
    off_end = off + qlen + dlen - 1;

    for (t = 0; t < qlen; ++t) qr[t] = nucleotideCodeSubstitutionMatrix.get_dna_acid_map(_dna_q[qlen - 1 - t]);
    //std::cout << "line 100" << std::endl;
    for (t = 0; t < dlen; ++t) sf[t] = nucleotideCodeSubstitutionMatrix.get_dna_acid_map(_dna_d[t]);
    //std::cout << "line 102" << std::endl;
    for (r = 0, last_st = last_en = -1; r < qlen + dlen - 1; ++r) { //r=i+j  t=i
        //std::cout << "line 104" << std::endl;
        int st = 0, en = dlen - 1, st0, en0, st_, en_;
        int8_t x1, v1;

        uint8_t *qrr = qr + (qlen - 1 - r);
        uint8_t *u8 = (uint8_t*)u, *v8 = (uint8_t*)v;
        __m128i x1_, v1_;
        // find the boundaries
        //std::cout << "line 245 st " << st << std::endl;
        if (st < r - qlen + 1) st = r - qlen + 1;
        //std::cout << "line 247 st " << st << std::endl;
        if (en > r) en = r;
        /*
        if (st < (r-w+1)>>1) st = (r-w+1)>>1; // take the ceil
        //std::cout << "line 250 st " << st << std::endl;
        if (en > (r+w)>>1) en = (r+w)>>1; // take the floor
*/
        st0 = st, en0 = en;
        st = st / 16 * 16, en = (en + 16) / 16 * 16 - 1;

        // set boundary conditions
        if (st > 0) {
            if (st - 1 >= last_st && st - 1 <= last_en)
                x1 = ((uint8_t*)x)[st - 1], v1 = v8[st - 1]; // (r-1,s-1) calculated in the last round
            else x1 = v1 = 0; // not calculated; set to zeros
        } else x1 = 0, v1 = r? q : 0;

        if (en >= r) ((uint8_t*)y)[r] = 0, u8[r] = r? q : 0;

        // loop fission: set scores first
        //std::cout << "r " << r << " st0 " << st0 << " en0 " << en0 << std::endl;
        for (t = st0; t <= en0; t += 16) { // this should be same with the another(else) part, but faster
            // all the snp adding scores the S matrix in the function
            __m128i sq, st, tmp, mask;
            sq = _mm_loadu_si128((__m128i*)&sf[t]);
            st = _mm_loadu_si128((__m128i*)&qrr[t]);
            mask = _mm_or_si128(_mm_cmpeq_epi8(sq, m1_), _mm_cmpeq_epi8(st, m1_)); //m1_ is N
            tmp = _mm_cmpeq_epi8(sq, st);
#ifdef __SSE4_1__
            tmp = _mm_blendv_epi8(sc_mis_, sc_mch_, tmp);
            tmp = _mm_blendv_epi8(tmp,     sc_N_,   mask);
#else
            tmp = _mm_or_si128(_mm_andnot_si128(tmp,  sc_mis_), _mm_and_si128(tmp,  sc_mch_));
				tmp = _mm_or_si128(_mm_andnot_si128(mask, tmp),     _mm_and_si128(mask, sc_N_));
#endif
            _mm_storeu_si128((__m128i*)((uint8_t*)s + t), tmp);
        }

        std::cout << "r " << r << " st0 " << st0 << " en0 " << en0 << std::endl;

        // core loop
        x1_ = _mm_cvtsi32_si128(x1);  //Copy 32-bit integer a to the lower elements of dst, and zero the upper elements of dst.
        v1_ = _mm_cvtsi32_si128(v1);
        st_ = st / 16, en_ = en / 16;
        assert(en_ - st_ + 1 <= n_col_);
//        std::cout << "line 240" << std::endl;
        __m128i *pr = p + (size_t)r * n_col_ - st_; //todo what pr is for baoxing song this is only important for cigar
        off[r] = st, off_end[r] = en;
//        std::cout << "line 243" << std::endl;
        for (t = st_; t <= en_; ++t) {
            __m128i d, z, a, b, xt1, vt1, ut, tmp;
            __dp_code_block1;
            d = _mm_and_si128(_mm_cmpgt_epi8(a, z), flag1_); // d = a > z? 1 : 0
#ifdef __SSE4_1__
            z = _mm_max_epi8(z, a);                          // z = z > a? z : a (signed)
            tmp = _mm_cmpgt_epi8(b, z);
            d = _mm_blendv_epi8(d, flag2_, tmp);             // d = b > z? 2 : d
#else // we need to emulate SSE4.1 intrinsics _mm_max_epi8() and _mm_blendv_epi8()
            z = _mm_and_si128(z, _mm_cmpgt_epi8(z, zero_));  // z = z > 0? z : 0;
				z = _mm_max_epu8(z, a);                          // z = max(z, a); this works because both are non-negative
				tmp = _mm_cmpgt_epi8(b, z);
				d = _mm_or_si128(_mm_andnot_si128(tmp, d), _mm_and_si128(tmp, flag2_)); // d = b > z? 2 : d; emulating blendv
#endif
            __dp_code_block2;
            tmp = _mm_cmpgt_epi8(a, zero_);
            _mm_store_si128(&x[t], _mm_and_si128(tmp, a));
            d = _mm_or_si128(d, _mm_and_si128(tmp, flag8_));  // d = a > 0? 0x08 : 0
            tmp = _mm_cmpgt_epi8(b, zero_);
            _mm_store_si128(&y[t], _mm_and_si128(tmp, b));
            d = _mm_or_si128(d, _mm_and_si128(tmp, flag16_)); // d = b > 0? 0x10 : 0
            _mm_store_si128(&pr[t], d);
        }
//        std::cout << "line 279" << std::endl;
        {
            int32_t max_H, max_t;
            // compute H[], max_H and max_t
            if (r > 0) {
                int32_t HH[4], tt[4], en1 = st0 + (en0 - st0) / 4 * 4, i;
                __m128i max_H_, max_t_, qe_;
                max_H = H[en0] = e;//en0 > 0? H[en0-1] + u8[en0] - qe : H[en0] + v8[en0] - qe; // special casing the last element
                max_t = en0;
                max_H_ = _mm_set1_epi32(max_H);
                max_t_ = _mm_set1_epi32(max_t);
                qe_    = _mm_set1_epi32(q + e);
                for (t = st0; t < en1; t += 4) { // this implements: H[t]+=v8[t]-qe; if(H[t]>max_H) max_H=H[t],max_t=t;
                    __m128i H1, tmp, t_;
                    H1 = _mm_loadu_si128((__m128i*)&H[t]);
                    t_ = _mm_setr_epi32(v8[t], v8[t+1], v8[t+2], v8[t+3]);
                    H1 = _mm_add_epi32(H1, t_);
                    H1 = _mm_sub_epi32(H1, qe_);

                    _mm_storeu_si128((__m128i*)&H[t], H1);
                    t_ = _mm_set1_epi32(t);
                    tmp = _mm_cmpgt_epi32(H1, max_H_);
#ifdef __SSE4_1__
                    max_H_ = _mm_blendv_epi8(max_H_, H1, tmp);
                    max_t_ = _mm_blendv_epi8(max_t_, t_, tmp);
#else
                    max_H_ = _mm_or_si128(_mm_and_si128(tmp, H1), _mm_andnot_si128(tmp, max_H_));
					max_t_ = _mm_or_si128(_mm_and_si128(tmp, t_), _mm_andnot_si128(tmp, max_t_));
#endif
                }
                _mm_storeu_si128((__m128i*)HH, max_H_);
                _mm_storeu_si128((__m128i*)tt, max_t_);

                for (i = 0; i < 4; ++i)
                    if (max_H < HH[i]) max_H = HH[i], max_t = tt[i] + i, std::cout << "i " << tt[i] + i << " j " << r-tt[i]-i << " max_H " << max_H << std::endl;
                for (; t < en0; ++t) { // for the rest of values that haven't been computed with SSE
                    H[t] += (int32_t)v8[t] - qe;
                    if (H[t] > max_H)
                        max_H = H[t], max_t = t;
                }
            } else H[0] = v8[0] - qe - qe, max_H = H[0], max_t = 0; // special casing r==0
            //std::cout << "r " << r << " max_H " << max_H << std::endl;
            // update ez
            if (en0 == dlen - 1 && H[en0] > ez->mte){
                ez->mte = H[en0], ez->mte_q = r - en0;
                std::cout << "line 333 i " << dlen - 1 << " j " << ez->mte_q << " highestScore " << ez->mte << std::endl;
            }
            if (r - st0 == qlen - 1 && H[st0] > ez->mqe){
                ez->mqe = H[st0], ez->mqe_t = st0;
                std::cout << "line 337 i " << ez->mqe_t << " j " <<qlen - 1 << " highestScore " << ez->mqe << std::endl;
            }
/*
            if (max_H > (int32_t)ez->max)
                ez->max = max_H, ez->max_t = max_t, ez->max_q = r - max_t;

            if (r == qlen + dlen - 2 && en0 == dlen - 1)
                ez->score = H[dlen - 1];
                */
        }
        last_st = st, last_en = en;
    }
    kfree(km, mem);
    kfree(km, H);

    { // backtrack
        int _dna_q_i=0;
        int _dna_d_i=0;
        if( ez->mte > ez->mqe ){
            std::cout << "line 352" << std::endl;
            ksw_backtrack(km, (uint8_t*)p, off, off_end, n_col_*16, dlen-1, ez->mte_q, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
            _dna_q_i=ez->mte_q;
            _dna_d_i=dlen-1;
        }else{
            std::cout << "line 357" << std::endl;
            ksw_backtrack(km, (uint8_t*)p, off, off_end, n_col_*16, ez->mqe_t, qlen-1, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
            _dna_q_i=qlen-1;
            _dna_d_i=ez->mqe_t;
        }
        std::cout << "_dna_q_i " << _dna_q_i << " _dna_d_i " << _dna_d_i << std::endl;
        // due to the order reason, it is difficult to put the result into vector directly
        for (int i = 0; i<ez->n_cigar; ++i){ // IMPORTANT: this gives the CIGAR in the aligned regions. NO soft/hard clippings!
            // the cigar is in reverse order, and that is what I want
            int cigarLength = ez->cigar[i]>>4;
            int cigarType = ez->cigar[i]&0xf;
            printf("%d%c", cigarLength, "MIDNSH"[cigarType]);
            if( cigarType == 0 ){ //M
                for ( int j=0; j <cigarLength; ++j ){
                    SQ.push(_dna_q[_dna_q_i]);
                    SD.push(_dna_d[_dna_d_i]);
                    --_dna_q_i;
                    --_dna_d_i;
                }
            }else if( cigarType == 1 ){ //insertion
                for ( int j=0; j <cigarLength; ++j ){
                    SQ.push(_dna_q[_dna_q_i]);
                    SD.push('-');
                    --_dna_q_i;
                }
            }else if( cigarType == 2 ){ //deletion
                for ( int j=0; j <cigarLength; ++j ){
                    SQ.push('-');
                    SD.push(_dna_d[_dna_d_i]);
                    --_dna_d_i;
                }
            }else{
                printf(" we could not work with this cigar %d%c", cigarLength, "MIDNSH"[cigarType]);
            }
        }
        std::cout << std::endl;
        kfree(km, mem2); kfree(km, off);
        kfree(km, ez1.cigar);
    }
}
