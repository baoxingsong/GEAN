//
// Created by Baoxing song on 2018-12-26.
//

/*
 * I ever read the source code of minimap2 to write the source code of this file
 * Some code fragments is borrowed from minimap2
 * https://github.com/lh3/minimap2/blob/master/ksw2_exts2_sse.c
 * */

#include "hengAlign.h"

void *km_init(void){
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

void km_destroy(void *_km){
    kmem_t *km = (kmem_t*)_km;
    header_t *p, *q;
    if (km == NULL) return;
    for (p = km->core_head; p != NULL;) {
        q = p->ptr;
        free(p);
        p = q;
    }
    free(km);
}

void mm_tbuf_destroy(mm_tbuf_t *b){
    if (b == 0) return;
    km_destroy(b->km);
    free(b);
}

static inline uint32_t *ksw_push_cigar(void *km, int *n_cigar, int *m_cigar, uint32_t *cigar, uint32_t op, int len){
    if (*n_cigar == 0 || op != (cigar[(*n_cigar) - 1]&0xf)) { // this is the first cigar or the cigar_op is different with last one
        if (*n_cigar == *m_cigar) { // it is the time to double the memory)
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
                  NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix){
    int qlen = _dna_q.length();
    int tlen = _dna_d.length();
    uint8_t *query = new uint8_t[qlen];
    uint8_t *target = new uint8_t[tlen];
    int ts;
    for (ts = 0; ts < qlen; ++ts) query[ts] = nucleotideCodeSubstitutionMatrix.get_dna_acid_map(_dna_q[ts]);
    for (ts = 0; ts < tlen; ++ts) target[ts] = nucleotideCodeSubstitutionMatrix.get_dna_acid_map(_dna_d[ts]);

    mm_tbuf_t *tbuf = mm_tbuf_init();
    void *km = tbuf->km;
    int8_t q = -(_open_gap_penalty - _extend_gap_penalty);
    int8_t e = - _extend_gap_penalty;
    //q = 2;
    int8_t q2 = q+e*5;
    //std::cout << "q " << (int)q << " e " << (int)e << std::endl;
    int8_t m = 5;

#define __dp_code_block1 \
	z = _mm_load_si128(&s[t]); \
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
	b = _mm_add_epi8(_mm_load_si128(&y[t]), ut);     /* b <- y[r-1][t..t+15] + u[r-1][t..t+15] */ \
	x2t1= _mm_load_si128(&x2[t]); \
	tmp = _mm_srli_si128(x2t1, 15); \
	x2t1= _mm_or_si128(_mm_slli_si128(x2t1, 1), x21_); \
	x21_= tmp; \
	a2  = _mm_add_epi8(x2t1, vt1); \
	a2a = _mm_add_epi8(a2, _mm_load_si128(&acceptor[t]));

#define __dp_code_block2 \
	_mm_store_si128(&u[t], _mm_sub_epi8(z, vt1));    /* u[r][t..t+15] <- z - v[r-1][t-1..t+14] */ \
	_mm_store_si128(&v[t], _mm_sub_epi8(z, ut));     /* v[r][t..t+15] <- z - u[r-1][t..t+15] */ \
	tmp = _mm_sub_epi8(z, q_); \
	a = _mm_sub_epi8(a, tmp); \
	b = _mm_sub_epi8(b, tmp); \
	a2= _mm_sub_epi8(a2, _mm_sub_epi8(z, q2_));

    int r, t, qe = q + e, n_col_, *off = 0, *off_end = 0, tlen_, qlen_, last_st, last_en, max_sc, min_sc, long_thres, long_diff;
    int with_cigar = 1, approx_max=0;
    int32_t *H = 0, H0 = 0, last_H0_t = 0;
    uint8_t *qr, *sf, *mem, *mem2 = 0;
    __m128i q_, q2_, qe_, zero_, sc_mch_, sc_mis_, sc_N_, m1_;
    __m128i *u, *v, *x, *y, *x2, *s, *p = 0, *donor, *acceptor;

    ksw_extz_t ez1;
    memset(&ez1, 0, sizeof(ksw_extz_t));
    ksw_extz_t * ez = &ez1;
    ksw_reset_extz(ez);

    if (m <= 1 || qlen <= 0 || tlen <= 0 || q2 <= q + e) return;
    zero_   = _mm_set1_epi8(0);
    q_      = _mm_set1_epi8(q);
    q2_     = _mm_set1_epi8(q2);
    qe_     = _mm_set1_epi8(q + e);
    sc_mch_ = _mm_set1_epi8(nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[0][0]);
    sc_mis_ = _mm_set1_epi8(nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[0][1]);
    sc_N_   = nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[4][4] == 0? _mm_set1_epi8(-e) : _mm_set1_epi8(nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[4][4]);
    m1_     = _mm_set1_epi8(m - 1); // wildcard

    tlen_ = (tlen + 15) / 16;
    n_col_ = ((qlen < tlen? qlen : tlen) + 15) / 16 + 1;
    qlen_ = (qlen + 15) / 16;
//    for (t = 1, max_sc = mat[0], min_sc = mat[1]; t < m * m; ++t) {
//        max_sc = max_sc > mat[t]? max_sc : mat[t];
//        min_sc = min_sc < mat[t]? min_sc : mat[t];
//    }
    min_sc = nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[0][0];
    max_sc = nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[0][1];
    //std::cout << " line 183" << std::endl;
    if (-min_sc > 2 * (q + e)) return; // otherwise, we won't see any mismatches
    //std::cout << " line 185" << std::endl;
    long_thres = (q2 - q) / e - 1;
    if (q2 > q + e + long_thres * e)
        ++long_thres;
    long_diff = long_thres * e - (q2 - q);
    //std::cout << "long_thres " << long_thres << std::endl;
    //std::cout << "long_diff " << long_diff << std::endl;

    mem = (uint8_t*)kcalloc(km, tlen_ * 9 + qlen_ + 1, 16);
    u = (__m128i*)(((size_t)mem + 15) >> 4 << 4); // 16-byte aligned
    v = u + tlen_, x = v + tlen_, y = x + tlen_, x2 = y + tlen_;
    donor = x2 + tlen_, acceptor = donor + tlen_;
    s = acceptor + tlen_, sf = (uint8_t*)(s + tlen_), qr = sf + tlen_ * 16;
    memset(u,  -q - e,  tlen_ * 16 * 4); // this set u, v, x, y (because they are in the same array)
    memset(x2, -q2,     tlen_ * 16);
    if (!approx_max) {
        H = (int32_t*)kmalloc(km, tlen_ * 16 * 4);
        for (t = 0; t < tlen_ * 16; ++t) H[t] = KSW_NEG_INF;
    }
    if (with_cigar) {
        mem2 = (uint8_t*)kmalloc(km, ((size_t)(qlen + tlen - 1) * n_col_ + 1) * 16);
        p = (__m128i*)(((size_t)mem2 + 15) >> 4 << 4);
        off = (int*)kmalloc(km, (qlen + tlen - 1) * sizeof(int) * 2);
        off_end = off + qlen + tlen - 1;
    }

    for (t = 0; t < qlen; ++t) qr[t] = query[qlen - 1 - t];
    memcpy(sf, target, tlen);

    /*
    // set the donor and acceptor arrays. TODO: this assumes 0/1/2/3 encoding!
    if (flag & (KSW_EZ_SPLICE_FOR|KSW_EZ_SPLICE_REV)) {
        int semi_cost = flag&KSW_EZ_SPLICE_FLANK? -noncan/2 : 0; // GTr or yAG is worth 0.5 bit; see PMID:18688272
        memset(donor, -noncan, tlen_ * 16);
        for (t = 0; t < tlen - 4; ++t) {
            int can_type = 0; // type of canonical site: 0=none, 1=GT/AG only, 2=GTr/yAG
            if ((flag & KSW_EZ_SPLICE_FOR) && target[t+1] == 2 && target[t+2] == 3) can_type = 1; // GTr...
            if ((flag & KSW_EZ_SPLICE_REV) && target[t+1] == 1 && target[t+2] == 3) can_type = 1; // CTr...
            if (can_type && (target[t+3] == 0 || target[t+3] == 2)) can_type = 2;
            if (can_type) ((int8_t*)donor)[t] = can_type == 2? 0 : semi_cost;
        }
        memset(acceptor, -noncan, tlen_ * 16);
        for (t = 2; t < tlen; ++t) {
            int can_type = 0;
            if ((flag & KSW_EZ_SPLICE_FOR) && target[t-1] == 0 && target[t] == 2) can_type = 1; // ...yAG
            if ((flag & KSW_EZ_SPLICE_REV) && target[t-1] == 0 && target[t] == 1) can_type = 1; // ...yAC
            if (can_type && (target[t-2] == 1 || target[t-2] == 3)) can_type = 2;
            if (can_type) ((int8_t*)acceptor)[t] = can_type == 2? 0 : semi_cost;
        }
    }
*/

    for (r = 0, last_st = last_en = -1; r < qlen + tlen - 1; ++r) {
        int st = 0, en = tlen - 1, st0, en0, st_, en_;
        int8_t x1, x21, v1, *u8 = (int8_t*)u, *v8 = (int8_t*)v;
        uint8_t *qrr = qr + (qlen - 1 - r);
        __m128i x1_, x21_, v1_;
        // find the boundaries
        if (st < r - qlen + 1) st = r - qlen + 1;
        if (en > r) en = r;
        st0 = st, en0 = en;
        st = st / 16 * 16, en = (en + 16) / 16 * 16 - 1;
        // set boundary conditions
        /*if (st > 0) {
            if (st - 1 >= last_st && st - 1 <= last_en)
                x1 = ((int8_t*)x)[st - 1], x21 = ((int8_t*)x2)[st - 1], v1 = v8[st - 1]; // (r-1,s-1) calculated in the last round
            else x1 = -q - e, x21 = -q2, v1 = -q - e;
        } else {
            x1 = -q - e, x21 = -q2;
            v1 = r == 0? -q - e : r < long_thres? -e : r == long_thres? long_diff : 0;
        }
        if (en >= r) {
            ((int8_t*)y)[r] = -q - e;
            u8[r] = r == 0? -q - e : r < long_thres? -e : r == long_thres? long_diff : 0;
        }
*/
        if (st > 0) {
            if (st - 1 >= last_st && st - 1 <= last_en)
                x1 = ((int8_t*)x)[st - 1], x21 = ((int8_t*)x2)[st - 1], v1 = v8[st - 1]; // (r-1,s-1) calculated in the last round
            else x1 = -q - e, x21 = -q2, v1 = -q - e;
        } else {
            x1 = -q - e, x21 = -q2;
            v1 = r == 0? -q - e : r < long_thres? -e : r == long_thres? long_diff : 0;
        }
        if (en >= r) {
            ((int8_t*)y)[r] = -q - e;
            u8[r] = r == 0? -q - e : r < long_thres? -e : r == long_thres? long_diff : 0;
        }
        // loop fission: set scores first
        //if (!(flag & KSW_EZ_GENERIC_SC)) {
        if( true ){
            for (t = st0; t <= en0; t += 16) {
                __m128i sq, st, tmp, mask;
                sq = _mm_loadu_si128((__m128i*)&sf[t]);
                st = _mm_loadu_si128((__m128i*)&qrr[t]);
                mask = _mm_or_si128(_mm_cmpeq_epi8(sq, m1_), _mm_cmpeq_epi8(st, m1_));
                tmp = _mm_cmpeq_epi8(sq, st);
                #ifdef __SSE4_1__
                tmp = _mm_blendv_epi8(sc_mis_, sc_mch_, tmp);
                tmp = _mm_blendv_epi8(tmp,     sc_N_,   mask);
                #else
                tmp = _mm_or_si128(_mm_andnot_si128(tmp,  sc_mis_), _mm_and_si128(tmp,  sc_mch_));
                tmp = _mm_or_si128(_mm_andnot_si128(mask, tmp),     _mm_and_si128(mask, sc_N_));
                #endif
                _mm_storeu_si128((__m128i*)((int8_t*)s + t), tmp);
            }
        } else {
            for (t = st0; t <= en0; ++t)
                ((uint8_t*)s)[t] = nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[sf[t]][qrr[t]];
        }
        // core loop
        x1_  = _mm_cvtsi32_si128((uint8_t)x1);
        x21_ = _mm_cvtsi32_si128((uint8_t)x21);
        v1_  = _mm_cvtsi32_si128((uint8_t)v1);
        st_ = st / 16, en_ = en / 16;
        assert(en_ - st_ + 1 <= n_col_);
        if (!with_cigar) { // score only
            std::cout << "score only" << std::endl;
            for (t = st_; t <= en_; ++t) {
                __m128i z, a, b, a2, a2a, xt1, x2t1, vt1, ut, tmp;
                __dp_code_block1;
#ifdef __SSE4_1__
                z = _mm_max_epi8(z, a);
                z = _mm_max_epi8(z, b);
                z = _mm_max_epi8(z, a2a);
                __dp_code_block2; // save u[] and v[]; update a, b and a2
                _mm_store_si128(&x[t],  _mm_sub_epi8(_mm_max_epi8(a,  zero_), qe_));
                _mm_store_si128(&y[t],  _mm_sub_epi8(_mm_max_epi8(b,  zero_), qe_));
                tmp = _mm_load_si128(&donor[t]);
                _mm_store_si128(&x2[t], _mm_sub_epi8(_mm_max_epi8(a2, tmp), q2_));
#else
                tmp = _mm_cmpgt_epi8(a,  z);
				z = _mm_or_si128(_mm_andnot_si128(tmp, z), _mm_and_si128(tmp, a));
				tmp = _mm_cmpgt_epi8(b,  z);
				z = _mm_or_si128(_mm_andnot_si128(tmp, z), _mm_and_si128(tmp, b));
				tmp = _mm_cmpgt_epi8(a2a, z);
				z = _mm_or_si128(_mm_andnot_si128(tmp, z), _mm_and_si128(tmp, a2a));
				__dp_code_block2;
				tmp = _mm_cmpgt_epi8(a, zero_);
				_mm_store_si128(&x[t],  _mm_sub_epi8(_mm_and_si128(tmp, a),  qe_));
				tmp = _mm_cmpgt_epi8(b, zero_);
				_mm_store_si128(&y[t],  _mm_sub_epi8(_mm_and_si128(tmp, b),  qe_));
				tmp = _mm_load_si128(&donor[t]); // TODO: check if this is correct
				tmp = _mm_cmpgt_epi8(a2, tmp);
				tmp = _mm_or_si128(_mm_andnot_si128(tmp, tmp), _mm_and_si128(tmp, a2));
				_mm_store_si128(&x2[t], _mm_sub_epi8(tmp, q2_));
#endif
            }
        //} else if (!(flag&KSW_EZ_RIGHT)) { // gap left-alignment
        } else if ( true ) {
            __m128i *pr = p + r * n_col_ - st_;
            off[r] = st, off_end[r] = en;
            for (t = st_; t <= en_; ++t) {
                __m128i d, z, a, b, a2, a2a, xt1, x2t1, vt1, ut, tmp, tmp2;
                __dp_code_block1;
#ifdef __SSE4_1__
                d = _mm_and_si128(_mm_cmpgt_epi8(a, z), _mm_set1_epi8(1));       // d = a  > z? 1 : 0
                z = _mm_max_epi8(z, a);
                d = _mm_blendv_epi8(d, _mm_set1_epi8(2), _mm_cmpgt_epi8(b,  z)); // d = b  > z? 2 : d
                z = _mm_max_epi8(z, b);
                d = _mm_blendv_epi8(d, _mm_set1_epi8(3), _mm_cmpgt_epi8(a2a, z)); // d = a2 > z? 3 : d
                z = _mm_max_epi8(z, a2a);
#else // we need to emulate SSE4.1 intrinsics _mm_max_epi8() and _mm_blendv_epi8()
                tmp = _mm_cmpgt_epi8(a,  z);
				d = _mm_and_si128(tmp, _mm_set1_epi8(1));
				z = _mm_or_si128(_mm_andnot_si128(tmp, z), _mm_and_si128(tmp, a));
				tmp = _mm_cmpgt_epi8(b,  z);
				d = _mm_or_si128(_mm_andnot_si128(tmp, d), _mm_and_si128(tmp, _mm_set1_epi8(2)));
				z = _mm_or_si128(_mm_andnot_si128(tmp, z), _mm_and_si128(tmp, b));
				tmp = _mm_cmpgt_epi8(a2a, z);
				d = _mm_or_si128(_mm_andnot_si128(tmp, d), _mm_and_si128(tmp, _mm_set1_epi8(3)));
				z = _mm_or_si128(_mm_andnot_si128(tmp, z), _mm_and_si128(tmp, a2a));
#endif
                __dp_code_block2;
                tmp = _mm_cmpgt_epi8(a, zero_);
                _mm_store_si128(&x[t],  _mm_sub_epi8(_mm_and_si128(tmp, a),  qe_));
                d = _mm_or_si128(d, _mm_and_si128(tmp, _mm_set1_epi8(0x08))); // d = a > 0? 1<<3 : 0
                tmp = _mm_cmpgt_epi8(b, zero_);
                _mm_store_si128(&y[t],  _mm_sub_epi8(_mm_and_si128(tmp, b),  qe_));
                d = _mm_or_si128(d, _mm_and_si128(tmp, _mm_set1_epi8(0x10))); // d = b > 0? 1<<4 : 0

                tmp2 = _mm_load_si128(&donor[t]);
                tmp = _mm_cmpgt_epi8(a2, tmp2);
#ifdef __SSE4_1__
                tmp2 = _mm_max_epi8(a2, tmp2);
#else
                tmp2 = _mm_or_si128(_mm_andnot_si128(tmp, tmp2), _mm_and_si128(tmp, a2));
#endif
                _mm_store_si128(&x2[t], _mm_sub_epi8(tmp2, q2_));
                d = _mm_or_si128(d, _mm_and_si128(tmp, _mm_set1_epi8(0x20)));
                _mm_store_si128(&pr[t], d);
            }
        } else { // gap right-alignment
            __m128i *pr = p + r * n_col_ - st_;
            off[r] = st, off_end[r] = en;
            for (t = st_; t <= en_; ++t) {
                __m128i d, z, a, b, a2, a2a, xt1, x2t1, vt1, ut, tmp, tmp2;
                __dp_code_block1;
#ifdef __SSE4_1__
                d = _mm_andnot_si128(_mm_cmpgt_epi8(z, a), _mm_set1_epi8(1));    // d = z > a?  0 : 1
                z = _mm_max_epi8(z, a);
                d = _mm_blendv_epi8(_mm_set1_epi8(2), d, _mm_cmpgt_epi8(z, b));  // d = z > b?  d : 2
                z = _mm_max_epi8(z, b);
                d = _mm_blendv_epi8(_mm_set1_epi8(3), d, _mm_cmpgt_epi8(z, a2a)); // d = z > a2? d : 3
                z = _mm_max_epi8(z, a2a);
#else // we need to emulate SSE4.1 intrinsics _mm_max_epi8() and _mm_blendv_epi8()
                tmp = _mm_cmpgt_epi8(z, a);
				d = _mm_andnot_si128(tmp, _mm_set1_epi8(1));
				z = _mm_or_si128(_mm_and_si128(tmp, z), _mm_andnot_si128(tmp, a));
				tmp = _mm_cmpgt_epi8(z, b);
				d = _mm_or_si128(_mm_and_si128(tmp, d), _mm_andnot_si128(tmp, _mm_set1_epi8(2)));
				z = _mm_or_si128(_mm_and_si128(tmp, z), _mm_andnot_si128(tmp, b));
				tmp = _mm_cmpgt_epi8(z, a2a);
				d = _mm_or_si128(_mm_and_si128(tmp, d), _mm_andnot_si128(tmp, _mm_set1_epi8(3)));
				z = _mm_or_si128(_mm_and_si128(tmp, z), _mm_andnot_si128(tmp, a2a));
#endif
                __dp_code_block2;
                tmp = _mm_cmpgt_epi8(zero_, a);
                _mm_store_si128(&x[t],  _mm_sub_epi8(_mm_andnot_si128(tmp, a),  qe_));
                d = _mm_or_si128(d, _mm_andnot_si128(tmp, _mm_set1_epi8(0x08))); // d = a > 0? 1<<3 : 0
                tmp = _mm_cmpgt_epi8(zero_, b);
                _mm_store_si128(&y[t],  _mm_sub_epi8(_mm_andnot_si128(tmp, b),  qe_));
                d = _mm_or_si128(d, _mm_andnot_si128(tmp, _mm_set1_epi8(0x10))); // d = b > 0? 1<<4 : 0

                tmp2 = _mm_load_si128(&donor[t]);
                tmp = _mm_cmpgt_epi8(tmp2, a2);
#ifdef __SSE4_1__
                tmp2 = _mm_max_epi8(tmp2, a2);
#else
                tmp2 = _mm_or_si128(_mm_andnot_si128(tmp, a2), _mm_and_si128(tmp, tmp2));
#endif
                _mm_store_si128(&x2[t], _mm_sub_epi8(tmp2, q2_));
                d = _mm_or_si128(d, _mm_andnot_si128(tmp, _mm_set1_epi8(0x20))); // d = a > 0? 1<<5 : 0
                _mm_store_si128(&pr[t], d);
            }
        }
        if (!approx_max) { // find the exact max with a 32-bit score array
            int32_t max_H, max_t;
            // compute H[], max_H and max_t
            if (r > 0) {
                int32_t HH[4], tt[4], en1 = st0 + (en0 - st0) / 4 * 4, i;
                __m128i max_H_, max_t_;
                max_H = H[en0] = en0 > 0? H[en0-1] + u8[en0] : H[en0] + v8[en0]; // special casing the last element
                max_t = en0;
                max_H_ = _mm_set1_epi32(max_H);
                max_t_ = _mm_set1_epi32(max_t);
                for (t = st0; t < en1; t += 4) { // this implements: H[t]+=v8[t]-qe; if(H[t]>max_H) max_H=H[t],max_t=t;
                    __m128i H1, tmp, t_;
                    H1 = _mm_loadu_si128((__m128i*)&H[t]);
                    t_ = _mm_setr_epi32(v8[t], v8[t+1], v8[t+2], v8[t+3]);
                    H1 = _mm_add_epi32(H1, t_);

                    int8_t * ta = ( int8_t*)(&H1);
//                    for ( i = 0; i < 4; ++i ) {
//                        std::cout << "298 i " << r - (t*4+i - st0) << " j " << t*4+i - st0 << " H[t] " << H[t] << std::endl;
//                    }

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
                    if (max_H < HH[i]) max_H = HH[i], max_t = tt[i] + i;
                for (; t < en0; ++t) { // for the rest of values that haven't been computed with SSE
                    H[t] += (int32_t)v8[t];
                    //std::cout << "319 i " << t << " r " << r << " j " << r-t << " H[t] " << H[t] << std::endl;
                    if (H[t] > max_H)
                        max_H = H[t], max_t = t;
                }
            } else H[0] = v8[0] - qe, max_H = H[0], max_t = 0; // special casing r==0
            // update ez
            if (en0 == tlen - 1 && H[en0] > ez->mte)
                ez->mte = H[en0], ez->mte_q = r - en;
            if (r - st0 == qlen - 1 && H[st0] > ez->mqe)
                ez->mqe = H[st0], ez->mqe_t = st0;
            // songbaoxing if (ksw_apply_zdrop(ez, 1, max_H, r, max_t, zdrop, 0)) break;
            if (r == qlen + tlen - 2 && en0 == tlen - 1)
                ez->score = H[tlen - 1];
        } else { // find approximate max; Z-drop might be inaccurate, too.
            if (r > 0) {
                if (last_H0_t >= st0 && last_H0_t <= en0 && last_H0_t + 1 >= st0 && last_H0_t + 1 <= en0) {
                    int32_t d0 = v8[last_H0_t];
                    int32_t d1 = u8[last_H0_t + 1];
                    if (d0 > d1) H0 += d0;
                    else H0 += d1, ++last_H0_t;
                } else if (last_H0_t >= st0 && last_H0_t <= en0) {
                    H0 += v8[last_H0_t];
                } else {
                    ++last_H0_t, H0 += u8[last_H0_t];
                }
            } else H0 = v8[0] - qe, last_H0_t = 0;
            // songbaoxing if ((flag & KSW_EZ_APPROX_DROP) && ksw_apply_zdrop(ez, 1, H0, r, last_H0_t, zdrop, 0)) break;
            if (r == qlen + tlen - 2 && en0 == tlen - 1)
                ez->score = H0;
        }
        last_st = st, last_en = en;
        //for (t = st0; t <= en0; ++t) printf("(%d,%d)\t(%d,%d,%d,%d)\t%d\n", r, t, ((int8_t*)u)[t], ((int8_t*)v)[t], ((int8_t*)x)[t], ((int8_t*)y)[t], H[t]); // for debugging
    }
    kfree(km, mem);
    if (!approx_max) kfree(km, H);
//    if (with_cigar) { // backtrack
//        int rev_cigar = !!(flag & KSW_EZ_REV_CIGAR);
//        if (!ez->zdropped && !(flag&KSW_EZ_EXTZ_ONLY))
//            ksw_backtrack(km, 1, rev_cigar, long_thres, (uint8_t*)p, off, off_end, n_col_*16, tlen-1, qlen-1, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
//        else if (ez->max_t >= 0 && ez->max_q >= 0)
//            ksw_backtrack(km, 1, rev_cigar, long_thres, (uint8_t*)p, off, off_end, n_col_*16, ez->max_t, ez->max_q, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
//        kfree(km, mem2); kfree(km, off);
//    }

    int _dna_q_i=0;
    int _dna_d_i=0;
    if( ez->mte > ez->mqe ){
        //std::cout << "line 352" << std::endl;
        ksw_backtrack(km, (uint8_t*)p, off, off_end, n_col_*16, tlen-1, ez->mte_q, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
        _dna_q_i=ez->mte_q;
        _dna_d_i=tlen-1;
        //std::cout << "_dna_q_i " << _dna_q_i << " _dna_d_i " << _dna_d_i << " ez->mte " << ez->mte << std::endl;
    }else{
        ksw_backtrack(km, (uint8_t*)p, off, off_end, n_col_*16, ez->mqe_t, qlen-1, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
        //std::cout << "line 357" << std::endl;
        _dna_q_i=qlen-1;
        _dna_d_i=ez->mqe_t;
        //std::cout << "_dna_q_i " << _dna_q_i << " _dna_d_i " << _dna_d_i << " ez->mqe " << ez->mqe << std::endl;
    }

    // due to the order reason, it is difficult to put the result into vector directly
    for (int i = 0; i<ez->n_cigar; ++i){ // IMPORTANT: this gives the CIGAR in the aligned regions. NO soft/hard clippings!
        // the cigar is in reverse order, and that is what I want
        int cigarLength = ez->cigar[i]>>4;
        int cigarType = ez->cigar[i]&0xf;
        //printf("%d%c", cigarLength, "MIDNSH"[cigarType]);
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
//    std::cout << std::endl;
    kfree(km, mem2); kfree(km, off);
    kfree(km, ez1.cigar);
    mm_tbuf_destroy(tbuf);
    delete[] query;
    delete[] target;
}

