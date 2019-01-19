//
// Created by Baoxing song on 2018-12-26.
//

#ifndef GEAN_HENGALIGN_H
#define GEAN_HENGALIGN_H

#include <string>
#include <stack>
#include <immintrin.h>
#include <map>
#include "../impl/impl.h"
#include "../util/nucleotideCodeSubstitutionMatrix.h"
#include "../util/parameters.h"
#include <stdio.h>


#define KSW_NEG_INF -0x40000000

#define kcalloc(km, count, size) calloc((count), (size))
#define krealloc(km, ptr, size) realloc((ptr), (size))
#define kmalloc(km, size) malloc((size))
#define kfree(km, ptr) free((ptr))

typedef struct {
    uint32_t max:31, zdropped:1;
    int max_q, max_t;      // max extension coordinate
    int mqe, mqe_t;        // max score when reaching the end of query
    int mte, mte_q;        // max score when reaching the end of target
    int score;             // max score reaching both ends; may be KSW_NEG_INF
    int m_cigar, n_cigar;
    int reach_end;
    uint32_t *cigar;
} ksw_extz_t;

static inline void ksw_reset_extz(ksw_extz_t *ez){
    ez->max_q = ez->max_t = ez->mqe_t = ez->mte_q = -1;
    ez->max = 0, ez->score = ez->mqe = ez->mte = KSW_NEG_INF;
    ez->n_cigar = 0, ez->zdropped = 0, ez->reach_end = 0;
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

void *km_init(void);

void hengAlign(const std::string& _dna_q, const std::string& _dna_d, std::stack<char> & SQ, std::stack<char> & SD,
               const int8_t & _open_gap_penalty, const int8_t & _extend_gap_penalty,
               NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix);
#endif //GEAN_HENGALIGN_H
