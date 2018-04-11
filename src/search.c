#define _GNU_SOURCE
#include <stdio.h>
#include <stddef.h>
#include <sys/types.h>

#include "divsufsort64.h"

#include "substitution.h"
#include "lexicographical.h"
#include "options.h"
#include "generic.h"
#include "seq.h"
#include "db.h"


#define min2(a,b) ((a) < (b) ? (a) : (b))

static
int
_compare(const char *T, saidx64_t Tsize,
         const char *P, saidx64_t Psize,
         saidx64_t suf, saidx64_t *match) {

    saidx64_t i, j;
    saint_t r;

    for(i = suf + *match, j = *match, r = 0; (i < Tsize) && (j < Psize) && ((r = T[i] - P[j]) == 0); ++i, ++j) {}
    *match = j;
    return (r == 0) ? -(j != Psize) : r;
}

int sa_lcp(db_t *db, seq_t *q, int qoffset, saidx64_t position) {
    saidx64_t idx = db->SA[position];
    char *T = &db->T[idx];
    char *Q = seq_ptr(q) + qoffset;
    int Qsize = seq_len(q) - qoffset;
    int length = 0;
    
    while(*T++ == *Q++) {
        if(++length >= Qsize) {
            break;
        }
    }

    return length;
}

int sa_gapless(db_t *db, seq_t *q, int qoffset, saidx64_t position, options_t* opt, int* gapless_score, int* gapless_length) {
    saidx64_t idx = db->SA[position];
    char *T = &db->T[idx];
    char *Q = seq_ptr(q) + qoffset;
    int max_length = min2(db->Tsize - idx, seq_len(q) - qoffset);
    int i, ind, rawscore = 0;
    int argmax_score = 0, argmax_len = 0;

    // first letter must match
    //if(*T != *Q)
    //    return 0;

    // require seed be at least 3 amino acids long
    // and an exact match
    if((max_length < 3) || (T[0] != Q[0])) // || (T[1] != Q[1]) || (T[2] != Q[2]))
        return 0;

    for(i = 0; i < max_length; ++i) {
        if(T[i] == '\n') {
            break;
        }

        // debug
        //if(T[i] != Q[i]) {
        //    break;
        //}

#ifdef USELEXICOGRAPHICAL
        ind = (aa_table[(int)T[i]] * 25) + aa_table[(int) upper_aa_only(Q[i])];
#else
        ind = (aa_table2[(int)T[i]] * 25) + aa_table2[(int) upper_aa_only(Q[i])];
#endif

        rawscore += (*opt->substitution)[ind];

        if(rawscore >= argmax_score) {
            argmax_score = rawscore;
            argmax_len = i;
        }
    }

    *gapless_score = argmax_score;
    *gapless_length = argmax_len;

    return 0;

    //return argmax_score < opt->min_hit_score ? 0 : argmax_score; // argmax_len;
}

int sa_hsp(db_t *db, seq_t *q, int qoffset, saidx64_t position, options_t* opt, int* gapless_score, int* gapless_length) {
    saidx64_t idx = db->SA[position];
    char *T = &db->T[idx];
    char *Q = seq_ptr(q) + qoffset;
    int max_length_right = min2(db->Tsize - idx, seq_len(q) - qoffset);
    int max_length_left = min2(idx, qoffset);
    int i, ind;
    int rawscore_right = 0, rawscore_left = 0;
    int argmax_score_right = 0, argmax_len_right = 0;
    int argmax_score_left = 0, argmax_len_left = 0;

    // require seed be at least 3 amino acids long
    // and an exact match
    //if((max_length_right < 3) || (T[0] != Q[0])) { // || (T[1] != Q[1]) || (T[2] != Q[2])) {
    if(T[0] != Q[0]) {
        return 0;
    }

    // go right ----------->
    for(i = 0; i < max_length_right; i++) {
        if(T[i] == '\n') {
            break;
        }

#ifdef USELEXICOGRAPHICAL
        ind = (aa_table[(int)T[i]] * 25) + aa_table[(int) upper_aa_only(Q[i])];
#else
        ind = (aa_table2[(int)T[i]] * 25) + aa_table2[(int) upper_aa_only(Q[i])]; // XXX
#endif

        rawscore_right += (*opt->substitution)[ind];

        if(rawscore_right >= argmax_score_right) {
            argmax_score_right = rawscore_right;
            argmax_len_right = i;
        }
/*
        if(rawscore_right < 1) {
            break;
        }
*/
    }

    // <------------ go left
    for(i = 1; i < max_length_left; i++) {
        if(T[-i] == '\n') {
            break;
        }

#ifdef USELEXICOGRAPHICAL
        ind = (aa_table[(int)T[-i]] * 25) + aa_table[(int) upper_aa_only(Q[-i])];
#else
        ind = (aa_table2[(int)T[-i]] * 25) + aa_table2[(int) upper_aa_only(Q[-i])]; // XXX
#endif

        rawscore_left += (*opt->substitution)[ind];
 
        if(rawscore_left >= argmax_score_left) {
            argmax_score_left = rawscore_left;
            argmax_len_left = i;
        }
/*
        if(rawscore_left < 1) {
            break;
        }
*/
    }


    *gapless_score = argmax_score_left + argmax_score_right;
    *gapless_length = argmax_len_left + argmax_len_right;

    return 0;
}

#ifdef USE_LCP
static
int _compare_lcp(const char *T, saidx64_t Tsize,
                 const char *P, saidx64_t Psize,
                 const saidx64_t *SA, saidx64_t *match,
                 const saidx64_t *LCPLC, const saidx64_t *LCPCR,
                 saidx64_t lmatch, saidx64_t rmatch,
                 saidx64_t index) {

    saidx64_t current_lcp;

    //fprintf(stderr, "index = %lld\n", index);

    if(lmatch == rmatch) {
        *match = lmatch;
        return _compare(T, Tsize, P, Psize, SA[index], match);
    }
    else if(lmatch > rmatch) {
        current_lcp = LCPLC[index];
        if(current_lcp > lmatch) {
            *match = lmatch;
            return -1;
        }
        if(current_lcp < lmatch) {
            *match = current_lcp;
            return 1;
        }
    }
    else {
        current_lcp = LCPCR[index];
        if(current_lcp > rmatch) {
            *match = rmatch;
            return 1;
        }
        if(current_lcp < rmatch) {
            *match = current_lcp;
            return -1;
        }
    }

    *match = lmatch;
    return _compare(T, Tsize, P, Psize, SA[index], match);
}
#endif

int sa_search(db_t *db, seq_t *q, saidx64_t offset, saidx64_t length, saidx64_t *idx, saidx64_t *num_matches, saidx64_t *last_position) {
    char *T;
    saidx64_t *SA;
    saidx64_t Tsize;
    
#ifdef USE_LCP
    saidx64_t *LCPLC;
    saidx64_t *LCPCR;
#endif

    char *P;
    saidx64_t Psize;

    saidx64_t size, lsize, rsize, half;
    saidx64_t match, lmatch, rmatch;
    saidx64_t llmatch, lrmatch, rlmatch, rrmatch;
    saidx64_t i, j, k;
    saint_t r;


    T = db->T;
    SA = db->SA;
    Tsize = db->Tsize;

#ifdef USE_LCP
    LCPLC = db->LCPLC;
    LCPCR = db->LCPCR;
#endif

    P = seq_ptr(q) + offset;
    // if the length is specified, use that
    // otherwise the length is the length of the suffix starting from offset
    Psize = length > 0 ? length : seq_len(q) - offset;


    if(idx != NULL) { 
        *idx = -1; 
    }

    if((T == NULL) || (P == NULL) || (SA == NULL) || (Tsize <= 0) || (Psize <= 0)) { 
        return -1;
    }

    for(i = j = k = 0, match = lmatch = rmatch = 0, size = Tsize, half = size >> 1;
        0 < size;
        size = half, half >>= 1) {

#ifndef USE_LCP
        match = MIN(lmatch, rmatch);
        r = _compare(T, Tsize, P, Psize, SA[i + half], &match);
#else
        r = _compare_lcp(T, Tsize, P, Psize, SA, &match, LCPLC, LCPCR, lmatch, rmatch, i + half);
#endif

        if(r < 0) { // upper half of range (i+half+1 : i+size)
            i += half + 1;
            half -= (size & 1) ^ 1;
            lmatch = match;
        } 
        else if(r > 0) { // lower half of range (i : i+half)
            rmatch = match;
        }
        else { // perfect match
            lsize = half, j = i, rsize = size - half - 1, k = i + half + 1;

            // search for lowest index that is a match
            for(llmatch = lmatch, lrmatch = match, half = lsize >> 1; 0 < lsize; lsize = half, half >>= 1) {

#ifndef USE_LCP
                lmatch = MIN(llmatch, lrmatch);
                r = _compare(T, Tsize, P, match, SA[j + half], &lmatch);
#else
                r = _compare_lcp(T, Tsize, P, match, SA, &lmatch, LCPLC, LCPCR, llmatch, lrmatch, j + half);
#endif

                if(r < 0) { // upper half
                    j += half + 1;
                    half -= (lsize & 1) ^ 1;
                    llmatch = lmatch;
                } 
                else { // lower half
                    lrmatch = lmatch;
                }
            }

            // search for highest index that is a match
            for(rlmatch = match, rrmatch = rmatch, half = rsize >> 1; 0 < rsize; rsize = half, half >>= 1) {
#ifndef USE_LCP
                rmatch = MIN(rlmatch, rrmatch);
                r = _compare(T, Tsize, P, match, SA[k + half], &rmatch);
#else
                r = _compare_lcp(T, Tsize, P, match, SA, &rmatch, LCPLC, LCPCR, rlmatch, rrmatch, k + half);
#endif
                if(r <= 0) { // upper half
                    k += half + 1;
                    half -= (rsize & 1) ^ 1;
                    rlmatch = rmatch;
                }
                else { // lower half
                    rrmatch = rmatch;
                }
            }

            break;
        }
    }


    // idx is pointer to first best match 
    //      if more than one match, point to first one in SA (j)
    //      otherwise exact match (i)
    if(idx != NULL) { 
        *idx = (0 < (k - j)) ? j : i; 
    }

    if(num_matches != NULL) {
        *num_matches = k - j;
    }

    if(last_position != NULL) {
        *last_position = match;
    }

    return 0;
}

