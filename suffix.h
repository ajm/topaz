#ifndef H_SUFFIX
#define H_SUFFIX

#include <stdint.h>

#include "seq.h"
#include "db.h"
#include "divsufsort64.h"
#include "options.h"


typedef struct {
    // suffix definition
    seq_t* seq;             // pointer to protein (read only)
    int offset;             // offset into protein (read only)

    // relating to the suffix array
    saidx64_t position;     // SA index

    // necessary for binary walk

    // necessary for linear walk
    int direction;          // SA walk direction (+1,-1) (read only)
    int gapless_length;         // length of LCP at current SA index
    int gapless_score;              // number of identical matches in SA

} suffix_t;

suffix_t* suffix_alloc(seq_t* seq, db_t *db, int offset, saidx64_t position, int direction, int gapless_score, int gapless_length);
suffix_t* suffix_duplicate(suffix_t* suf);
void suffix_free(suffix_t* suf);
int suffix_active(suffix_t* suf, int min_gapless_length);
int suffix_increment(suffix_t *suf, db_t* db, options_t* opt);
uint64_t suffix_priority(suffix_t *suf);

/*
#define SIZE_GRANULARITY 10

typedef struct {
    uint64_t page_offset;       // page offset is the key for grouping (unique)
    uint64_t size_priority;     // size of suffix_list divided by SIZE_GRANULARITY (non-unique)
    vector_t* suffix_list;      // length is the priority

    // input key is pages (i.e. range of memory)
    // output key is length of suffix list + whether page is already in memory

} access_t;

access_t* access_alloc();
void access_free(access_t* acc);
void access_add(access_t* acc, suffix_t* suf);
*/

#endif

