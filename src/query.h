#ifndef H_QUERY
#define H_QUERY

#ifdef USE_PARASAIL
#include <parasail.h>
#endif

#include "seq.h"
#include "priority_queue.h"
#include "hit.h"
#include "db.h"
#include "options.h"
#include "ssw.h"


typedef enum {
    QUERY_FIND_SUFFIXES,
    QUERY_GET_ALIGNMENTS,
    QUERY_DO_ALIGNMENTS,
    QUERY_DO_STATISTICS,
    QUERY_DONE
} query_state_t;

typedef struct {
    seq_t* query;
    db_t* db;
    options_t* opt;
    pqueue_t* pq;
    hitlist_t* hl;
    int hit_counter;
    int alignment_counter;
    int min_alignment_score;
    query_state_t state;

    int iterations;

    s_profile* ssw_profile;

#ifdef USE_PARASAIL
    parasail_profile_t* parasail_profile;
    int use_parasail;
#endif


} query_t;

query_t* query_alloc(seq_t* s, db_t* db, options_t* opt);
void query_free(query_t* q);
int query_done(query_t* q);
int query_find_suffixes(query_t* q);
int query_fetch_alignment_batch(query_t* q);
int query_perform_alignments(query_t* q);
int query_perform_statistics(query_t* q);
int query_stop_looking(query_t* q);

#ifdef USE_PARASAIL
int perform_alignment_with_profile(hit_t *h, s_profile *profile, parasail_profile_t* parasail_profile, db_t *db, options_t *opt);
#else
int perform_alignment_with_profile(hit_t *h, s_profile *profile, db_t *db, options_t *opt);
#endif

#endif

