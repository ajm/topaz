#include <stdlib.h>

#include "suffix.h"
#include "priority_queue.h"


void _suffix_deallocator(void *payload) {
    suffix_free((suffix_t*) payload);
}

pqueue_t* pqueue_alloc(db_t *db, options_t *opt) {
    pqueue_t *pq;

    pq = (pqueue_t*) malloc(sizeof(pqueue_t));
    pq->t = rb_alloc();
    pq->db = db;
    pq->opt = opt;

    return pq;
}

void pqueue_free(pqueue_t *pq) {
    if(pq) {
        if(pq->t) {
            rb_free(pq->t, _suffix_deallocator);
        }

        free(pq);
    }
}

void pqueue_add(pqueue_t *pq, suffix_t *s) {
    rb_put(pq->t, suffix_priority(s), (void*) s);
}

int pqueue_next(pqueue_t *pq, saidx64_t *sa_pos, saidx64_t *query_offset, int *gapless_score, int *gapless_length) {
    suffix_t *suf;

    if((suf = (suffix_t*) rb_popmax(pq->t)) == NULL) {  // get suffix with max LCP
        return -1;       
    }

    *sa_pos = suf->position;
    *query_offset = suf->offset;
    *gapless_score = suf->gapless_score;
    *gapless_length = suf->gapless_length;

    if(suffix_increment(suf, pq->db, pq->opt) == 0 && suffix_active(suf, pq->opt->min_hit_length)) {   // increment index in the correct direction
        rb_put(pq->t, suffix_priority(suf), (void*) suf);
    }
    else {
        suffix_free(suf);
    }

    return 0;
}

