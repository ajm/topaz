#ifndef H_PRIORITY_QUEUE
#define H_PRIORITY_QUEUE

#include "rbtree.h"
#include "db.h"
#include "options.h"
#include "divsufsort64.h"
#include "suffix.h"


typedef struct {
    rbtree_t *t;
    db_t *db;
    options_t *opt;

} pqueue_t;


pqueue_t* pqueue_alloc(db_t *db, options_t *opt);
void pqueue_free(pqueue_t *pq);
void pqueue_add(pqueue_t *pq, suffix_t *s);
int pqueue_next(pqueue_t *pq, 
        saidx64_t *sa_pos, 
        saidx64_t *query_offset, 
        int *gapless_score, 
        int *gapless_length);

#endif

