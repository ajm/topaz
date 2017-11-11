#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include "hit.h"
#include "rbtree.h"
#include "vector.h"
#include "divsufsort64.h"
#include "options.h"


hit_t* hit_alloc() {
    hit_t *h;

    h = (hit_t*) malloc(sizeof(hit_t));

    h->db_index = -1;
    //h->exactmatches = vector_alloc();
    h->max_score = 0;

    h->is_ready = 0;
    h->is_aligned = 0;
    h->alignment_rank = -1;
    h->rawscore = -1;
    h->mismatch = -1;
    h->length = -1;
    h->gapopen = -1;
    h->identity = INFINITY;
    h->evalue = INFINITY;
    h->bitscore = INFINITY;
    h->qstart = -1;
    h->qend = -1;
    h->sstart = -1;
    h->send = -1;
    h->effective_length = -1;
    h->suffix_counter = 0;

    return h;
}

void hit_free(hit_t *h) {
    if(h) {
        //vector_free(h->exactmatches);
        free(h);
    }
}

int hit_maxscore(hit_t *h) {
    return h->max_score;
}

int _compare_exactmatch(void *item1, void *item2) {
    exactmatch_t *em1, *em2;

    em1 = (exactmatch_t*) item1;
    em2 = (exactmatch_t*) item2;

    if(em1->subject_stop < em2->subject_stop)
        return 1;

    if(em1->subject_stop > em2->subject_stop)
        return -1;

    return 0;
}

int _exactmatch_subsumed(exactmatch_t *lhs, exactmatch_t *rhs) {
    return (lhs->query_stop == rhs->query_stop) \
        && (lhs->subject_stop == rhs->subject_stop);
}

exactmatch_t* exactmatch_alloc(saidx64_t query_offset, saidx64_t subject_offset, saidx64_t lcp_length) {
    exactmatch_t *e;
    
    e = (exactmatch_t*) malloc(sizeof(exactmatch_t));

    e->query_start = query_offset;
    e->query_stop = query_offset + lcp_length;

    e->subject_start = subject_offset;
    e->subject_stop = subject_offset + lcp_length;

    e->length = lcp_length;
    //e->rawscore = rawscore;

    return e;
}

void exactmatch_free(exactmatch_t *e) {
    if(e) {
        free(e);
    }
}

int hitlist_len(hitlist_t* hl) {
    assert(rb_size(hl->id_tree) == rb_size(hl->priority_tree));
    return rb_size(hl->id_tree);
}

void _hit_deallocator(void *payload) {
    hit_free((hit_t*) payload);
}

hitlist_t* hitlist_alloc(options_t *opt) {
    hitlist_t* h;

    h = (hitlist_t*) malloc(sizeof(hitlist_t));
    h->id_tree = rb_alloc();
    h->priority_tree = rb_alloc();
    h->opt = opt;

    return h;
}

void hitlist_free(hitlist_t *hl) {
    if(hl) {
        if(hl->id_tree)
            rb_free(hl->id_tree, _hit_deallocator);

        if(hl->priority_tree)
            rb_free(hl->priority_tree, NULL);

        free(hl);
    }
}

hit_t* hitlist_add(hitlist_t* hl, db_t* db, saidx64_t SAindex, int gapless_score, int gapless_length, saidx64_t suffix_offset) {
    hit_t *h, *tmp;
    int max_hitlist_length = hl->opt->number_of_alignments;
    saidx64_t index;

    if(rb_size(hl->id_tree) >= max_hitlist_length) {    // are we at capacity?
        h = (hit_t*) rb_getmin(hl->priority_tree);      // get worst hit by score
        if(h->gapless_score >= gapless_score) {         // check if score is worse than worst in list
            return h;                                   // i don't want to return NULL, because this is legit
        }
    }

    if(db_index(db, db->SA[SAindex], &index) < 0) {
        fprintf(stderr, "Error: internal error in protein lookup\n");
        _exit(EXIT_FAILURE);
    }

    if((h = (hit_t*) rb_get(hl->id_tree, index)) == NULL) {
        // protein not in hitlist
/*
        if(rb_size(hl->id_tree) >= max_hitlist_length) {    // are we at capacity?
            h = (hit_t*) rb_getmin(hl->priority_tree);      // get worst hit by score
            if(h->gapless_score >= gapless_score) {         // check if score is worse than worst in list
                return h;                                   // i don't want to return NULL, because this is legit
            }
        }
*/
        // not at capacity or improvement on worst

        h = hit_alloc();

        h->db_index = index;

        h->gapless_score = gapless_score;
        h->gapless_length = gapless_length;
        h->suffix_counter = 1;
        h->suffix_offset = suffix_offset;
        h->suffix_multiple = 0;

        rb_put(hl->id_tree, index, (void*) h);
        h->priority_node = rb_put(hl->priority_tree, gapless_score, (void*) h);

        if(rb_size(hl->id_tree) > max_hitlist_length) {     // are we over capacity?
            tmp = (hit_t*) rb_popmin(hl->priority_tree);    // pop worst result by score
            rb_rm(hl->id_tree, tmp->db_index);              // delete rb_node by db_index
            hit_free(tmp);                                  // free hit
        }

        return h;
    }

    // protein already in hitlist
    // test if score for protein has gotten better, re-insert into rb_tree
    if(h->gapless_score < gapless_score) {
        h->gapless_score = gapless_score;
        h->gapless_length = gapless_length;
        rb_reput(hl->priority_tree, gapless_score, h->priority_node);
    }

    h->suffix_counter++;
    if(!h->suffix_multiple) {
        if(suffix_offset != h->suffix_offset) {
            h->suffix_multiple = 1;
        }
    }

    return NULL;
}

/*
void hitlist_print(hitlist_t *hl) {
    traversal_t* tr;
    hit_t* h;

    if(rb_empty(hl->t)) {
        printf("No hits!\n");
        return;
    }

    tr = traversal_alloc(hl->t);

    while((h = (hit_t*) traversal_next(tr)) != NULL) {
        printf("%llu %d\n", h->db_index, h->votes);
    }

    traversal_free(tr);
}
*/

traversal_t* hitlist_traversal(hitlist_t* hl) {
    return traversal_alloc(hl->priority_tree);
}

