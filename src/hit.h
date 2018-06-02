#ifndef H_HIT
#define H_HIT

#include <stdint.h>

#include "divsufsort64.h"
#include "db.h"
#include "rbtree.h"
#include "options.h"


typedef struct {
    saidx64_t query_start;
    saidx64_t query_stop;

    saidx64_t subject_start;
    saidx64_t subject_stop;

    saidx64_t length;
    int rawscore;

} exactmatch_t;

exactmatch_t* exactmatch_alloc(saidx64_t query_offset, saidx64_t subject_offset, saidx64_t lcp_length); //, int rawscore);
void exactmatch_free(exactmatch_t *e);

typedef struct {
    // sa search
    uint64_t db_index;  // index of protein in db->index (array of position_t objects)

    // alignment status
    int is_aligned;     // already aligned?
    int alignment_rank;

    // alignment stats
    int rawscore;
    int mismatch;
    int length;
    int gapopen;
    double identity;
    double bitscore;
    double evalue;
    int qstart;
    int qend;
    int sstart;
    int send;
    int effective_length;

    int gapless_score;
    int gapless_length;
    int suffix_counter;
    int suffix_offset;
    int suffix_multiple;

    double time_stats;
    double time_align;

    rbnode_t *priority_node;

} hit_t;

hit_t* hit_alloc();
void hit_free(hit_t* h);
int hit_add_exactmatch(hit_t *h, exactmatch_t *e, int min_distance, int max_distance);
int hit_add_exactmatch_no_merge(hit_t *h, exactmatch_t *e, int min_distance, int max_distance);
int hit_maxscore(hit_t *h);

typedef struct {
    rbtree_t *id_tree;
    rbtree_t *priority_tree;
    options_t *opt;
    int max_length;

} hitlist_t;

hitlist_t* hitlist_alloc(options_t *opt);
void hitlist_free(hitlist_t* hl);
//void hitlist_prune(hitlist_t* hl, int min_votes, int max_hits);
//void hitlist_increment(hitlist_t* hl, uint64_t db_index);
//hit_t* hitlist_increment2(hitlist_t* hl, saidx64_t db_index, exactmatch_t *e);

//hit_t* hitlist_add(hitlist_t* hl, saidx64_t db_index, int gapless_score, int gapless_length, saidx64_t suffix_offset);
hit_t* hitlist_add(hitlist_t* hl, db_t* db, saidx64_t SAindex, int gapless_score, int gapless_length, saidx64_t suffix_offset);

//void hitlist_print(hitlist_t* hl);
traversal_t* hitlist_traversal(hitlist_t* hl);
int hitlist_len(hitlist_t* hl);

#endif

