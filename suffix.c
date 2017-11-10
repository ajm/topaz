#include <stdlib.h>
#include <string.h>

#include "suffix.h"
#include "db.h"
#include "search.h"


suffix_t* suffix_alloc(seq_t* seq, db_t *db, int offset, saidx64_t position, int direction, int gapless_score, int gapless_length) {
    suffix_t* suf;

    suf = (suffix_t*) malloc(sizeof(suffix_t));
    suf->seq = seq;
    suf->offset = offset;
    suf->position = position;
    suf->direction = direction;
    suf->gapless_score = gapless_score;
    suf->gapless_length = gapless_length;

    return suf;
}

suffix_t* suffix_duplicate(suffix_t* suf) {
    suffix_t* tmp;

    tmp = (suffix_t*) malloc(sizeof(suffix_t));
    
    memcpy(tmp, suf, sizeof(suffix_t));
    tmp->direction *= -1;

    return tmp;
}

void suffix_free(suffix_t *suf) {
    if(suf) {
        free(suf);
    }
}

int suffix_increment(suffix_t *suf, db_t* db, options_t* opt) {
    int lcp_length;
    saidx64_t num_matches;

    suf->position += suf->direction;

    if(!db_valid_index(db, suf->position))
        return -1;

    sa_hsp(db, 
            suf->seq, 
            suf->offset, 
            suf->position, 
            opt, 
            &suf->gapless_score, 
            &suf->gapless_length);

    return 0;
}

int suffix_active(suffix_t *suf, int min_length) {
    return (suf->gapless_length >= min_length);
}

// unfortunately this has to be an integer at the moment...
uint64_t suffix_priority(suffix_t *suf) {
    return suf->gapless_score;
}

