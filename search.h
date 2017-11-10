#ifndef H_SASEARCH
#define H_SASEARCH

#include <sys/types.h>

#include "divsufsort64.h"

#include "db.h"
#include "seq.h"
#include "options.h"

int sa_lcp(db_t *db, seq_t *q, int qoffset, saidx64_t position);
int sa_gapless(db_t *db, seq_t *q, int qoffset, saidx64_t position, options_t* opt, int* gapless_score, int* gapless_length);
int sa_hsp(db_t *db, seq_t *q, int qoffset, saidx64_t position, options_t* opt, int* gapless_score, int* gapless_length);
int sa_search(db_t *db, seq_t *q, saidx64_t offset, saidx64_t length, saidx64_t *idx, saidx64_t *num_matches, saidx64_t *last_position);

#endif

