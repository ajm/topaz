#ifndef H_DB
#define H_DB

#include <stdint.h>

#include "divsufsort64.h"
#include "seq.h"


typedef struct {
    uint64_t header_start;
    uint64_t protein_start;

    int header_len;
    int protein_len;

} position_t;

typedef struct {
    char *prefix;       // db files have same prefix

    char *headers;      // protein headers
    uint64_t Hsize;     // header text length, bytes

    position_t *index;
    uint64_t Psize;     // number of proteins

    char *T;            // protein text, newline delimited
    uint64_t Tsize;     // protein text length, bytes

    saidx64_t *SA;      // suffix array
    uint64_t SAsize;    // suffix array length, bytes

#ifdef USE_LCP
    saidx64_t *LCPLC;   // LCP arrays
    saidx64_t *LCPCR;
#endif

} db_t;


db_t* db_alloc(char *prefix);
void db_free(db_t *db);
int db_load(db_t *db);
int db_close(db_t *db);
int _load_index(db_t *db);
int db_index(db_t *db, uint64_t pos, saidx64_t *index);

int db_proteinid(db_t* db, uint64_t i, char** id, size_t* len);
seq_t* db_seq(db_t* db, uint64_t i, seq_t** s);
int db_valid_index(db_t *db, saidx64_t ind);

#endif

