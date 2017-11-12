#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <inttypes.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/types.h>

#include "divsufsort64.h"

#include "utils.h"
#include "db.h"
#include "seq.h"


db_t* db_alloc(char *prefix) {
    db_t *db;

    db = (db_t*) malloc(sizeof(db_t));
    db->prefix = strdup(prefix);

    return db;
}

void db_free(db_t *db) {
    if(db) {
        db_close(db);

        free(db->prefix);
        free(db);
    }
}

int db_load(db_t *db) {
    char *T_fname, 
         *SA_fname,
         *phr_fname;

#ifdef USE_LCP
    char *LCPLC_fname, *LCPCR_fname;
    uint64_t LCPsize;
#endif

    phr_fname = mk_fname(db->prefix, ".phr");
    T_fname = mk_fname(db->prefix, ".psq");
    SA_fname = mk_fname(db->prefix, ".SA");

    if((db->T = (char*) fmmap_ro(T_fname, &db->Tsize)) == NULL) {
        return -1;
    }

    fprintf(stderr, "read T (%s), file size is %" PRIu64 " bytes\n", T_fname, db->Tsize);

    if((db->SA = (saidx64_t*) fmmap_ro(SA_fname, &db->SAsize)) == NULL) {
        return -1;
    }

    fprintf(stderr, "read SA (%s), file size is %" PRIu64 " bytes\n", SA_fname, db->SAsize);


    if(db->SAsize != (db->Tsize * sizeof(saidx64_t))) {
        fprintf(stderr, "Error: SA (%s, %" PRIu64 " bytes) should be %lu x T (%s, %" PRIu64 " bytes)!\nExiting...\n", SA_fname, db->SAsize, sizeof(saidx64_t), T_fname, db->Tsize);
        exit(EXIT_FAILURE);
    }


#ifdef USE_LCP
    LCPLC_fname = mk_fname(db->prefix, ".LCPLC");
    LCPCR_fname = mk_fname(db->prefix, ".LCPCR");

    if((db->LCPLC = (saidx64_t*) fmmap_ro(LCPLC_fname, &LCPsize)) == NULL) {
        return -1;
    }

    fprintf(stderr, "read LCPLC (%s), file size is %" PRIu64 " bytes\n", LCPLC_fname, LCPsize);

    if((db->LCPCR = (saidx64_t*) fmmap_ro(LCPCR_fname, &LCPsize)) == NULL) {
        return -1;
    }

    fprintf(stderr, "read LCPCR (%s), file size is %" PRIu64 " bytes\n", LCPCR_fname, LCPsize);

    free(LCPLC_fname);
    free(LCPCR_fname);
#endif

    if((db->headers = (char*) fmmap_ro(phr_fname, &db->Hsize)) == NULL) {
        return -1;
    }

    fprintf(stderr, "read protein headers (%s), file size is %" PRIu64 " bytes\n", phr_fname, db->Hsize);

    // Psize = number of proteins
    db->Psize = strncount(db->headers, '\n', db->Hsize);
    fprintf(stderr, "%" PRIu64 " proteins in database\n", db->Psize);
    

    if(_load_index(db) != 0) {
        return -1;
    }   

    free(T_fname);
    free(SA_fname);
    free(phr_fname);

    return 0;
}

int _load_index(db_t *db) {
    char *fname;
    uint64_t Isize;

    fname = mk_fname(db->prefix, ".ind");

    if((db->index = (position_t*) fmmap_ro(fname, &Isize)) == NULL) {
        fprintf(stderr, "Error: could not open %s, %s\n", fname, strerror(errno));
        return -1;
    }

    if(Isize != db->Psize * sizeof(position_t)) {
        fprintf(stderr, "Error: index is the wrong size (read %" PRIu64 ", expected %" PRIu64 ")", Isize, db->Psize * sizeof(position_t));
        exit(EXIT_FAILURE);
    }

    return 0;
}

int _load_index_manual(db_t *db) {
    size_t linecap;
    ssize_t linelen;
    uint64_t pstart, hstart, linenum = 0, length;
    int plen, hlen, nmatch;
    char *fname, *line = NULL;
    FILE *f;

    fname = mk_fname(db->prefix, ".pin");
    
    if((f = fopen(fname, "r")) == NULL) {
        return -1;
    }

    db->index = (position_t *) malloc(db->Psize * sizeof(position_t));

    while((linelen = getline(&line, &linecap, f)) > 0) {
        if((nmatch = sscanf(line, "%" PRIu64 " %d %" PRIu64 " %d\n", &pstart, &plen, &hstart, &hlen)) != 4) {
            fprintf(stderr, "Error: line %" PRIu64 " of %s contained %d elements (expected 4)\n", linenum + 1, fname, nmatch);
            exit(EXIT_FAILURE);
        }

        db->index[linenum].header_start = hstart;
        db->index[linenum].protein_start = pstart;
        db->index[linenum].header_len = hlen;
        db->index[linenum].protein_len = plen;

        linenum++;
    }

    if(linenum != db->Psize) {
        fprintf(stderr, "Error: expected %" PRIu64 " index entries, but read %" PRIu64 " from %s\n", db->Psize, linenum, fname);
        exit(EXIT_FAILURE);
    }

    length = db->index[db->Psize-1].protein_start + db->index[db->Psize-1].protein_len;

    if(db->Tsize != (length + 1)) {
        fprintf(stderr, "Error: index length is %" PRIu64 ", but T length is %" PRIu64 " ...\n", length, db->Tsize);
        exit(EXIT_FAILURE);
    }

    fclose(f);
    free(fname);
    free(line);
    return 0;
}

int db_index(db_t *db, uint64_t pos, saidx64_t *index) {
    uint64_t i, half, size;
    position_t *tmp;

    if(pos >= db->Tsize) {
        return -1;
    }

    for(i = 0, size = db->Psize, half = size >> 1; 
        0 < size; 
        size = half, half >>= 1) {

        tmp = db->index + i + half;

        if(pos < tmp->protein_start) {
            continue;
        }
        
        if(pos > (tmp->protein_start + tmp->protein_len)) {
            i += half + 1;
            half -= (size & 1) ^ 1;
            continue;
        }

        break;
    }

    *index = i + half;
    return 0;
}

int db_close(db_t *db) {
    return munmap(db->T, db->Tsize) == 0
        || munmap(db->SA, db->SAsize) == 0
        || munmap(db->headers, db->Hsize) == 0
        || munmap(db->index, db->Psize * sizeof(position_t)) == 0
#ifdef USE_LCP
        || munmap(db->LCPLC, db->SAsize) == 0
        || munmap(db->LCPCR, db->SAsize) == 0
#endif
    ;
}

int db_proteinid(db_t* db, uint64_t i, char** id, size_t* len) {
    position_t* pos;

    if(i >= db->Psize)
        return -1;
    
    pos = db->index + i;

    *id = db->headers + pos->header_start;
    *len = pos->header_len;

    return 0;
}

seq_t* db_seq(db_t* db, uint64_t i, seq_t **s) {
    seq_t* seq;
    position_t* pos;

    if(i >= db->Psize)
        return NULL;

    seq = s != NULL ? *s : seq_alloc();
    pos = db->index + i;

    seq_id( seq, db->headers + pos->header_start,  pos->header_len);
    seq_seq(seq, db->T       + pos->protein_start, pos->protein_len);

    return seq;
}

int db_valid_index(db_t *db, saidx64_t ind) {
    return ind >= 0 && ind < db->SAsize;
}

