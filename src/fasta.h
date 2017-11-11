#ifndef H_FASTA
#define H_FASTA

#include <pthread.h>
#include <stdint.h>
#include "seq.h"

typedef struct {
    char *T;
    char *fname;
    uint64_t size;
    uint64_t pos;
    uint64_t count;

    pthread_mutex_t read_mutex;

} fasta_t;

fasta_t* fasta_alloc(char *fname);
void fasta_free(fasta_t *fbuf);
seq_t* fasta_next(fasta_t *fbuf, seq_t *sbuf);


#endif

