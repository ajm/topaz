#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <unistd.h>
#include <pthread.h>

#include "utils.h"
#include "fasta.h"
#include "seq.h"


fasta_t* fasta_alloc(char *fname) {
    fasta_t *fbuf;
    int error_number;

    fbuf = (fasta_t*) malloc(sizeof(fasta_t));

    fbuf->fname = strdup(fname);
    fbuf->pos = 0;
    fbuf->count = 0;

    if((error_number = pthread_mutex_init(&fbuf->read_mutex, NULL)) != 0) {
        fprintf(stderr, "Error: pthread_mutex_init failed: %s\n", strerror(error_number));
        exit(EXIT_FAILURE);
    }

    if((fbuf->T = (char*) fmmap_ro(fname, &fbuf->size)) == NULL) {
        fprintf(stderr, "Error: problem mapping %s, %s\n", fname, strerror(errno));
        exit(EXIT_FAILURE);
    }

    if(madvise(fbuf->T, fbuf->size, MADV_SEQUENTIAL) < 0) {
        fprintf(stderr, "Error: madvise failed on %s, %s\n", fname, strerror(errno));
        exit(EXIT_FAILURE);
    }

    if(fbuf->T[0] != '>') {
        fprintf(stderr, "Error: %s does not appear to be a FASTA file (first char not '>')\n", fname);
        exit(EXIT_FAILURE);
    }

    return fbuf;
}

void fasta_free(fasta_t *fbuf) {
    if(fbuf) {
        if(munmap(fbuf->T, fbuf->size) < 0) {
            fprintf(stderr, "Error: problem unmapping %s, %s\n", fbuf->fname, strerror(errno));
            exit(EXIT_FAILURE);
        }

        free(fbuf->fname);
        free(fbuf);
    }
}

int _create_seq(seq_t* s, char* data, size_t len) {
    char *pos;
    size_t ind;

    if((pos = memchr(data, '\n', len)) == NULL)
        return -1;

    ind = pos - data;

    if(ind > len)
        return -1;

    seq_id(s, data, ind);
    seq_seq(s, pos + 1, len - ind, 1);

    return 0;
}

seq_t* fasta_next(fasta_t *fbuf, seq_t *sbuf) {
    seq_t *tmp;
    uint64_t i;
    int error_number;

    if(fbuf->pos == fbuf->size)
        return NULL;

    tmp = sbuf ? sbuf : seq_alloc();

    if((error_number = pthread_mutex_lock(&fbuf->read_mutex)) != 0) {
        fprintf(stderr, 
                "Error: pthread_mutex_lock failed: %s\n", 
                strerror(error_number));
        _exit(EXIT_FAILURE);
    }

    for(i = fbuf->pos + 1; i < fbuf->size + 1; ++i) {
        if(i == fbuf->size || (fbuf->T[i-1] == '\n' && fbuf->T[i] == '>')) {
            if(_create_seq(tmp, fbuf->T + fbuf->pos, i - fbuf->pos - 1) < 0) {
                fprintf(stderr, 
                        "Error: encountered a problem parsing FASTA file (%" PRIu64 " sequences were read without error)\n", 
                        fbuf->count);
                exit(EXIT_FAILURE);
            }
            fbuf->pos = i;
            break;
        }
    }

    fbuf->count++;
    
    if((error_number = pthread_mutex_unlock(&fbuf->read_mutex)) != 0) {
        fprintf(stderr, 
                "Error: pthread_mutex_unlock failed: %s\n", 
                strerror(error_number));
        _exit(EXIT_FAILURE);
    }

    return tmp;
}

