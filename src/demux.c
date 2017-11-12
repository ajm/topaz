#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <ctype.h>
#include <errno.h>
#include <sys/mman.h>

#include "seq.h"
#include "utils.h"
#include "fasta.h"


int demultiplex_fasta(char* fname, char* prefix) {
    FILE *pin, *psq, *phr;
    seq_t* s;
    fasta_t* f;

    pin = fopen_pf(prefix, ".pin");
    psq = fopen_pf(prefix, ".psq");
    phr = fopen_pf(prefix, ".phr");

    s = seq_alloc();
    f = fasta_alloc(fname);

    while((s = fasta_next(f, s)) != NULL) {

        fprintf(pin, "%ld %d %ld %d\n", ftell(psq), s->seq_len, ftell(phr), s->id_len);

#ifndef USELEXICOGRAPHICAL
        seq_2internal(s);
#endif

        seq_write(s, SEQ_ID, phr);
        seq_write(s, SEQ_SEQ, psq);
    }

    fclose(phr);
    fclose(psq);
    fclose(pin);

    fprintf(stderr, "Read %" PRIu64 " sequences from %s\n", f->count, fname);

    fasta_free(f);
    seq_free(s);

    return 0;
}

