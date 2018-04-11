#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "seg_wrapper.h"

#include "seg/genwin.h"
#include "seg/seg.h"


void initialise_seg() {
    seg_init();
}

void run_seg(seq_t* s1) {
    struct Sequence *seq = (struct Sequence *) malloc(sizeof(struct Sequence));
    struct Segment *segs = (struct Segment *) NULL;
    struct Segment  *tmp;
    int i;

    // copy sequence
    seq->length = seq_len(s1);
    seq->seq = seq_ptr(s1);
    //seq->seq = (char*) malloc(seq->length);
    //memcpy(seq->seq, seq_ptr(s1), seq->length);

    segseq(seq, &segs, 0);
    mergesegs(seq, segs);

    //fprintf(stderr, "PRE-SEG: %.*s\n", seq->length, seq->seq);

    for (tmp = segs; tmp != NULL; tmp = tmp->next) {
        // HARD MASK
        memset(seq->seq + tmp->begin, 'X', tmp->end - tmp->begin + 1);
        /*
        // SOFT MASK
        for(i = tmp->begin; i < tmp->end - tmp->begin + 1; ++i) {
            seq->seq[i] = seq->seq[i] + 32; // convert to lowercase
        }
        */
    }

    //fprintf(stderr, "POSTSEG: %.*s\n", seq->length, seq->seq);

    freesegs(segs);
    //closeseq(seq);
    free(seq);
}

