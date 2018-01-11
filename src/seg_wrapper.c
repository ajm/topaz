#include <stdlib.h>
#include <string.h>

#include "seg_wrapper.h"

#include "seg/genwin.h"
#include "seg/seg.h"


void initialise_seg() {
    genwininit();
    entropy_init(12); // window length from seg/seg.c
}

void run_seg(seq_t* s1) {
    struct Sequence *seq = (struct Sequence *) malloc(sizeof(struct Sequence));
    struct Segment *segs = (struct Segment *) NULL;
    struct Segment  *tmp;

    // copy sequence
    seq->length = seq_len(s1);
    seq->seq = (char*) malloc(seq->length);
    strncpy(seq_ptr(s1), seq->seq, seq->length);

    segseq(seq, &segs, 0);
    mergesegs(seq, segs);

    fprintf(stderr, "PRE-SEG: %.*s\n", seq_len(s1), seq_ptr(s1));

    // HARD MASK
    for (tmp = segs; tmp != NULL; tmp = tmp->next) {
        memset(seq_ptr(s1) + tmp->begin, 'X', tmp->end - tmp->begin + 1);
    }

    fprintf(stderr, "POSTSEG: %.*s\n", seq_len(s1), seq_ptr(s1));

    freesegs(segs);
    closeseq(seq);
}

