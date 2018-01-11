#ifndef H_SEG
#define H_SEG

struct Segment
  {
   int begin;
   int end;
   struct Segment *next;
  };

void segseq(struct Sequence *seq, struct Segment **segs, int offset);
void mergesegs(struct Sequence *seq, struct Segment *segs);
void freesegs(struct Segment *segs);


#endif

