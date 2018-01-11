#ifndef H_OPTIONS
#define H_OPTIONS

#include <stdint.h>

#define DEFAULT_THREADS 1
//#define DEFAULT_WINDOW 10
#define DEFAULT_HITS 100
//#define DEFAULT_VOTES 4
//#define DEFAULT_VOTELIST DEFAULT_HITS
#define DEFAULT_SMITHWATERMAN 0
#define DEFAULT_EVALUE 1.0
#define DEFAULT_GAPOPEN 11
#define DEFAULT_GAPEXTEND 1
#define DEFAULT_SUBSTITUTION "BLOSUM62"
#define DEFAULT_LAMBDA 0.267
#define DEFAULT_K 0.041
#define DEFAULT_H 0.140
#define DEFAULT_DBSIZE 0
#define DEFAULT_DBPROTEINS 0

#define DEFAULT_MINLENGTH 3
#define DEFAULT_MINSCORE 1
#define DEFAULT_LOOKAHEADS 200000
#define DEFAULT_ALIGNMENTS 5000
//#define DEFAULT_SCORETHRESHOLD 30
// liisa suggested 10-40
//#define DEFAULT_MINDISTANCE 0
//#define DEFAULT_MAXDISTANCE 50
//#define DEFAULT_ALIGN_NOW 12


typedef struct {
//    int suffix_window;
//    int min_votes;
//    int votelist_len;
    int num_threads;

    int num_hits;
    double evalue;

    int smithwaterman;
    int gap_open;
    int gap_extend;

    char *sub_name;
    int8_t (*substitution)[625];
    uint64_t db_size;
    uint64_t db_proteins;

    double lambda;
    double k;
    double H;

//    int align_now_threshold;
    int min_hit_length;
    int min_hit_score;
    int number_of_alignments;
    int number_of_lookaheads;
//    int min_score_threshold;
//    int min_distance;
//    int max_distance;

    int seg_enabled;

    int superfast;

} options_t;

#endif

