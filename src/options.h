#ifndef H_OPTIONS
#define H_OPTIONS

#include <stdint.h>

#define DEFAULT_THREADS 1
#define DEFAULT_HITS 100
#define DEFAULT_EVALUE 1e-3
#define DEFAULT_GAPOPEN 11
#define DEFAULT_GAPEXTEND 1
#define DEFAULT_SUBSTITUTION "BLOSUM62"
#define DEFAULT_LAMBDA 0.267
#define DEFAULT_K 0.041
#define DEFAULT_H 0.140
#define DEFAULT_DBSIZE 0
#define DEFAULT_DBPROTEINS 0

#define DEFAULT_MINLENGTH 5
#define DEFAULT_MINSCORE 1
#define DEFAULT_SEEDS 300000
#define DEFAULT_ALIGNMENTS 5000
#define FAST_SEEDS 100000
#define FAST_ALIGNMENTS 1500


typedef struct {
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

    int min_hit_length;
    int min_hit_score;
    int number_of_alignments;
    int number_of_seeds;

    int seg_enabled;
    int superfast;

} options_t;

#endif

