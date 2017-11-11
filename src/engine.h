#ifndef H_ENGINE
#define H_ENGINE

#include <stdio.h>
#include <sys/time.h>
#include <pthread.h>

#include "fasta.h"
#include "db.h"
#include "options.h"
#include "work_queue.h"


typedef struct {
    db_t *db;
    fasta_t *f;
    options_t *opt;
    FILE *timing;

    int queries_done;
    int queries_added;
    int queries_completed;

    pthread_mutex_t mutex;

    //wqueue_t* q_search;
    //wqueue_t* q_alignment;
    wqueue2_t* q_search;
    wqueue2_t* q_alignment;

    int thread_number;
    struct timeval start_time;
    double total_suffix_time;
    double total_alignment_time;

} search_t;

search_t* search_alloc(char *queries, char *prefix, options_t *opt);
void search_free(search_t *s);

int search_openmp(search_t *s);
int search_pthread(search_t *s);


#endif

