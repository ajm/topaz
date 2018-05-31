#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <errno.h>
#include <unistd.h>
#include <sys/time.h>

#include <omp.h>
#include <pthread.h>

#include "divsufsort64.h"

#include "options.h"
#include "seq.h"
#include "fasta.h"
#include "db.h"
#include "search.h"
#include "utils.h"
#include "engine.h"
#include "hit.h"
#include "alignment.h"
#include "ssw_wrapper.h"
#include "ssw.h"
#include "vector.h"
#include "suffix.h"
#include "query.h"
#include "work_queue.h"
#include "seg_wrapper.h"

//#include <gperftools/profiler.h>

//#define NO_OUTPUT 1
#define DEFAULT_QUEUE_LENGTH 100

search_t* search_alloc(char *queries, char *prefix, options_t *opt) {
    search_t *s;
    int error_num;

    s = (search_t *) malloc(sizeof(search_t));
    s->f = fasta_alloc(queries);
    s->opt = opt;

    s->db = db_alloc(prefix);
    if(db_load(s->db) < 0) {
        fprintf(stderr, "Error: could not load database, %s\n", strerror(errno));
        exit(EXIT_FAILURE);
    }

    if(opt->db_size == 0) {
        opt->db_size = s->db->Tsize - s->db->Psize; // - Psize because of delimiters
        opt->db_proteins = s->db->Psize;

        fprintf(stderr, "db size = %" PRIu64 "\n", opt->db_size);
        fprintf(stderr, "db proteins = %" PRIu64 "\n", opt->db_proteins);
    }

    s->timing = fopen("timing.txt", "w");
    fprintf(s->timing, "id\tlength\tattempts\tsuccesses\ttime1\ttime2\ttime3\n");

    //s->q_search   = wqueue_alloc(DEFAULT_QUEUE_LENGTH);
    //s->q_alignment = wqueue_alloc(DEFAULT_QUEUE_LENGTH);
    s->q_search   = wqueue2_alloc();
    s->q_alignment = wqueue2_alloc();

    s->queries_done = 0;
    s->queries_added = 0;
    s->queries_completed = 0;

    s->thread_number = 0;
    gettimeofday(&s->start_time, NULL);
    s->total_suffix_time = 0.0;
    s->total_alignment_time = 0.0;


    if((error_num = pthread_mutex_init(&s->mutex, NULL)) != 0) {
        fprintf(stderr, "Error: pthread_mutex_init failed: %s\n", strerror(error_num));
        _exit(EXIT_FAILURE);
    }

    return s;
}

//void _query_deallocator(void* payload) {
//    query_free((query_t*) payload);
//}

void search_free(search_t *s) {
    if(s) {
        db_free(s->db);
        fasta_free(s->f);
        fclose(s->timing);
        //wqueue_free(s->q_search, _query_deallocator);
        //wqueue_free(s->q_alignment, _query_deallocator);
        wqueue2_free(s->q_search);
        wqueue2_free(s->q_alignment);
        free(s);
    }
}

double elapsed_time(struct timeval* start, struct timeval* stop) {
    double tmp;

    tmp = (stop->tv_sec - start->tv_sec) * 1000.0;
    tmp += (stop->tv_usec - start->tv_usec) / 1000.0;

    return tmp;
}

void _progress(search_t* s) {
    struct timeval current_time;

    pthread_mutex_lock(&s->mutex);
    
    gettimeofday(&current_time, NULL);
/*
    fprintf(stderr, 
        "\r[%d/%d%s] suffix(%d) alignment(%d), %.1f sec      ", 
        s->queries_completed, 
        s->queries_added, 
        s->queries_done ? "" : "*",
        wqueue2_length(s->q_search),
        wqueue2_length(s->q_alignment),
        elapsed_time(&s->start_time, &current_time) / 1000);
*/
    fprintf(stderr, "\rPROGRESS: %d/%d (%.1f sec)", 
            s->queries_completed, 
            s->queries_added, 
            elapsed_time(&s->start_time, &current_time) / 1000);

    pthread_mutex_unlock(&s->mutex);
}

int _compare_bitscore(const void *l, const void *r) {
    double lbs = (*(hit_t**) l)->bitscore;
    double rbs = (*(hit_t**) r)->bitscore;

    if(lbs == rbs) {
        return 0;
    }
    else if(lbs < rbs) {
        return 1;
    }

    return -1;
}

void hit_print(hit_t* h, seq_t* q, seq_t *s, int fastmode) {

    if(!fastmode) {
        printf("%.*s\t%.*s\t%.2f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.3e\t%d\n",
        seq_idlen(q), seq_idptr(q), 
        seq_idlen(s), seq_idptr(s),
        h->identity * 100,
        h->length,
        h->mismatch,
        h->gapopen,
        h->qstart + 1,
        h->qend + 1,
        h->sstart + 1,
        h->send + 1,
        h->evalue,
        (int) h->bitscore);
    }
    else {
        printf("%.*s\t%.*s\t%.3e\t%d\n",
        seq_idlen(q), seq_idptr(q), 
        seq_idlen(s), seq_idptr(s),
        h->evalue,
        (int) h->bitscore);
    }

#if 0
    // qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    printf("%.*s\t%.*s\t%.2f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.3e\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
        seq_idlen(q), seq_idptr(q), 
        seq_idlen(s), seq_idptr(s),
        h->identity * 100,
        h->length,
        h->mismatch,
        h->gapopen,
        h->qstart,
        h->qend,
        h->sstart,
        h->send,
        h->evalue,
        (int) h->bitscore,

        h->alignment_rank,
        h->gapless_score,
        h->gapless_length,
        seq_len(q),
        seq_len(s),
        h->suffix_counter,
        h->suffix_offset,
        h->suffix_multiple
    );
#endif
}

void output_top_hits(hitlist_t* hl, seq_t* query, db_t* db, options_t* opt) {
    traversal_t* t;
    hit_t* h;
    hit_t** all_hits;
    seq_t *tmp, *s;
    int i = 0;
    size_t len = hitlist_len(hl);
    size_t aligned_hits = 0;
    int count = 0;

    all_hits = malloc(sizeof(hit_t*) * len);

    t = hitlist_traversal(hl);
    traversal_start(t);

    while((h = traversal_next(t)) != NULL) {
        if(h->is_aligned) {
            all_hits[aligned_hits++] = h;
        }
    }

    //fprintf(stderr, "[HITLIST] %zu elements, %zu aligned\n", len, aligned_hits);

    qsort(all_hits, aligned_hits, sizeof(hit_t*), _compare_bitscore);

    aligned_hits = aligned_hits < opt->num_hits ? aligned_hits : opt->num_hits;

    s = seq_alloc();

    for(i = 0; i < aligned_hits; ++i) {
        h = all_hits[i];

        if(h->evalue <= opt->evalue) {
            if((tmp = db_seq(db, h->db_index, &s)) == NULL) {
                fprintf(stderr, "Error: protein at index %" PRIu64 " not found!", h->db_index);
                exit(EXIT_FAILURE);
            }

            hit_print(h, query, tmp, opt->superfast);
            ++count;
        }
    }

    seq_free(s);
    traversal_free(t);
    free(all_hits);

#ifndef NO_OUTPUT
    if(count == 0) {
        fprintf(stderr, "Warning: %.*s generated zero hits!\n", seq_idlen(query), seq_idptr(query));
    }
#endif
}

void read_data(search_t* s) {
    seq_t* seq;
    query_t* query;

    if(s->opt->seg_enabled) {
        initialise_seg();
    }

    while((seq = fasta_next(s->f, NULL)) != NULL) {

        // if this is going to happen, it has to be now
        if(s->opt->seg_enabled) {
            run_seg(seq);
        }

#ifndef USELEXICOGRAPHICAL
        seq_2internal(seq);
#endif
        query = query_alloc(seq, s->db, s->opt);
        //wqueue_prepend_blocking(s->q_search, (void*) query);
        wqueue2_push(s->q_search, query);
        s->queries_added++;

        //seq_write(seq, SEQ_ID | SEQ_SEQ, stderr);
    }

    s->queries_done = 1;
}

/*
void *pthread_read_data_worker(void *arg) {
    search_t* s = (search_t*) arg;
    seq_t* seq;
    query_t* query;

    while((seq = fasta_next(s->f, NULL)) != NULL) {
        query = query_alloc(seq, s->db, s->opt);
        //wqueue_prepend_blocking(s->q_search, (void*) query);
        wqueue2_push(s->q_search, query);
        s->queries_added++;
    }

    s->queries_done = 1;

    pthread_exit(NULL);
}
*/

int _incomplete(search_t* s) {
    return !(s->queries_done && (s->queries_added == s->queries_completed));
}

void _increment_complete(search_t* s) {
    pthread_mutex_lock(&s->mutex);
    s->queries_completed++;
    pthread_mutex_unlock(&s->mutex);
}

void _search(search_t* s, query_t* query) {
    if(query->state == QUERY_FIND_SUFFIXES) {
        query_find_suffixes(query);
    }

    query_fetch_alignment_batch(query);

    if(! query_done(query)) {
        //wqueue_prepend(s->q_alignment, (void*) query);
        wqueue2_push(s->q_alignment, query);
    }
    else {
        output_top_hits(query->hl, query->query, query->db, query->opt);
        query_free(query);
        _increment_complete(s);
    }
}

void _alignment(search_t* s, query_t* query) {
    query_perform_alignments(query);

    if(! query_done(query)) {
        //wqueue_prepend(s->q_search, (void*) query);
        wqueue2_push(s->q_search, query);
    }
    else {
        output_top_hits(query->hl, query->query, query->db, query->opt);
        query_free(query);
        _increment_complete(s);
    }
}

void _full_search(search_t* s, query_t* query, double* suffix_time, double* alignment_time) {
    struct timeval start_time, end_time;

    while(query->state != QUERY_DONE) {
        switch(query->state) {
            case QUERY_FIND_SUFFIXES :
                gettimeofday(&start_time, NULL);
                query_find_suffixes(query);
                gettimeofday(&end_time, NULL);
                *suffix_time += elapsed_time(&start_time, &end_time);
                break;
                
            case QUERY_GET_ALIGNMENTS :
                gettimeofday(&start_time, NULL);
                query_fetch_alignment_batch(query);
                gettimeofday(&end_time, NULL);
                *suffix_time += elapsed_time(&start_time, &end_time);
                break;

            case QUERY_DO_ALIGNMENTS :
                gettimeofday(&start_time, NULL);
                query_perform_alignments(query);
                gettimeofday(&end_time, NULL);
                *alignment_time += elapsed_time(&start_time, &end_time);
                break;

            default :
                break;
        }
    }

    output_top_hits(query->hl, query->query, query->db, query->opt);
    query_free(query);
    _increment_complete(s);
}

void *pthread_suffix_worker(void *arg) {
    search_t *s = (search_t*) arg;
    query_t *query;
    int thread_id;

    pthread_mutex_lock(&s->mutex);
    thread_id = s->thread_number++;
    pthread_mutex_unlock(&s->mutex);

    fprintf(stderr, " -> i am thread %d (suffix)\n", thread_id);

    while(_incomplete(s)) {
        if((query = wqueue2_pop(s->q_search, 0)) != NULL) {
            _search(s, query);
        }
        else {
            if(_incomplete(s) && ((query = wqueue2_pop(s->q_alignment, 0)) != NULL)) {
                //fprintf(stderr, "S %d %.*s\n", thread_id, seq_idlen(query->query), seq_idptr(query->query));
                _alignment(s, query);
            }
        }

        //_progress(s);
    }

    fprintf(stderr, " -> thread %d exiting (suffix)\n", thread_id);

    wqueue2_wakeup(s->q_alignment);
    pthread_exit(NULL);
}

void *pthread_generic_worker(void *arg) {
    search_t *s = (search_t*) arg;
    query_t *query;
    //int thread_id;
    double suffix_time, alignment_time;

    suffix_time = 0.0;
    alignment_time = 0.0;
/*
    pthread_mutex_lock(&s->mutex);
    thread_id = s->thread_number++;
    pthread_mutex_unlock(&s->mutex);
*/
    //fprintf(stderr, " -> i am thread %d (generic)\n", thread_id);

    while(_incomplete(s)) {
        if((query = wqueue2_pop(s->q_search, 0)) != NULL) {

//            ProfilerStart("topaz.log");

            _full_search(s, query, &suffix_time, &alignment_time);

//            ProfilerStop();
        }

        //_progress(s);
    }

    pthread_mutex_lock(&s->mutex);
    s->total_suffix_time += suffix_time;
    s->total_alignment_time += alignment_time;
    pthread_mutex_unlock(&s->mutex);

    pthread_exit(NULL);
}

void *pthread_alignment_worker(void *arg) {
    search_t *s = (search_t*) arg;
    query_t *query;
    int thread_id;

    pthread_mutex_lock(&s->mutex);
    thread_id = s->thread_number++;
    pthread_mutex_unlock(&s->mutex);

    fprintf(stderr, " -> i am thread %d (alignment)\n", thread_id);

    while(_incomplete(s)) {
        if((query = wqueue2_pop(s->q_alignment, 1)) != NULL) {
            //fprintf(stderr, "A %d %.*s\n", thread_id, seq_idlen(query->query), seq_idptr(query->query));
            _alignment(s, query);
        }

        //_progress(s);
    }

    fprintf(stderr, " -> thread %d exiting (alignment)\n", thread_id);

    wqueue2_wakeup(s->q_alignment);
    pthread_exit(NULL);
}

void stop_threads(pthread_t* threads, int count) {
    int i;

    for(i = 0; i < count; ++i) {
        pthread_cancel(threads[i]);
    }
}

int search_pthread(search_t *s) {
    pthread_t* threads;
    void *value_ptr;
    int error_number;
    int i;
    struct timeval end_time;

    printf("qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore alignrank gaplessscore gaplesslength qlength slength suffixcounter offset multiple\n");

    threads = (pthread_t*) malloc(sizeof(pthread_t) * s->opt->num_threads);

    fprintf(stderr, "spawning %d threads...\n", s->opt->num_threads);

    for(i = 0; i < s->opt->num_threads; ++i) {
#define _NONGENERIC
#ifdef NONGENERIC
        if(i < 16) {
            if((error_number = pthread_create(&threads[i], NULL, pthread_suffix_worker, s)) != 0) {
                fprintf(stderr, "Error: pthread_create failed: %s\n", strerror(error_number));
                stop_threads(threads, i-1);
                break;
            }
        }
        else {
            if((error_number = pthread_create(&threads[i], NULL, pthread_alignment_worker, s)) != 0) {
                fprintf(stderr, "Error: pthread_create failed: %s\n", strerror(error_number));
                stop_threads(threads, i-1);
                break;
            }
        }
#else
        if((error_number = pthread_create(&threads[i], NULL, pthread_generic_worker, s)) != 0) {
            fprintf(stderr, "Error: pthread_create failed: %s\n", strerror(error_number));
            stop_threads(threads, i-1);
            break;
        }
#endif
        //fprintf(stderr, "created thread-%d\n", i);
    }

    read_data(s);

    //pthread_create(&threads[0], NULL, pthread_read_data_worker, s);
    //pthread_generic_worker((void*) s);

    for(i = 0; i < s->opt->num_threads; ++i) {
        pthread_join(threads[i], &value_ptr);
        //fprintf(stderr, "joined thread-%d\n", i);
    }

    gettimeofday(&end_time, NULL);
/*
    fprintf(stderr, 
            "\nsuffix: %.3f ms\nalignment: %.3f ms\ntotal: %.3f\n", 
            s->total_suffix_time, 
            s->total_alignment_time, 
            elapsed_time(&s->start_time, &end_time));
*/
    fprintf(stderr, "done!\n");

    return 0;
}

