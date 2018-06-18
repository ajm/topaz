#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <unistd.h>
#include <math.h>
#include <sys/time.h>

#include "divsufsort64.h"

#include "query.h"
#include "hit.h"
#include "priority_queue.h"
#include "seq.h"
#include "search.h"
#include "ssw_wrapper.h"
#include "ssw.h"
#include "utils.h"
#include "lexicographical.h"

#ifdef USE_PARASAIL
#include "parasail_wrapper.h"
#include <parasail.h>
#endif

double get_elapsed_time(struct timeval* start, struct timeval* stop) {
    double tmp;

    tmp = (stop->tv_sec - start->tv_sec) * 1000.0;
    tmp += (stop->tv_usec - start->tv_usec) / 1000.0;

    return tmp;
}

query_t* query_alloc(seq_t* s, db_t* db, options_t* opt) {
    query_t* q;
    double correction, corrected_q, corrected_db;

    q = (query_t*) malloc(sizeof(query_t));

    q->query = s;
    q->db = db;
    q->opt = opt;
    q->pq = pqueue_alloc(db, opt);
    q->hl = hitlist_alloc(opt, seq_len(s));
    q->ssw_profile = NULL;
    q->hit_counter = 0;
    q->alignment_counter = 0;
    q->state = QUERY_FIND_SUFFIXES;
    q->iterations = 0;

#ifdef USE_PARASAIL
    q->parasail_profile = NULL;
    q->use_parasail = 1;
#endif

    correction = log(q->opt->k * seq_len(q->query) * q->opt->db_size) / q->opt->H;
    corrected_q = fmax(seq_len(q->query) - correction, 1.0); // / q->opt->k);
    corrected_db = fmax(q->opt->db_size - (correction * q->opt->db_proteins), 1.0);
    q->min_alignment_score = (int) floor((-log(q->opt->evalue / corrected_q / corrected_db) + log(q->opt->k)) / q->opt->lambda);

    //fprintf(stderr, "%.*s %d\n", seq_idlen(q->query), seq_idptr(q->query), q->min_alignment_score);

    return q;
}

void query_free(query_t* q) {
    if(q) {
        seq_free(q->query);
        pqueue_free(q->pq);
        hitlist_free(q->hl);

        if(q->ssw_profile)
            init_destroy2(q->ssw_profile);

#ifdef USE_PARASAIL
        if(q->parasail_profile)
            parasail_profile_free(q->parasail_profile);
#endif

        free(q);
    }
}

int query_done(query_t* q) {
    return q->state == QUERY_DONE;
}

int _bad_prefix(char *s) {
    return (upper_aa_only(s[0]) == 'V'); // THIS MEANS 'X' IN THE INTERNAL REPRESENTATION OF AMINOACIDS
        //|| ((s[0] == s[2]) && ((s[0] == s[1]) || (s[0] == s[4])));
    //return (s[0] == 'V') || ((s[0] == s[2]) && ((s[0] == s[1]) || (s[0] == s[4])));
}

int query_find_suffixes(query_t* q) {
    int i, gapless_length, gapless_score, num_suffixes;
    saidx64_t pos, num_mat, last_pos;

    num_suffixes = seq_len(q->query) - q->opt->min_hit_length;

#ifdef CHATTY
    fprintf(stderr, "starting %.*s ...\n", seq_idlen(q->query), seq_idptr(q->query));
#endif

    for(i = 0; i < num_suffixes; ++i) {
        pos = -1;

        if(_bad_prefix(seq_ptr(q->query) + i))
            continue;

        if(sa_search(q->db, q->query, i, 0, &pos, &num_mat, &last_pos) != 0) {
            fprintf(stderr, "Error: internal error in suffix array search\n");
            _exit(EXIT_FAILURE);
        }

        if(pos != -1) {

            sa_hsp(q->db, q->query, i, pos, q->opt, &gapless_score, &gapless_length);

            if(gapless_length >= q->opt->min_hit_length) {
                pqueue_add(q->pq, suffix_alloc(q->query, q->db, i, pos, 1, gapless_score, gapless_length));
            }

            if(pos != 0) {
                    
                sa_hsp(q->db, q->query, i, pos-1, q->opt, &gapless_score, &gapless_length);

                if(gapless_length >= q->opt->min_hit_length) {
                    pqueue_add(q->pq, suffix_alloc(q->query, q->db, i, pos-1, -1, gapless_score, gapless_length));
                }
            }
        }
    }

    q->state = QUERY_GET_ALIGNMENTS;

    return 0;
}

int query_fetch_alignment_batch(query_t* q) {
    saidx64_t SAindex, Qoffset;
    int gapless_score, gapless_length;
    hit_t *h;
    int batch_size = q->opt->number_of_lookaheads;
    int counter = 0;

#ifdef CHATTY
    fprintf(stderr, "looking for hits to align...\n");
    fprintf(stderr, "starting %.*s (iteration %d)...\n", seq_idlen(q->query), seq_idptr(q->query), q->iterations);
#endif

    while(pqueue_next(q->pq, &SAindex, &Qoffset, &gapless_score, &gapless_length) == 0) {
        if((h = hitlist_add(q->hl, q->db, SAindex, gapless_score, gapless_length, Qoffset)) != NULL) {
            if(++counter == batch_size) {
                break;
            }
        }
    }

#ifdef CHATTY
    fprintf(stderr, "got %d possible hits to align...\n", hitlist_len(q->hl));
#endif

    // QUERY_DONE here means that we should abandon traversing the SA
    // because the min. length and assorted criteria were exceeded
    if(counter == 0) {
        q->state = q->opt->superfast ? QUERY_DONE : QUERY_DO_STATISTICS;
    }
    else {
        q->state = QUERY_DO_ALIGNMENTS;
    }
    //q->state = (counter == 0) ? QUERY_DONE : QUERY_DO_ALIGNMENTS;

    return 0;
}

double alignment_time;

int query_perform_alignments(query_t* q) {
    traversal_t* tr;
    hit_t* h;
    
    alignment_time = 0.0;

#ifdef CHATTY
    fprintf(stderr, "doing alignments (%d)...\n", hitlist_len(q->hl));
#endif

#ifdef USE_PARASAIL
    if(q->use_parasail && q->parasail_profile == NULL) {
        q->parasail_profile = get_parasail_profile(q->query, q->opt);
    }
    if(!q->use_parasail && q->ssw_profile == NULL) {
#else
    if(q->ssw_profile == NULL) {
#endif
        q->ssw_profile = get_ssw_profile(q->query, q->opt);
    }

    tr = hitlist_traversal(q->hl);
    traversal_start(tr);

    while((h = (hit_t*) traversal_next(tr)) != NULL) {
        if(h->is_aligned) {
            continue;
        }

        h->alignment_rank = q->alignment_counter++;

        //gettimeofday(&start_time, NULL);

#ifdef USE_PARASAIL
        if(perform_alignment_with_profile(h, q->ssw_profile, q->parasail_profile, q->db, q->opt)) {
#else
        if(perform_alignment_with_profile(h, q->ssw_profile, q->db, q->opt)) {
#endif
            q->hit_counter++;
        }

        //gettimeofday(&end_time, NULL);
        //alignment_time += get_elapsed_time(&start_time, &end_time);
/*
        // only align top 1000
        if(q->hit_counter >= 1000) {
            fprintf(stderr, "done %d alignments, stopping...\n", q->hit_counter);
            break;
        }
*/
    }

    traversal_end(tr);
    traversal_free(tr);

#ifdef CHATTY
    fprintf(stderr, "alignments done (%d)\n", q->hit_counter);
    fprintf(stderr, "alignment time = %f seconds\n", alignment_time / 1000.0);
#endif

    if(query_stop_looking(q)) {
        q->state = q->opt->superfast ? QUERY_DONE : QUERY_DO_STATISTICS;
    }
    else {
        q->state = QUERY_GET_ALIGNMENTS;
    }

    //q->state = query_stop_looking(q) ? QUERY_DONE : QUERY_GET_ALIGNMENTS;

    return 0;
}

int query_perform_statistics(query_t* q) {
    traversal_t* t;
    hit_t* h;
    hit_t** all_hits;
    int i = 0;
    size_t len = hitlist_len(q->hl);
    size_t aligned_hits = 0;
    seq_t *tmp, *s;
    struct timeval start_time, end_time;

    all_hits = malloc(sizeof(hit_t*) * len);
    s = seq_alloc();

    t = hitlist_traversal(q->hl);
    traversal_start(t);

    while((h = traversal_next(t)) != NULL) {
        if(h->is_aligned) {
            all_hits[aligned_hits++] = h;
        }
    }

    qsort(all_hits, aligned_hits, sizeof(hit_t*), compare_bitscore);

    aligned_hits = aligned_hits < q->opt->num_hits ? aligned_hits : q->opt->num_hits;

    for(i = 0; i < aligned_hits; ++i) {
        h = all_hits[i];

        if(h->evalue <= q->opt->evalue) {
            if((tmp = db_seq(q->db, h->db_index, &s)) == NULL) {
                fprintf(stderr, "Error: protein at index %" PRIu64 " not found!", h->db_index);
                exit(EXIT_FAILURE);
            }
            
            gettimeofday(&start_time, NULL);
            stats_using_ssw_profile(q->ssw_profile, tmp, h, q->opt);
            gettimeofday(&end_time, NULL);

            h->time_stats = get_elapsed_time(&start_time, &end_time);
        }
    }

    seq_free(s);
    traversal_free(t);
    free(all_hits);

    q->state = QUERY_DONE;
    return 0;
}

#ifdef USE_PARASAIL
int perform_alignment_with_profile(hit_t *h, s_profile *profile, parasail_profile_t* parasail_profile, db_t *db, options_t *opt) {
#else
int perform_alignment_with_profile(hit_t *h, s_profile *profile, db_t *db, options_t *opt) {
#endif
    seq_t *s, *tmp;
    struct timeval start_time, end_time;

    s = seq_alloc();

    if((tmp = db_seq(db, h->db_index, &s)) == NULL) {
        fprintf(stderr, "Error: protein at index %" PRIu64 " not found!", h->db_index);
        _exit(EXIT_FAILURE);
    }

    gettimeofday(&start_time, NULL);

#ifdef USE_PARASAIL
    if(parasail_profile) {
        align_using_parasail_profile(parasail_profile, tmp, h, opt);
    }
    else if(profile) {
#endif
        align_using_ssw_profile(profile, tmp, h, opt);
#ifdef USE_PARASAIL
    }
    else {
        fprintf(stderr, "Error: no profile object found\nExiting...\n");
        abort();
    }
#endif

    gettimeofday(&end_time, NULL);
    //alignment_time += get_elapsed_time(&start_time, &end_time);
    h->time_align = get_elapsed_time(&start_time, &end_time);;

    seq_free(s);

    return h->evalue < opt->evalue;
}

int _compare_double(const void* d1, const void* d2) {
    if(d1 == d2) {
        return 0;
    }
    else if(d1 < d2) {
        return 1;
    }
    return -1;
}

int query_stop_looking(query_t *q) {
    return 1;
    return (q->hit_counter >= (q->opt->num_hits / 4)) || (q->iterations++ >= 10);
}

