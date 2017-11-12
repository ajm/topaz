#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <unistd.h>
#include <math.h>

#include "divsufsort64.h"

#include "query.h"
#include "hit.h"
#include "priority_queue.h"
#include "seq.h"
#include "search.h"
#include "ssw_wrapper.h"
#include "ssw.h"

#ifdef USE_PARASAIL
#include "parasail_wrapper.h"
#include <parasail.h>
#endif

query_t* query_alloc(seq_t* s, db_t* db, options_t* opt) {
    query_t* q;
    double correction, corrected_q, corrected_db;

    q = (query_t*) malloc(sizeof(query_t));

    q->query = s;
    q->db = db;
    q->opt = opt;
    q->pq = pqueue_alloc(db, opt);
    q->hl = hitlist_alloc(opt);
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
    //return (s[0] == 'V'); // THIS MEANS 'X' IN THE INTERNAL REPRESENTATION OF AMINOACIDS
        //|| ((s[0] == s[2]) && ((s[0] == s[1]) || (s[0] == s[4])));
    return (s[0] == 'V') || ((s[0] == s[2]) && ((s[0] == s[1]) || (s[0] == s[4])));
}

int query_find_suffixes(query_t* q) {
    int i, gapless_length, gapless_score, num_suffixes;
    saidx64_t pos, num_mat, last_pos;

    num_suffixes = seq_len(q->query) - q->opt->min_hit_length;

    //fprintf(stderr, "starting %.*s ...\n", seq_idlen(q->query), seq_idptr(q->query));

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

    //fprintf(stderr, "looking for hits to align...\n");
    //fprintf(stderr, "starting %.*s (iteration %d)...\n", seq_idlen(q->query), seq_idptr(q->query), q->iterations);

    while(pqueue_next(q->pq, &SAindex, &Qoffset, &gapless_score, &gapless_length) == 0) {
/*
        if(db_index(q->db, q->db->SA[SAindex], &Sindex) < 0) {
            fprintf(stderr, "Error: internal error in protein lookup\n");
            _exit(EXIT_FAILURE);
        }
*/
        //Soffset = q->db->SA[SAindex] - q->db->index[Sindex].protein_start;

        //if((h = hitlist_add(q->hl, Sindex, gapless_score, gapless_length, Qoffset)) != NULL) {
        if((h = hitlist_add(q->hl, q->db, SAindex, gapless_score, gapless_length, Qoffset)) != NULL) {
            // the gapless score alone is sufficient to get a significant E-value
//            if(gapless_score >= q->min_alignment_score) {
//                h->is_ready = 1;
//            }

            //debug
            h->is_ready = 1;

            //fprintf(stderr, "  added %d\n", counter);

            if(h->is_ready && (++counter == batch_size)) {
                break;
            }
        }
        //else {
        //    fprintf(stderr, "nope\n");
        //}
    }

    //fprintf(stderr, "got %d possible hits to align...\n", hitlist_len(q->hl));

    // QUERY_DONE here means that we should abandon traversing the SA
    // because the min. length and assorted criteria were exceeded
    q->state = (counter == 0) ? QUERY_DONE : QUERY_DO_ALIGNMENTS;

    return 0;
}

int query_perform_alignments(query_t* q) {
    traversal_t* tr;
    hit_t* h;

    //fprintf(stderr, "doing alignments (%d)...\n", hitlist_len(q->hl));

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

        if(h->is_ready) {
            h->alignment_rank = q->alignment_counter++;
        }

#ifdef USE_PARASAIL
        if(h->is_ready && perform_alignment_with_profile(h, q->ssw_profile, q->parasail_profile, q->db, q->opt)) {
#else
        if(h->is_ready && perform_alignment_with_profile(h, q->ssw_profile, q->db, q->opt)) {
#endif
            q->hit_counter++;
        }
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

    //fprintf(stderr, "alignments done (%d)\n", q->hit_counter);

    q->state = query_stop_looking(q) ? QUERY_DONE : QUERY_GET_ALIGNMENTS;

    return 0;
}

#ifdef USE_PARASAIL
int perform_alignment_with_profile(hit_t *h, s_profile *profile, parasail_profile_t* parasail_profile, db_t *db, options_t *opt) {
#else
int perform_alignment_with_profile(hit_t *h, s_profile *profile, db_t *db, options_t *opt) {
#endif
    seq_t *s, *tmp;

    s = seq_alloc();

    if((tmp = db_seq(db, h->db_index, &s)) == NULL) {
        fprintf(stderr, "Error: protein at index %" PRIu64 " not found!", h->db_index);
        _exit(EXIT_FAILURE);
    }

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
    // ASAP
    //return (q->hit_counter >= q->opt->num_hits) || ((q->opt->max_alignments != 0) && (q->alignment_counter >= q->opt->max_alignments));
    //return (q->hit_counter >= 250) || (q->iterations++ >= 50);
    return (q->hit_counter >= (q->opt->num_hits / 4)) || (q->iterations++ >= 50);
    //return 1;

    // use up max_alignments
    //return (q->opt->max_alignments != 0) && (q->alignment_counter >= q->opt->max_alignments);
/*
    // be cleverer, or some approximation thereof...
    traversal_t* tr;
    hit_t* h;
    double* all_hits;
    int index = 0;
    double mean = 0;
    double median = 0;

    // if we have not got enough hits, just check whether we have exceeded max. alignments
    if(q->hit_counter < q->opt->num_hits) {
        return (q->opt->max_alignments != 0) && (q->alignment_counter >= q->opt->max_alignments);
    }

    // fail early if we already have enough hits and 
    // this batch did not improve the kth bitscore in the hitlist
    if(highest_bitscore < q->lowest_bitscore) {
        return 1;
    }

    // kth bitscore is improved, so update the lowest_bitscore field in the query_t object
    tr = hitlist_traversal(q->hl);
    all_hits = malloc(hitlist_len(q->hl) * sizeof(double));

    while((h = (hit_t*) traversal_next(tr)) != NULL) {
        all_hits[index++] = h->bitscore;
        mean += h->bitscore;
    }

    qsort(all_hits, index, sizeof(double), _compare_double);

    q->lowest_bitscore = all_hits[q->opt->num_hits];

    mean /= index;
    median = all_hits[index / 2];

    q->lowest_bitscore = median;

    traversal_free(tr);
    free(all_hits);

    return 0;
*/
}

