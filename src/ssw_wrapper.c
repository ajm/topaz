#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>

#include "seq.h"
#include "hit.h"
#include "options.h"
#include "substitution.h"
#include "lexicographical.h"

#include "ssw.h"


//	Print the BLAST like output.
static void ssw_write (const s_align* a,
			const char* ref_seq,
			const char* read_seq,
			const int8_t* table) {

	fprintf(stdout, "optimal_alignment_score: %d\tsub-optimal_alignment_score: %d\t", a->score1, a->score2);
	if (a->ref_begin1 + 1) fprintf(stdout, "target_begin: %d\t", a->ref_begin1 + 1);
	fprintf(stdout, "target_end: %d\t", a->ref_end1 + 1);
	if (a->read_begin1 + 1) fprintf(stdout, "query_begin: %d\t", a->read_begin1 + 1);
	fprintf(stdout, "query_end: %d\n\n", a->read_end1 + 1);
	if (a->cigar) {
		int32_t c = 0, left = 0, e = 0, qb = a->ref_begin1, pb = a->read_begin1;
		uint32_t i;
		while (e < a->cigarLen || left > 0) {
			int32_t count = 0;
			int32_t q = qb;
			int32_t p = pb;
			fprintf(stdout, "Target: %8d    ", q + 1);
			for (c = e; c < a->cigarLen; ++c) {
				char letter = cigar_int_to_op(a->cigar[c]);
				uint32_t length = cigar_int_to_len(a->cigar[c]);
				uint32_t l = (count == 0 && left > 0) ? left: length;
				for (i = 0; i < l; ++i) {
					if (letter == 'I') fprintf(stdout, "-");
					else {
						fprintf(stdout, "%c", *(ref_seq + q));
						++ q;
					}
					++ count;
					if (count == 60) goto step2;
				}
			}
step2:
			fprintf(stdout, "    %d\n                    ", q);
			q = qb;
			count = 0;
			for (c = e; c < a->cigarLen; ++c) {
				char letter = cigar_int_to_op(a->cigar[c]);
				uint32_t length = cigar_int_to_len(a->cigar[c]);
				uint32_t l = (count == 0 && left > 0) ? left: length;
				for (i = 0; i < l; ++i){
					if (letter == 'M') {
						if (table[(int)*(ref_seq + q)] == table[(int)*(read_seq + p)])fprintf(stdout, "|");
						else fprintf(stdout, "*");
						++q;
						++p;
					} else {
						fprintf(stdout, "*");
						if (letter == 'I') ++p;
						else ++q;
					}
					++ count;
					if (count == 60) {
						qb = q;
						goto step3;
					}
				}
			}
step3:
			p = pb;
			fprintf(stdout, "\nQuery:  %8d    ", p + 1);
			count = 0;
			for (c = e; c < a->cigarLen; ++c) {
				char letter = cigar_int_to_op(a->cigar[c]);
				uint32_t length = cigar_int_to_len(a->cigar[c]);
				uint32_t l = (count == 0 && left > 0) ? left: length;
				for (i = 0; i < l; ++i) {
					if (letter == 'D') fprintf(stdout, "-");
					else {
						fprintf(stdout, "%c", *(read_seq + p));
						++p;
					}
					++ count;
					if (count == 60) {
						pb = p;
						left = l - i - 1;
						e = (left == 0) ? (c + 1) : c;
						goto end;
					}
				}
			}
			e = c;
			left = 0;
end:
			fprintf(stdout, "    %d\n\n", p);
		}
	}
}

s_profile* get_ssw_profile(seq_t* s1, options_t* opt) {
    s_profile* profile;
    int8_t *num;
    int i;

    num = (int8_t*) malloc(seq_len(s1));

    //seq_write(s1, SEQ_ID | SEQ_SEQ, stderr);
    //fprintf(stderr, "\n");

    for(i = 0; i < seq_len(s1); ++i) {
#ifdef USELEXICOGRAPHICAL
        num[i] = aa_table[(int) upper_aa(seq_char(s1, i))];
#else
        num[i] = aa_table2[(int) upper_aa(seq_char(s1, i))];
#endif
    }

    profile = ssw_init(num, 
                       seq_len(s1), 
                       opt->substitution, 
                       25, 
                       2);

    fprintf(stderr, "%.*s\n", seq_idlen(s1), seq_idptr(s1));

    return profile;
}

void align_using_ssw_profile(s_profile* profile, seq_t* s2, hit_t* hit, options_t* opt) {
//    s_profile* profile;
    s_align* result;
    double correction, corrected_q, corrected_db;
//    int8_t *num, *ref;
    int8_t *ref;
    int i; //, c, identical;
    //int ref_offset, read_offset;
    //uint32_t cigar_int, cigar_len;
    //char cigar_char;
    int min_alignment_score;

    //seq_write(s2, SEQ_ID | SEQ_SEQ, stderr);
    //fprintf(stderr, "\n");


    fprintf(stderr, "    %.*s\n", seq_idlen(s2), seq_idptr(s2));

    correction = log(opt->k * profile->readLen * opt->db_size) / opt->H;
    corrected_q = fmax(profile->readLen - correction, 1.0); // / opt->k);
    corrected_db = fmax(opt->db_size - (correction * opt->db_proteins), 1.0);

    min_alignment_score = (int) floor((-log(opt->evalue / corrected_q / corrected_db) + log(opt->k)) / opt->lambda);

//    num = (int8_t*) malloc(seq_len(s1));
    ref = (int8_t*) malloc(seq_len(s2));

//    for(i = 0; i < seq_len(s1); ++i)
//        num[i] = aa_table[(int)seq_char(s1, i)]; //seq_char(s1, i) - 'A';

    for(i = 0; i < seq_len(s2); ++i) {
#ifdef USELEXICOGRAPHICAL
        //ref[i] = aa_table[(int) upper_aa(seq_char(s2, i))]; //seq_char(s2, i) - 'A';
        ref[i] = aa_table[(int) seq_char(s2, i)];
#else
        //ref[i] = aa_table2[(int) upper_aa(seq_char(s2, i))]; // XXX
        ref[i] = aa_table2[(int) seq_char(s2, i)];
#endif
    }

//    profile = ssw_init(num, 
//                       seq_len(s1), 
//                       opt->substitution, 
//                       25, 
//                       2);

    result = ssw_align(profile, 
                       ref, 
                       seq_len(s2), 
                       opt->gap_open, 
                       opt->gap_extend, 
                       0, //opt->superfast ? 0 : 2,
                       min_alignment_score, 
                       0, 
                       15);
    
    //ssw_write(result, seq_ptr(s1), seq_ptr(s2), aa_table);

    //printf("alignment done\n");


    hit->sstart = result->ref_begin1;
    hit->send = result->ref_end1;
    hit->qstart = result->read_begin1;
    hit->qend = result->read_end1;
    hit->rawscore = result->score1;

    hit->length = 0;
    hit->gapopen = 0;
    hit->mismatch = 0;
    hit->identity = 0.0;

    hit->bitscore = (opt->lambda * hit->rawscore - log(opt->k)); // / M_LN2;
    hit->evalue = corrected_q * corrected_db * exp(-hit->bitscore);
    hit->evalue = hit->evalue > 1.0e-180 ? hit->evalue : 0.0;
    hit->bitscore /= M_LN2; // nats --> bits

    hit->is_aligned = 1;
    hit->effective_length = profile->readLen - correction;


    ssw_align_stats (profile, ref, seq_len(s2), 
                    opt->gap_open, opt->gap_extend, 15, 
                    hit->qend, hit->send, hit->rawscore, 
                    &(hit->qstart), &(hit->sstart), &(hit->length), 
                    &(hit->gapopen), &(hit->mismatch), &(hit->identity));

    printf("%.*s %d %f %f %d-%d %d-%d\n", seq_idlen(s2), seq_idptr(s2), hit->rawscore, hit->bitscore, hit->evalue, hit->qstart, hit->qend, hit->sstart, hit->send);

    //hit->q_length = seq_len(s2);
    //hit->s_length = profile->readLen;
/*
    identical = 0;
    ref_offset = result->ref_begin1;
    read_offset = result->read_begin1;

    for (c = 0; c < result->cigarLen; ++c) {
        cigar_int = result->cigar[c];
        cigar_len = (uint32_t)(cigar_int >> BAM_CIGAR_SHIFT);
        cigar_char = (cigar_int & 0xfU) > 8 ? 'M': MAPSTR[cigar_int & 0xfU];

        hit->length += cigar_len;
        
        if(cigar_char == 'I') {
            hit->gapopen++;
            read_offset += cigar_len;
        }
        else if(cigar_char == 'D') {
            hit->gapopen++;
            ref_offset += cigar_len;
        }
        else {
            for(i = 0; i < cigar_len; ++i) {
                if((ref_offset < seq_len(s2)) && (read_offset < profile->readLen)) {
                    if(ref[ref_offset] == profile->read[read_offset]) {
                        identical++;
                    }
                    else {
                        hit->mismatch++;
                    }
                    ref_offset++;
                    read_offset++;
                }
            }
        }
    }

    hit->identity = identical / (float)hit->length;
*/
/*
    hit->rawscore = result->score1;
    hit->bitscore = (opt->lambda * hit->rawscore - log(opt->k)); // / M_LN2;
    hit->evalue = corrected_q * corrected_db * exp(-hit->bitscore);
    hit->evalue = hit->evalue > 1.0e-180 ? hit->evalue : 0.0;
    hit->bitscore /= M_LN2; // nats --> bits
*/
/*
    if(correction > hit->length) {
        hit->evalue = INFINITY;
    }
*/
/*
    hit->is_aligned = 1;
    hit->effective_length = profile->readLen - correction;
*/
//    free(num);
    free(ref);
    align_destroy(result);
//    init_destroy(profile);
//
    //return hit->rawscore > min_alignment_score;
}

void align_using_ssw(seq_t* s1, seq_t* s2, hit_t* hit, options_t* opt) {
    s_profile* profile;
    s_align* result;
    double correction, corrected_q, corrected_db;
    int8_t *num, *ref;
    int i, c, identical;
    int ref_offset, read_offset;
    uint32_t cigar_int, cigar_len;
    char cigar_char;


    correction = log(opt->k * seq_len(s1) * opt->db_size) / opt->H;
    corrected_q = fmax(seq_len(s1) - correction, 1.0); // / opt->k);
    corrected_db = fmax(opt->db_size - (correction * opt->db_proteins), 1.0);


    num = (int8_t*) malloc(seq_len(s1));
    ref = (int8_t*) malloc(seq_len(s2));

    for(i = 0; i < seq_len(s1); ++i) {
#ifdef USELEXICOGRAPHICAL
        num[i] = aa_table[(int) upper_aa(seq_char(s1, i))]; //seq_char(s1, i) - 'A';
#else
        num[i] = aa_table2[(int) upper_aa(seq_char(s1, i))]; // XXX
#endif
    }

    for(i = 0; i < seq_len(s2); ++i) {
#ifdef USELEXICOGRAPHICAL
        ref[i] = aa_table[(int) upper_aa(seq_char(s2, i))]; //seq_char(s2, i) - 'A';
#else
        ref[i] = aa_table2[(int) upper_aa(seq_char(s2, i))]; // XXX
#endif
    }

    profile = ssw_init(num, 
                       seq_len(s1), 
                       opt->substitution, 
                       25, 
                       2);

    result = ssw_align(profile, 
                       ref, 
                       seq_len(s2), 
                       opt->gap_open, 
                       opt->gap_extend, 
                       1, 
                       0, 
                       0, 
                       15);
    
    //ssw_write(result, seq_ptr(s1), seq_ptr(s2), aa_table);


    hit->sstart = result->ref_begin1 + 1;
    hit->send = result->ref_end1 + 1;
    hit->qstart = result->read_begin1 + 1;
    hit->qend = result->read_end1 + 1;

    hit->length = 0;
    hit->gapopen = 0;
    hit->mismatch = 0;
    hit->identity = 0.0;

    identical = 0;
    ref_offset = result->ref_begin1;
    read_offset = result->read_begin1;

    for (c = 0; c < result->cigarLen; ++c) {
        cigar_int = result->cigar[c];
        cigar_len = (uint32_t)(cigar_int >> BAM_CIGAR_SHIFT);
        cigar_char = (cigar_int & 0xfU) > 8 ? 'M': MAPSTR[cigar_int & 0xfU];

        hit->length += cigar_len;
        
        if(cigar_char == 'I') {
            hit->gapopen++;
            read_offset += cigar_len;
        }
        else if(cigar_char == 'D') {
            hit->gapopen++;
            ref_offset += cigar_len;
        }
        else {
            for(i = 0; i < cigar_len; ++i) {
                if(ref[ref_offset] == num[read_offset]) {
                    identical++;
                }
                else {
                    hit->mismatch++;
                }
                ref_offset++;
                read_offset++;
            }
        }
    }

    hit->identity = identical / (float)hit->length;

    hit->rawscore = result->score1;
    hit->bitscore = (opt->lambda * hit->rawscore - log(opt->k)); // / M_LN2;
    hit->evalue = corrected_q * corrected_db * exp(-hit->bitscore);
    hit->evalue = hit->evalue > 1.0e-180 ? hit->evalue : 0.0;
    hit->bitscore /= M_LN2; // nats --> bits
    
    //if(correction > hit->length) {
    //    hit->evalue = INFINITY;
    //}

    hit->is_aligned = 1;
    hit->effective_length = seq_len(s1) - correction;

    free(num);
    free(ref);
    align_destroy(result);
    init_destroy(profile);
}


#define TESTING_TESTINGXXX
#ifdef TESTING_TESTING
#include "fasta.h"

int main(int argc, char** argv) {
    seq_t *s1, *s2;
    options_t par;
    hit_t hit;
    s_profile* prof;
    int count = 0;

    int s1_start, s2_start, length;

    fasta_t* f;

    f = fasta_alloc("bad_human.fasta");
    s1 = fasta_next(f, NULL);
    fasta_free(f);

    seq_2internal(s1);

/*
    f = fasta_alloc("Q8WZ42.fasta");
    s2 = fasta_next(f, NULL);
    fasta_free(f);
*/
/*
    f = fasta_alloc("A0A0C2DC85.fasta");
    s2 = fasta_next(f, NULL);
    fasta_free(f);
*/

    par.gap_open = DEFAULT_GAPOPEN;
    par.gap_extend = DEFAULT_GAPEXTEND;
    par.substitution = &BLOSUM62;
    par.lambda = DEFAULT_LAMBDA;
    par.k = DEFAULT_K;
    par.H = DEFAULT_H;
    par.superfast = 1;
    par.db_size = 26172347515;

    prof = get_ssw_profile(s1, &par);

    f = fasta_alloc("bad.fa");
    while((s2 = fasta_next(f, NULL)) != NULL) {

        seq_2internal(s2);

        //lcp(s1, s2, &s1_start, &s2_start, &length);

        //printf("%.*s %d-%d %d-%d %d\n", seq_idlen(s2), seq_idptr(s2), s1_start, s1_start + length, s2_start, s2_start + length, length);
        //printf("%d\n", length);

        //align_using_ssw(s1, s2, &hit, &par);
        align_using_ssw_profile(prof, s2, &hit, &par);
        //printf("%.*s %f %f\n", seq_idlen(s2), seq_idptr(s2), hit.evalue, hit.bitscore);

        printf("%.*s %d %f %f %d-%d %d-%d\n", seq_idlen(s2), seq_idptr(s2), hit.rawscore, hit.bitscore, hit.evalue, hit.qstart, hit.qend, hit.sstart, hit.send);
        
        seq_free(s2);

        //if(++count == 200)
        //    break;
    }

    seq_free(s1);
    //seq_free(s2);

    return EXIT_SUCCESS;
}
#endif


/*
int main(int argc, char** argv) {
    seq_t *s1, *s2;
    options_t par;
    hit_t hit;

    s1 = seq_alloc();
    seq_id(s1, "protein1", 8);
    seq_seq(s1, "HEAGAWGHEE", 10);

    s2 = seq_alloc();
    seq_id(s2, "protein2", 8);
    seq_seq(s2, "PAWHEAE", 7);

    seq_write(s1, SEQ_FASTA, stdout);
    seq_write(s2, SEQ_FASTA, stdout);

    par.gap_open = 10;
    par.gap_extend = 1;
    par.substitution = &BLOSUM50;
    par.lambda = DEFAULT_LAMBDA;
    par.k = DEFAULT_K;
    par.db_size = 1e3;

    align_using_ssw(s1, s2, &hit, &par);

    printf("raw score: %d\n", hit.rawscore);
    printf("bit score: %f\n", hit.bitscore);
    printf("E value: %f\n", hit.evalue);
    printf("qto: %d, qfrom: %d\n", hit.qend, hit.qstart);
    printf("sto: %d, sfrom: %d\n", hit.send, hit.sstart);
    printf("identity: %f\n", hit.identity);
    printf("mismatches: %d\n", hit.mismatch);

    seq_free(s1);
    seq_free(s2);

    return EXIT_SUCCESS;
}
*/

