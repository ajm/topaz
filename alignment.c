#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>

#include "seq.h"
#include "alignment.h"
#include "options.h"
#include "hit.h"
#include "substitution.h"
#include "lexicographical.h"

// 2D -> 1D conversion
// ------------------> X = max_x
// |
// |    (x,y)
// |
// V
// Y
#define POS(x, y, max_x) (((y) * (max_x)) + (x))
#define MAX2(x, y) (((x) > (y)) ? (x) : (y))
#define MAX3(x, y, z) MAX2(MAX2(x, y), z)
#define MAX4(x, y, z, w) MAX2(MAX2(x, y), MAX2(z, w))


void printdp32(int* m, size_t dim1, size_t dim2) {
    int i, j, pos = 0;

    for(i = 0; i < dim2; ++i) {         // rows
        for(j = 0; j < dim1; ++j) {     // columns
            printf("%3d ", m[pos++]);
        }
        printf("\n");
    }
}

void printdp8(uint8_t* m, size_t dim1, size_t dim2) {
    int i, j, pos = 0;

    for(i = 0; i < dim2; ++i) {         // rows
        for(j = 0; j < dim1; ++j) {     // columns
            printf("%3d ", m[pos++]);
        }
        printf("\n");
    }
}

int _forward_pass(  int8_t* seq1, int8_t* seq2,
                    size_t seq1_len, size_t seq2_len,
                    int* E, uint16_t* trE,
                    int* F, uint16_t* trF,
                    int* H, uint16_t* trH,
                    hit_t* hit, options_t* par) {
    
    int i, j;
    int8_t c1, c2;
    int tmp1, tmp2;

    size_t max_x = seq1_len + 1;
    size_t max_y = seq2_len + 1;
    int pos = seq1_len + 2;
    int max_score = 0;
    int max_pos = 0;

    for(i = 1; i < max_y; ++i) {      // rows
        //c1 = seq_char(s2, i-1);
        c1 = seq2[i-1];

        for(j = 1; j < max_x; ++j) {  // columns
            //c2 = seq_char(s1, j-1);
            c2 = seq1[j-1];

            // E matrix
            tmp1 = E[pos-1] - par->gap_extend;
            tmp2 = H[pos-1] - par->gap_open; // - par->gap_extend;
            if(tmp1 > tmp2) {
                E[pos] = tmp1;
                trE[pos] = SW_E_E_LEFT;
            }
            else {
                E[pos] = tmp2;
                trE[pos] = SW_E_H_LEFT;
            }

            // F matrix
            tmp1 = F[pos-max_x] - par->gap_extend;
            tmp2 = H[pos-max_x] - par->gap_open; // - par->gap_extend;
            if(tmp1 > tmp2) {
                F[pos] = tmp1;
                trF[pos] = SW_F_F_UP;
            }
            else {
                F[pos] = tmp2;
                trF[pos] = SW_F_H_UP;
            }

            // H matrix
            tmp1 = H[pos-max_x-1] + (*par->substitution)[(c1 * 25) + c2];

            H[pos] = trH[pos] = 0;

            if(tmp1 > H[pos]) {
                H[pos] = tmp1;
                trH[pos] = (c1 == c2) ? SW_H_H_DIAG_MATCH : SW_H_H_DIAG_MISS;
            }

            if(E[pos] > H[pos]) {
                H[pos] = E[pos];
                trH[pos] = SW_H_E_JUMP;
            }
            
            if(F[pos] > H[pos]) {
                H[pos] = F[pos];
                trH[pos] = SW_H_F_JUMP;
            }


            if(H[pos] > max_score) {
                max_score = H[pos];
                max_pos = pos;
            }

            ++pos;
        }

        ++pos;
    }

    hit->rawscore = max_score;
    hit->send = max_pos / max_x;
    hit->qend = max_pos % max_x;

/*
    printf("========= E ========\n");
    printdp32(E, max_x, max_y);
    printf("========= F ========\n");
    printdp32(F, max_x, max_y);
    printf("========= H ========\n");
    printdp32(H, max_x, max_y);
    printf("=========trE========\n");
    printdp8(trE, max_x, max_y); 
    printf("=========trF========\n");
    printdp8(trF, max_x, max_y); 
    printf("=========trH========\n");
    printdp8(trH, max_x, max_y);
*/

    return max_pos;
}

void _backtrack(uint16_t* trE, uint16_t* trF, uint16_t* trH, int pos, int max_x, hit_t* hit) {
    int tmp = trH[pos];
    char matrix = 'h';
    int total = 0;
    int identical = 0;
    int mismatch = 0;
    int gapopen = 0;
    
    do {
        switch(matrix) {
            case 'h' :
                switch(tmp) {
                    case SW_H_H_DIAG_MATCH :
                        pos -= (max_x + 1);
                        identical++;
                        total++;
                        break;
                    case SW_H_H_DIAG_MISS :
                        pos -= (max_x + 1);
                        mismatch++;
                        total++;
                        break;
                    case SW_H_E_JUMP :
                        gapopen++;
                        matrix = 'e';
                        break;
                    case SW_H_F_JUMP :
                        gapopen++;
                        matrix = 'f';
                        break;
                    default :
                        goto backtrack_error;
                }
                break;

            case 'e' :
                pos--;
                total++;

                if(tmp == SW_E_H_LEFT)
                    matrix = 'h';

                break;

            case 'f' :
                pos -= max_x;
                total++;

                if(tmp == SW_F_H_UP)
                    matrix = 'h';
                
                break;

            default :
                goto backtrack_error;
        }


        switch(matrix) {
            case 'h' :
                tmp = trH[pos];
                break;
            case 'e' :
                tmp = trE[pos];
                break;
            case 'f' :
                tmp = trF[pos];
                break;
            default :
                goto backtrack_error;
        }

    } while(tmp != SW_NULL);

    // pos at start of function tells us sto and qto
    // pos at end of function tells us sfrom and qfrom

    hit->sstart = (pos / max_x) + 1;
    hit->qstart = (pos % max_x) + 1;
    hit->mismatch = mismatch;
    hit->gapopen = gapopen;
    hit->length = total;
    hit->identity = identical / (float) total;

    return;

backtrack_error:
    fprintf(stderr, "Error: backtrack failed, illegal condition\n");
    abort();
}

void _initialise(seq_t* s1, seq_t* s2, int* E, int* F, int* H, uint16_t* trE, uint16_t* trF, uint16_t* trH, options_t* par) {
    int i, tmp;
    size_t max_x = seq_len(s1) + 1;
    size_t max_y = seq_len(s2) + 1;

    E[0] = F[0] = H[0] = 0;
    trE[0] = trF[0] = trH[0] = 0;

    for(i = 1; i < max_x; ++i) {
        E[i] = F[i] = H[i] = 0;
        trE[i] = trF[i] = trH[i] = 0;
    }

    for(i = 1; i < max_y; ++i) {
        tmp = i * max_x;
        E[tmp] = F[tmp] = H[tmp] = 0;
        trE[tmp] = trF[tmp] = trH[tmp] = 0;
    }
}

// TODO make more uniform interface, for example, gapped_rawscore()
//      and upgapped_rawscore(), with wrappers for different types + adding
//      results to hit_t structs
int ungapped_rawscore(char* s1, char* s2, uint32_t length, options_t* opt) {
    int i, rawscore = 0;
    int c1, c2;
    int ind;

    //fprintf(stderr, "[1] %.*s\n[2] %.*s\n\n", length, s1, length, s2);

    for(i = 0; i < length; ++i) {
        c1 = (int) s1[i];
        c2 = c1 ; //(int) s2[i];

#ifdef USELEXICOGRAPHICAL
        ind = (aa_table[c1] * 25) + aa_table[c2];
#else
        ind = (aa_table2[c1] * 25) + aa_table2[c2];
#endif

        rawscore += (*opt->substitution)[ind];
    }

    return rawscore;
}

int ungapped_length(char* s1, char* s2, uint32_t max_length, options_t* opt) {
    int i, rawscore = 0;
    int c1, c2, ind;

    for(i = 0; i > max_length; ++i) {
        c1 = (int) s1[i];
        c2 = (int) s2[i];

#ifdef USELEXICOGRAPHICAL
        ind = (aa_table[c1] * 25) + aa_table[c2];
#else
        ind = (aa_table2[c1] * 25) + aa_table2[c2];
#endif

        rawscore += (*opt->substitution)[ind];

        if(rawscore < 0) {
            return i-1;
        }

        if(rawscore == 0) {
            return i;
        }
    }

    return max_length;
}

int gotoh(seq_t* s1, seq_t* s2, hit_t* hit, options_t* opt) {
    int *E, *F, *H;
    uint16_t *trE, *trF, *trH;
    int msize = (seq_len(s1) + 1) * (seq_len(s2) + 1);
    int max_pos, i;
    double correction, corrected_q, corrected_db;
    int8_t *seq1, *seq2;

    // allocate
    seq1 = (int8_t*) malloc(seq_len(s1));
    seq2 = (int8_t*) malloc(seq_len(s2));

    E = (int*) malloc(msize * sizeof(*E));
    F = (int*) malloc(msize * sizeof(*F));
    H = (int*) malloc(msize * sizeof(*H));

    trE = (uint16_t*) malloc(msize * sizeof(*trE));
    trF = (uint16_t*) malloc(msize * sizeof(*trF));
    trH = (uint16_t*) malloc(msize * sizeof(*trH));

    // initialise
    _initialise(s1, s2, E, F, H, trE, trF, trH, opt);
    
    for(i = 0; i < seq_len(s1); ++i) {
#ifdef USELEXICOGRAPHICAL
        seq1[i] = aa_table[(int) seq_char(s1, i)];
#else
        seq1[i] = aa_table2[(int) seq_char(s1, i)];
#endif
    }

    for(i = 0; i < seq_len(s2); ++i) {
#ifdef USELEXICOGRAPHICAL
        seq2[i] = aa_table[(int) seq_char(s2, i)];
#else
        seq2[i] = aa_table2[(int) seq_char(s2, i)];
#endif
    }

    // alignment
    max_pos = _forward_pass(seq1, seq2, 
                            seq_len(s1), seq_len(s2),
                            E, trE, 
                            F, trF, 
                            H, trH, 
                            hit, opt);

    _backtrack(trE, trF, trH, max_pos, seq_len(s1)+1, hit);
    
    // statistics
    correction = log(opt->k * seq_len(s1) * opt->db_size) / opt->H;
    corrected_q = fmax(seq_len(s1) - correction, 1.0); // / opt->k);
    corrected_db = fmax(opt->db_size - (correction * opt->db_proteins), 1.0);

    hit->bitscore = (opt->lambda * hit->rawscore - log(opt->k)); // / M_LN2;
    hit->evalue = corrected_q * corrected_db * exp(-hit->bitscore);
    hit->evalue = hit->evalue > 1.0e-180 ? hit->evalue : 0.0;
    hit->bitscore /= M_LN2; // nats --> bits

    if(correction > hit->length) {
        hit->evalue = INFINITY;
    }

    hit->is_aligned = 1;
    hit->effective_length = seq_len(s1) - correction;

    // clean up
    free(trE);
    free(trF);
    free(trH);

    free(E);
    free(F);
    free(H);

    free(seq1);
    free(seq2);

    return 0;
}

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

    gotoh(s1, s2, &hit, &par);

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

