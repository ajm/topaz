#ifndef H_ALIGNMENT
#define H_ALIGNMENT

#include "seq.h"
#include "hit.h"
#include "options.h"

// TODO IMPLEMENT DEFAULT OPTIONS FROM THESE WEBPAGES
// blast options + defaults https://www.arabidopsis.org/Blast/BLASToptions.jsp
// blast sub matrices for length http://www.ncbi.nlm.nih.gov/blast/html/sub_matrix.html

// pointer matrix 
#define SW_NULL             0
#define SW_H_H_DIAG_MATCH   1
#define SW_H_H_DIAG_MISS    8
#define SW_H_E_JUMP         2
#define SW_H_F_JUMP         3
#define SW_E_E_LEFT         4
#define SW_E_H_LEFT         5
#define SW_F_F_UP           6
#define SW_F_H_UP           7


int _forward_pass(  int8_t* s1, int8_t* s2,
                    size_t s1_len, size_t s2_len,
                    int* E, uint16_t* trE,
                    int* F, uint16_t* trF,
                    int* H, uint16_t* trH,
                    hit_t* hit, options_t* opt);
void _backtrack(uint16_t* trE, uint16_t* trF, uint16_t* trH, int pos, int max_x, hit_t* hit);
int gotoh(seq_t* s1, seq_t* s2, hit_t* hit, options_t* p);
int ungapped_rawscore(char* s1, char* s2, uint32_t length, options_t* opt);

#endif

