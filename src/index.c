#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <errno.h>
#include <sys/mman.h>

#include <divsufsort64.h>

#include "generic.h"
#include "utils.h"
#include "demux.h"
#include "db.h"


#ifdef USE_LCP
static int create_lcp(char *prefix, char *T, saidx64_t Tsize, saidx64_t *SA, saidx64_t *ISA, saidx64_t SAsize);
static void build_lcp(char *T, saidx64_t *SA, saidx64_t *ISA, saidx64_t *LCP, saidx64_t n);
static saidx64_t build_lcplr(saidx64_t *LCP, saidx64_t *LCPLC, saidx64_t *LCPCR, saidx64_t start, saidx64_t size);
#endif
static int create_index(char* prefix, uint64_t Tsize);
//static void print_sa(char *T, saidx64_t *SA, saidx64_t *LCP, saidx64_t *LCPLC, saidx64_t *LCPCR, saidx64_t n);

void serif_index(char *fname, char *prefix) {
    char *T;
    char *T_fname, 
         *SA_fname, 
         *ISA_fname;
    saidx64_t *SA, *ISA;
    uint64_t Tsize, SAsize;


    demultiplex_fasta(fname, prefix);

    T_fname = mk_fname(prefix, ".psq");
    SA_fname = mk_fname(prefix, ".SA");
    ISA_fname = mk_fname(prefix, ".ISA");

    // mmap text
    if((T = (char*) fmmap_ro(T_fname, &Tsize)) == NULL) {
        fprintf(stderr, "Error: problem mmapping %s, %s\n", T_fname, strerror(errno));
        exit(EXIT_FAILURE);
    }

    fprintf(stderr, "mmapped T (%s), file size is %lld bytes\n", T_fname, Tsize);

    SAsize = Tsize * sizeof(saidx64_t);

    // mmap suffix array
    if((SA = (saidx64_t*) fmmap_rw(SA_fname, SAsize)) == NULL) {
        fprintf(stderr, "Error: problem mmapping %s, %s\n", SA_fname, strerror(errno));
        exit(EXIT_FAILURE);
    }

    fprintf(stderr, "mmapped SA (%s), file size is %lld bytes\n", SA_fname, SAsize);

    // mmap inverse suffix array
    if((ISA = (saidx64_t*) fmmap_rw(ISA_fname, SAsize)) == NULL) {
        fprintf(stderr, "Error: problem mmapping %s, %s\n", ISA_fname, strerror(errno));
        exit(EXIT_FAILURE);
    }

    fprintf(stderr, "mmapped ISA (%s), file size is %lld bytes\n", ISA_fname, SAsize);

    /* build suffix array */
    if(divsufsort64((unsigned char*) T, SA, Tsize) < 0) {
        fprintf(stderr, "divsufsort64 failed\n");
        exit(EXIT_FAILURE);
    }

#ifdef USE_LCP
    if(create_lcp(prefix, T, Tsize, SA, ISA, SAsize) != 0) {
        fprintf(stderr, "Error: LCP arrays could not be built, %s\n", strerror(errno));
        exit(EXIT_FAILURE);
    }
#endif

    fprintf(stderr, "creating protein offset index...\n");
    if(create_index(prefix, Tsize) != 0) {
        fprintf(stderr, "Error: could not create protein index\n");
        exit(EXIT_FAILURE);
    }

    fprintf(stderr, "finishing up...\n");

    funmap(SA, SAsize, SA_fname);
    funmap(ISA, SAsize, ISA_fname);

    free(SA_fname);
    free(ISA_fname);

    fprintf(stderr, "Done!\n");
}

uint64_t _num_proteins(char* prefix) {
    char *phr_fname, *headers;
    uint64_t Hsize, Psize;

    phr_fname = mk_fname(prefix, ".phr");

    if((headers = (char*) fmmap_ro(phr_fname, &Hsize)) == NULL) {
        fprintf(stderr, "Error: could not open %s, %s\n", phr_fname, strerror(errno));
        exit(EXIT_FAILURE);
    }   

    Psize = strncount(headers, '\n', Hsize);

    funmap(headers, Hsize, phr_fname);
    free(phr_fname);

    return Psize;
}

static
int create_index(char* prefix, uint64_t Tsize) {
    size_t linecap;
    ssize_t linelen;
    uint64_t pstart, hstart, linenum = 0, length;
    int plen, hlen, nmatch;
    char *pin_fname, *ind_fname, *line = NULL;
    FILE *f;
    position_t* index;
    uint64_t Psize, Isize;
    

    Psize = _num_proteins(prefix);
    Isize = Psize * sizeof(position_t);


    pin_fname = mk_fname(prefix, ".pin");
    ind_fname = mk_fname(prefix, ".ind");

    if((f = fopen(pin_fname, "r")) == NULL) {
        return -1;
    }

    index = (position_t*) fmmap_rw(ind_fname, Isize);

    while((linelen = getline(&line, &linecap, f)) > 0) {
        if((nmatch = sscanf(line, "%lld %d %lld %d\n", &pstart, &plen, &hstart, &hlen)) != 4) {
            fprintf(stderr, "Error: line %lld of %s contained %d elements (expected 4)\n", linenum + 1, pin_fname, nmatch);
            exit(EXIT_FAILURE);
        }

        index[linenum].header_start = hstart;
        index[linenum].protein_start = pstart;
        index[linenum].header_len = hlen;
        index[linenum].protein_len = plen;

        linenum++;
    }

    if(linenum != Psize) {
        fprintf(stderr, "Error: expected %lld index entries, but read %lld from %s\n", Psize, linenum, pin_fname);
        exit(EXIT_FAILURE);
    }

    length = index[Psize-1].protein_start + index[Psize-1].protein_len;

    if(Tsize != (length + 1)) {
        fprintf(stderr, "Error: index length is %lld, but T length is %lld ...\n", length, Tsize);
        exit(EXIT_FAILURE);
    }

    fclose(f);

    funmap(index, Isize, ind_fname);
    free(pin_fname);
    free(ind_fname);
    free(line);
    return 0;
}

#ifdef USE_LCP
static
int create_lcp(char *prefix, char *T, saidx64_t Tsize, saidx64_t *SA, saidx64_t *ISA, saidx64_t SAsize) {
    char *LCP_fname, *LCPLC_fname, *LCPCR_fname;
    saidx64_t *LCP, *LCPLC, *LCPCR;

    LCP_fname = mk_fname(prefix, ".LCP");
    LCPLC_fname = mk_fname(prefix, ".LCPLC");
    LCPCR_fname = mk_fname(prefix, ".LCPCR");

    if((LCP = (saidx64_t*) fmmap_rw(LCP_fname, SAsize)) == NULL) {
        return -1;
    }
   
    fprintf(stderr, "mmapped LCP (%s), file size is %lld bytes\n", LCP_fname, SAsize);

    if((LCPLC = (saidx64_t*) fmmap_rw(LCPLC_fname, SAsize)) == NULL) {
        return -1;
    }

    fprintf(stderr, "mmapped LCP-LC (%s), file size is %lld bytes\n", LCPLC_fname, SAsize);

    if((LCPCR = (saidx64_t*) fmmap_rw(LCPCR_fname, SAsize)) == NULL) {
        return -1;
    }

    fprintf(stderr, "mmapped LCP-CR (%s), file size is %lld bytes\n", LCPCR_fname, SAsize);

    
    build_lcp(T, SA, ISA, LCP, Tsize);
    build_lcplr(LCP, LCPLC, LCPCR, 0, Tsize);

    //print_sa(T, SA, LCP, LCPLC, LCPCR, Tsize);
    
    funmap(LCP, SAsize, LCP_fname);
    funmap(LCPLC, SAsize, LCPLC_fname);
    funmap(LCPCR, SAsize, LCPCR_fname);
   
    free(LCP_fname);
    free(LCPLC_fname);
    free(LCPCR_fname);

    return 0;
}

static 
void build_lcp(char *T, saidx64_t *SA, saidx64_t *ISA, saidx64_t *LCP, saidx64_t n) {
    saidx64_t i, k, h, maxh = 0;

    for(i = 0; i < n; ++i) {
        ISA[SA[i]] = i;
    }

    for(i = 0, h = 0; i < n; ++i) {
        if(ISA[i] > 0) {
            k = SA[ISA[i] - 1];

            while(T[i + h] == T[k + h]) {
                h++;
            }

            LCP[ISA[i]] = h;

            if(h > maxh) {
                maxh = h;
            }

            if(h > 0) {
                h--;
            }
        }
    }

    if(maxh > LCP_MAX) {
        fprintf(stderr, "Critical Error: values in the LCP array cannot exceed %d (max. LCP = %lld), please _recompile_ with greater precision\n", LCP_MAX, maxh);
        exit(EXIT_FAILURE);
    }
}

static 
saidx64_t build_lcplr(saidx64_t *LCP, saidx64_t *LCPLC, saidx64_t *LCPCR, saidx64_t start, saidx64_t size) {
    saidx64_t half;

    if(size == 0)
        return LCP[start];

    half = size >> 1;

    LCPLC[start + half] = build_lcplr(LCP, LCPLC, LCPCR, start, half);
    LCPCR[start + half] = build_lcplr(LCP, LCPLC, LCPCR, start + half + 1, half - ((size & 1) ^ 1));

    return MIN(LCPLC[start + half], LCPCR[start + half]);
}
#endif

/*
static 
void print_sa(char *T, saidx64_t *SA, saidx64_t *LCP, saidx64_t *LCPLC, saidx64_t *LCPCR, saidx64_t n) {
    saidx64_t i;

    printf("  i LCP LCP-LC LCP-CR SA T\n");

    for(i = 0; i < n; ++i) {
        printf("%8lld%8lld%8lld%8lld%8lld  %s", i, LCP[i], LCPLC[i], LCPCR[i], SA[i], T + SA[i]);
    }
}
*/

