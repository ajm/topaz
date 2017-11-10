#ifndef H_UTILS
#define H_UTILS

#include <stdint.h>

void *fmmap_ro(char *fname, uint64_t* size);
void *fmmap_rw(char *fname, uint64_t size);
void funmap(void *addr, size_t len, char *fname);
char* mk_fname(char *prefix, char *ext);
FILE* fopen_pf(char* prefix, char* ext);
//void rtrim(char* s, ssize_t* len);
uint64_t strncount(char *s, char c, size_t n);

#endif

