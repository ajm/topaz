#define _GNU_SOURCE
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <string.h>
#include <errno.h>
#include <stdint.h>

#include "hit.h"


void *fmmap_ro(char *fname, uint64_t* size) {
    int fd;
    struct stat sbuf;
    void *ptr;

    if((fd = open(fname, O_RDONLY)) < 0) {
        return NULL;
    }

    if(fstat(fd, &sbuf) != 0) {
        return NULL;
    }

    *size = sbuf.st_size;

    if((ptr = mmap((caddr_t) 0, *size, PROT_READ, MAP_FILE | MAP_PRIVATE, fd, 0)) == MAP_FAILED) {
        return NULL;
    }
/*
    if(madvise(ptr, *size, MADV_SEQUENTIAL) != 0) {
        return NULL;
    }
*/
    close(fd);

    return ptr;
}

void *fmmap_rw(char *fname, uint64_t size) {
    int fd;
    void *ptr;

    if((fd = open(fname, O_RDWR | O_CREAT, 0666)) < 0) {
        return NULL;
    }

    if(ftruncate(fd, size) != 0) {
        return NULL;
    }

    if((ptr = mmap((caddr_t) 0, size, PROT_READ | PROT_WRITE, MAP_FILE | MAP_SHARED, fd, 0)) == MAP_FAILED) {
        return NULL;
    }

    close(fd);

    return ptr;
}

void funmap(void *addr, size_t len, char *fname) {

    fprintf(stderr, "Syncing %s to disk...\n", fname);

    if(msync(addr, len, MS_SYNC) != 0) {
        fprintf(stderr, "ERROR: msync %s failed, %s\n", fname, strerror(errno));
        exit(EXIT_FAILURE);
    }

    if(munmap(addr, len) != 0) {
        fprintf(stderr, "ERROR: %s munmap failed, %s\n", fname, strerror(errno));
        exit(EXIT_FAILURE);
    }
}

char* mk_fname(char *prefix, char *ext) {
    char* s;
    ssize_t plen, elen;

    plen = strlen(prefix);
    elen = strlen(ext);

    s = (char*) malloc(plen + elen + 1);
    
    memcpy(s, prefix, plen);
    memcpy(s + plen, ext, elen);

    s[plen + elen] = '\0';

    return s;
}

FILE* fopen_pf(char* prefix, char* ext) {
    FILE *tmp;
    char *buf;

    buf = mk_fname(prefix, ext);
    tmp = fopen(buf, "w");

    free(buf);
    return tmp;
}

/*
void rtrim(char* s, ssize_t* len) {
    int i, end;

    for(i = *len - 1, end = 0; i > -1; --i) {
        if(! isspace(s[i]))
            break;
        
        end++;
    }

    s[*len-end] = '\0';
    *len -= end;
}
*/

uint64_t strncount(char *s, char c, size_t n) {
    uint64_t i;
    uint64_t counter;

    for(i = 0, counter = 0; i < n; ++i) {
        if(s[i] == '\0')
            break;

        if(s[i] == c)
            counter++;
    }

    return counter;
}

int compare_bitscore(const void *l, const void *r) {
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

