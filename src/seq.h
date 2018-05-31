#ifndef H_SEQ
#define H_SEQ

#include <stdio.h>
#include <stddef.h>


#define SEQ_ID      1
#define SEQ_SEQ     2
#define SEQ_FASTA   4

typedef struct {
    char *buf;
    int len;

} buffer_t;

void buffer_alloc(buffer_t* buf, size_t len);
void buffer_resize(buffer_t* buf, size_t len);
void buffer_copy(buffer_t* buf, char* cbuf, size_t len);
char* buffer_ptr(buffer_t* buf);
void buffer_free(buffer_t* buf);

typedef struct {
    buffer_t id_buf;
    buffer_t seq_buf;

    int id_len;
    int seq_len;

} seq_t;

seq_t* seq_alloc();
void seq_free(seq_t* s);
int seq_id(seq_t* s, char* buf, size_t len);
int seq_seq(seq_t* s, char* buf, size_t len, int check);
int seq_2internal(seq_t* s);
int seq_len(seq_t* s);
char* seq_ptr(seq_t* s);
char seq_char(seq_t* s, size_t pos);
void seq_write(seq_t* s, int bitmask, FILE* stream);

char* seq_idptr(seq_t* s);
int seq_idlen(seq_t* s);

#endif

