#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "seq.h"
#include "lexicographical.h"


void buffer_alloc(buffer_t* buf, size_t len) {
    buf->buf = (char*) malloc(len);
    buf->len = len;
}

void buffer_resize(buffer_t* buf, size_t len) {
    if(buf->len < len) {
        buf->buf = (char*) realloc(buf->buf, len);
        buf->len = len;
    }
}

void buffer_copy(buffer_t* buf, char* cbuf, size_t len) {
    buffer_resize(buf, len);
    memcpy(buf->buf, cbuf, len);
}

char* buffer_ptr(buffer_t* buf) {
    return buf->buf;
}

void buffer_free(buffer_t* buf) {
    if(buf) {
        if(buf->buf) {
            free(buf->buf);
        }
        //free(buf);
    }
}

seq_t* seq_alloc() {
    seq_t *s = (seq_t*) malloc(sizeof(seq_t));
    buffer_alloc(&s->id_buf, BUFSIZ);
    buffer_alloc(&s->seq_buf, BUFSIZ);
    s->id_len = s->seq_len = 0;
    return s;
}

void seq_free(seq_t* s) {
    if(s) {
        buffer_free(&s->id_buf);
        buffer_free(&s->seq_buf);
        free(s);
    }
}

int seq_id(seq_t* s, char* buf, size_t len) {
    int i, start;
    char* sbuf;

    if(len == 0) {
        s->id_len = 0;
        return 0;
    }

    buffer_resize(&s->id_buf, len);
    sbuf = buffer_ptr(&s->id_buf);

    start = buf[0] == '>' ? 1 : 0;

    for(i = start; i < len; ++i) {
        if((buf[i] == ' ') || (buf[i] == '\t'))
            break;

        sbuf[i-start] = buf[i];
    }

    s->id_len = i - start;
    return s->id_len;
}

int seq_seq(seq_t* s, char* buf, size_t len) {
    int i, j;
    char* sbuf;

    buffer_resize(&s->seq_buf, len);
    sbuf = buffer_ptr(&s->seq_buf);

    for(i = 0, j = 0; i < len; ++i) {
        if(buf[i] != '\n') {
            sbuf[j++] = buf[i];
        }
    }

    s->seq_len = j;
    return j;
}

int seq_2internal(seq_t* s) {
    int i;
    char *buf;

    buf = buffer_ptr(&s->seq_buf);

    for(i = 0; i < s->seq_len; ++i)
        buf[i] = aa2internal(buf[i]);

    return 0;
}

int seq_len(seq_t* s) {
    return s->seq_len;
}

char* seq_ptr(seq_t* s) {
    return buffer_ptr(&s->seq_buf);
}

char seq_char(seq_t* s, size_t pos) {
    return buffer_ptr(&s->seq_buf)[pos];
}

void seq_write(seq_t* s, int bitmask, FILE* stream) {

    if(bitmask & SEQ_FASTA) {
        fprintf(stream, ">%.*s\n%.*s\n", 
            (int) s->id_len, buffer_ptr(&s->id_buf), 
            (int) s->seq_len, buffer_ptr(&s->seq_buf));
    }
    
    if(bitmask & SEQ_ID)
        fprintf(stream, "%.*s\n", (int) s->id_len, buffer_ptr(&s->id_buf));
    
    if(bitmask & SEQ_SEQ)
        fprintf(stream, "%.*s\n", (int) s->seq_len, buffer_ptr(&s->seq_buf));
}

char* seq_idptr(seq_t *s) {
    return buffer_ptr(&s->id_buf);
}

int seq_idlen(seq_t* s) {
    return s->id_len;
}

