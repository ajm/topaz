#ifndef H_VECTOR
#define H_VECTOR

#include <stdint.h>

#define VECTOR_SCALE_FACTOR 2

typedef struct {
    void** contents;
    size_t size;
    size_t length;

} vector_t;

vector_t* vector_alloc();
void vector_free(vector_t* v);
void vector_push(vector_t* v, void* item);
void vector_insert(vector_t* v, void* item, uint32_t pos);
uint32_t vector_index(vector_t* v, void* item, int (*compare)(void* item1, void* item2));
void* vector_get(vector_t* v, uint32_t pos);
void* vector_pop(vector_t* v);
size_t vector_length(vector_t* v);
void vector_reset(vector_t* v);
void vector_delete(vector_t* v, uint32_t index);
int vector_empty(vector_t* v);

#endif

