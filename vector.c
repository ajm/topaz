#include <stdio.h>
#include <stdlib.h>
#include "vector.h"

vector_t* vector_alloc() {
    vector_t* tmp;

    tmp = (vector_t*) malloc(sizeof(vector_t));

    tmp->size = BUFSIZ / sizeof(void**);
    tmp->contents = malloc(sizeof(void**) * tmp->size);
    tmp->length = 0;

    return tmp;
}

// current vector implementation can only hold pointers to
// structs without pointers that need to be freed (i.e. I just
// make a call to free())
void vector_free(vector_t* v) {
    int i;

    if(v) {
        if(v->contents) {
            for(i = 0; i < v->length; ++i) {
                free(v->contents[i]);
            }

            free(v->contents);
        }

        free(v);
    }
}

void _increase_allocation(vector_t *v) {
    void **new_contents;
    new_contents = realloc(v->contents, sizeof(void**) * v->size * VECTOR_SCALE_FACTOR);
    v->contents = new_contents;
    v->size *= VECTOR_SCALE_FACTOR;
}

void vector_push(vector_t* v, void* item) {
    if(v->length == v->size) 
        _increase_allocation(v);

    v->contents[v->length++] = item;
}

void* vector_get(vector_t *v, uint32_t pos) {
    return pos < v->length ? v->contents[pos] : NULL;
}

// return index that item would be at if we assume that
// it would be inserted in sorted order
uint32_t vector_index(vector_t *v, void *item, int (*compare)(void* item1, void* item2)) {
    int i, size, half, comp;
    void *tmp;

    for(i = 0, size = v->length, half = size >> 1;
        0 < size;
        size = half, half >>= 1) {

        tmp = v->contents[i + half];
        comp = compare(item, tmp);

        if(comp == -1) {
            continue;
        }

        if(comp == 1) {
            i += half + 1;
            half -= (size & 1) ^ 1;
            continue;
        }

        break;
    }

    return i + half;
}

void vector_insert(vector_t* v, void* item, uint32_t pos) {
    int i;

    if(v->length == v->size) {
        _increase_allocation(v);
    }

    // shuffle everything up one
    for(i = v->length; i > pos; --i) {
        v->contents[i] = v->contents[i-1];
    }

    // add item at position
    v->contents[pos] = item;
    v->length++;
}

void* vector_pop(vector_t* v) {
    return v->length ? v->contents[--(v->length)] : NULL;
}

size_t vector_length(vector_t* v) {
    return v->length;
}

int vector_empty(vector_t* v) {
    return vector_length(v) == 0;
}

void vector_reset(vector_t* v) {
    v->length = 0;
}

void vector_delete(vector_t* v, uint32_t index) {
    int i;

    free(v->contents[index]);

    for(i = index; i < (v->length - 1); ++i) {
        v->contents[i] = v->contents[i+1];
    }

    v->length--;
}

/*
int main(int argc, char** argv) {
    vector_t* v;
    int a = 1, b = 2, c = 3;
    int* tmp;

    v = vector_alloc();

    vector_push(v, &a);
    vector_push(v, &b);
    vector_push(v, &c);

    while((tmp = (int*) vector_pop(v)) != NULL) {
        printf("pop %d\n", *tmp);
    }

    vector_free(v);

    return EXIT_SUCCESS;
}
*/

