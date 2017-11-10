#ifndef H_WORK_QUEUE
#define H_WORK_QUEUE

#include <pthread.h>
#include "rbtree.h"
#include "query.h"


typedef struct witem {
    void* item;
    struct witem* next;
    struct witem* prev;

} witem_t;

typedef struct {
    witem_t* head;
    witem_t* tail;
    int length;
    int max_length;
    pthread_mutex_t mutex;
    pthread_cond_t cond;

} wqueue_t;

typedef struct {
    rbtree_t* t;
    pthread_mutex_t mutex;
    pthread_cond_t cond;

} wqueue2_t;


wqueue_t* wqueue_alloc();
void wqueue_free(wqueue_t* wq, void (*payload_free)(void*));
int wqueue_full(wqueue_t* wq);
int wqueue_empty(wqueue_t* wq);

// non-blocking
void wqueue_prepend(wqueue_t* wq, void* item);
void wqueue_append(wqueue_t* wq, void* item);
void* wqueue_pop(wqueue_t* wq);

// blocking (respect max_length)
void wqueue_prepend_blocking(wqueue_t* wq, void* item);
void wqueue_append_blocking(wqueue_t* wq, void* item);
void* wqueue_pop_blocking(wqueue_t* wq);

int wqueue_length(wqueue_t* wq);
int wqueue_length_head(wqueue_t* wq);

void wqueue_wakeup(wqueue_t* wq);

void _lock(pthread_mutex_t* m);
void _unlock(pthread_mutex_t* m);
void _wait(pthread_cond_t* c, pthread_mutex_t* m);
int _timedwait(pthread_cond_t* c, pthread_mutex_t* m);

// prioritise by length
wqueue2_t* wqueue2_alloc();

void wqueue2_free(wqueue2_t* wq);

void wqueue2_push(wqueue2_t* wq, query_t* query);
query_t* wqueue2_pop(wqueue2_t* wq, int blocking);
int wqueue2_length(wqueue2_t* wq);
int wqueue2_empty(wqueue2_t* wq);
void wqueue2_wakeup(wqueue2_t* wq);

#endif

