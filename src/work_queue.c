#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>

#include "work_queue.h"
#include "query.h"


wqueue_t* wqueue_alloc(int max_length) {
    wqueue_t* wq;
    int error_num;

    wq = (wqueue_t*) malloc(sizeof(wqueue_t));

    wq->head = NULL;
    wq->tail = NULL;
    wq->length = 0;
    wq->max_length = max_length;

    if((error_num = pthread_mutex_init(&wq->mutex, NULL)) != 0) {
        fprintf(stderr, "Error: pthread_mutex_init failed: %s\n", strerror(error_num));
        _exit(EXIT_FAILURE);
    }

    if((error_num = pthread_cond_init(&wq->cond, NULL)) != 0) {
        fprintf(stderr, "Error: pthread_cond_init failed: %s\n", strerror(error_num));
        _exit(EXIT_FAILURE);
    }

    return wq;
}

wqueue2_t* wqueue2_alloc() {
    wqueue2_t* wq;
    int error_num;

    wq = (wqueue2_t*) malloc(sizeof(wqueue2_t));
    wq->t = rb_alloc();

    if((error_num = pthread_mutex_init(&wq->mutex, NULL)) != 0) {
        fprintf(stderr, "Error: pthread_mutex_init failed: %s\n", strerror(error_num));
        _exit(EXIT_FAILURE);
    }

    if((error_num = pthread_cond_init(&wq->cond, NULL)) != 0) {
        fprintf(stderr, "Error: pthread_cond_init failed: %s\n", strerror(error_num));
        _exit(EXIT_FAILURE);
    }

    return wq;
}

void wqueue_free(wqueue_t* wq, void (*payload_free)(void*)) {
    witem_t* tmp;
    int error_num;

    if(wq) {
        while(wq->length) {
            tmp = wqueue_pop(wq);

            if(payload_free) {
                payload_free(tmp);
            }
        }

        if((error_num = pthread_mutex_destroy(&wq->mutex)) != 0) {
            fprintf(stderr, "Error: pthread_mutex_destroy failed: %s\n", strerror(error_num));
            _exit(EXIT_FAILURE);
        }

        if((error_num = pthread_cond_destroy(&wq->cond)) != 0) {
            fprintf(stderr, "Error: pthread_cond_destroy failed: %s\n", strerror(error_num));
            _exit(EXIT_FAILURE);
        }

        free(wq);
    }
}

void _query_deallocator(void *q) {
    query_free((query_t*) q);
}

void wqueue2_free(wqueue2_t* wq) {
    int error_num;

    if(wq) {
        if(wq->t) {
            rb_free(wq->t, _query_deallocator);
        }

        if((error_num = pthread_mutex_destroy(&wq->mutex)) != 0) {
            fprintf(stderr, "Error: pthread_mutex_destroy failed: %s\n", strerror(error_num));
            _exit(EXIT_FAILURE);
        }

        if((error_num = pthread_cond_destroy(&wq->cond)) != 0) {
            fprintf(stderr, "Error: pthread_cond_destroy failed: %s\n", strerror(error_num));
            _exit(EXIT_FAILURE);
        }

        free(wq);
    }
}

void _lock(pthread_mutex_t* m) {
    int error_num;

    if((error_num = pthread_mutex_lock(m)) != 0) {
        fprintf(stderr, "Error: pthread_mutex_lock failed: %s\n", strerror(error_num));
        _exit(EXIT_FAILURE);
    }
}

void _unlock(pthread_mutex_t* m) {
    int error_num;

    if((error_num = pthread_mutex_unlock(m)) != 0) {
        fprintf(stderr, "Error: pthread_mutex_unlock failed: %s\n", strerror(error_num));
        _exit(EXIT_FAILURE);
    }
}

void _wait(pthread_cond_t* c, pthread_mutex_t* m) {
    int error_num;

    if((error_num = pthread_cond_wait(c, m)) != 0) {
        fprintf(stderr, "Error: pthread_cond_wait failed: %s\n", strerror(error_num));
        _exit(EXIT_FAILURE);
    }
}

int _timedwait(pthread_cond_t* c, pthread_mutex_t* m) {
//    int error_num;
    struct timespec ts;

    clock_gettime(CLOCK_REALTIME, &ts);
    ts.tv_sec += 1;

    return pthread_cond_timedwait(c, m, &ts);
/*
    if((error_num = pthread_cond_timedwait(c, m, &ts)) != 0) {
        fprintf(stderr, "Error: pthread_cond_timedwait failed: %s\n", strerror(error_num));
        _exit(EXIT_FAILURE);
    }
*/
}

void _signal(pthread_cond_t* c) {
    int error_num;

    if((error_num = pthread_cond_signal(c)) != 0) {
        fprintf(stderr, "Error: pthread_cond_signal failed: %s\n", strerror(error_num));
        _exit(EXIT_FAILURE);
    }
}

int wqueue_full(wqueue_t* wq) {
    return wq->length >= wq->max_length;
}

int wqueue_empty(wqueue_t* wq) {
    return wq->length == 0;
}

int wqueue_length(wqueue_t* wq) {
    return wq->length;
}

int wqueue_length_head(wqueue_t* wq) {
    int count = 0;
    witem_t* tmp = wq->head;

    if(tmp == NULL)
        return 0;

    while((tmp = tmp->next) != NULL) {
        ++count;
    }

    return count;
}

void _prepend(wqueue_t* wq, void* item, int blocking) {
    witem_t* tmp;

    tmp = (witem_t*) malloc(sizeof(witem_t));
    tmp->item = item;
    tmp->prev = NULL;

    _lock(&wq->mutex);

    while(blocking && wqueue_full(wq)) {
        _wait(&wq->cond, &wq->mutex);
    }

    tmp->next = wq->head;

    if(!wq->head) {
        wq->head = wq->tail = tmp;
    }
    else {
        wq->head->prev = tmp;
        wq->head = tmp;
    }

    wq->length++;

    _unlock(&wq->mutex);
    _signal(&wq->cond);
}

void wqueue_prepend_blocking(wqueue_t* wq, void* item) {
    _prepend(wq, item, 1);
}

void wqueue_prepend(wqueue_t* wq, void* item) {
    _prepend(wq, item, 0);
}

void _append(wqueue_t* wq, void* item, int blocking) {
    witem_t* tmp;

    tmp = (witem_t*) malloc(sizeof(witem_t));
    tmp->item = item;
    tmp->next = NULL;

    _lock(&wq->mutex);

    while(blocking && wqueue_full(wq)) {
        _wait(&wq->cond, &wq->mutex);
    }

    tmp->prev = wq->tail;

    if(!wq->head) {
        wq->head = wq->tail = tmp;
    }
    else {
        wq->tail->next = tmp;
        wq->tail = tmp;
    }

    wq->length++;

    _unlock(&wq->mutex);
    _signal(&wq->cond);
}

void wqueue_append_blocking(wqueue_t* wq, void* item) {
    _append(wq, item, 1);
}

void wqueue_append(wqueue_t* wq, void* item) {
    _append(wq, item, 0);
}

void* _pop(wqueue_t* wq, int blocking) {
    witem_t* wi = NULL;
    void* tmp = NULL;

    _lock(&wq->mutex);

    while(blocking && wqueue_empty(wq)) {
        if(_timedwait(&wq->cond, &wq->mutex) != 0) {
            _unlock(&wq->mutex);
            return NULL;
        }
    }

    if(wq->tail) {
        //fprintf(stderr, "pop: tail is not NULL\n");

        wi = wq->tail;
        wq->tail = wi->prev;

        if(wq->tail) {
            wq->tail->next = NULL;
        }
        else {
            wq->head = NULL;
        }

        wq->length--;
    }
    //else {
    //    fprintf(stderr, "pop: tail is NULL\n");
    //}

    //fprintf(stderr, "wqueue: %d\n", wq->length);

    _unlock(&wq->mutex);
    _signal(&wq->cond);

    if(wi) {
        tmp = wi->item;
        free(wi);
    }

    return tmp;
}

void wqueue_wakeup(wqueue_t* wq) {
    _signal(&wq->cond);
}

void* wqueue_pop_blocking(wqueue_t* wq) {
    return _pop(wq, 1);
}

void* wqueue_pop(wqueue_t* wq) {
    return _pop(wq, 0);
}

void wqueue2_push(wqueue2_t* wq, query_t* query) {
    //fprintf(stderr, "push\n");
    _lock(&wq->mutex);
    rb_put(wq->t, seq_len(query->query), (void*) query);
    _unlock(&wq->mutex);
    _signal(&wq->cond);
}

int wqueue2_length(wqueue2_t* wq) {
    return rb_size(wq->t);
}

int wqueue2_empty(wqueue2_t* wq) {
    return rb_size(wq->t) == 0;
}

query_t* wqueue2_pop(wqueue2_t* wq, int blocking) {
    query_t* q = NULL;

    //fprintf(stderr, "pop\n");
    _lock(&wq->mutex);

    while(blocking && wqueue2_empty(wq)) {
        if(_timedwait(&wq->cond, &wq->mutex) != 0) {
            _unlock(&wq->mutex);
            return NULL;
        }
    }

    q = (query_t*) rb_popmax(wq->t);

    _unlock(&wq->mutex);
    _signal(&wq->cond);

    return q;
}

void wqueue2_wakeup(wqueue2_t* wq) {
    _signal(&wq->cond);
}

/*
int main() {
    wqueue_t *wq;
    void *tmp;
    int i = 1;
    int j = 2;
    int k = 3;

    wq = wqueue_alloc(5);

    wqueue_append(wq, (void*) &i);
    wqueue_pop(wq);
    wqueue_append(wq, (void*) &j);
    wqueue_append(wq, (void*) &k);

    while((tmp = wqueue_pop(wq)) != NULL) {
        fprintf(stderr, "popped %d\n", *((int*)tmp));
    }

    wqueue_free(wq, NULL);

    return 0;
}
*/

