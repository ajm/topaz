#ifndef H_RBTREE
#define H_RBTREE

#include <sys/types.h>
#include <stdint.h>

#include "vector.h"


enum colour { RED, BLACK };

typedef struct rbnode {
    uint64_t key;
    void *payload;

    struct rbnode *left, *right, *parent;
    enum colour colour;

} rbnode_t;

typedef struct {
    rbnode_t *root;
    rbnode_t *nil;
    int count;

} rbtree_t;

typedef struct {
    rbtree_t *tree;
    vector_t *stack;
    rbnode_t *node;

} traversal_t;


rbtree_t *rb_alloc(); // public
void rb_free(rbtree_t *t, void (*payload_free)(void*)); // public
void rb_node_free(rbtree_t *t, rbnode_t *x, void (*payload_free)(void*));
void* rb_get(rbtree_t *t, uint64_t key); // public
rbnode_t* rb_put(rbtree_t *t, uint64_t key, void *payload); // public
rbnode_t* rb_reput(rbtree_t *t, uint64_t key, rbnode_t *tmp); // public
void rb_prune(rbtree_t *t, int threshold, int (*prune_test)(void*, int), void (*payload_free)(void*)); // public
int rb_nth_score(rbtree_t *t, int num_items, int (*get_value)(void*)); // public
int rb_size(rbtree_t* t); // public
int rb_size_traverse(rbtree_t* t);
int rb_empty(rbtree_t* t); // public
void* rb_popmax(rbtree_t* t);
void* rb_popmin(rbtree_t* t); // public
void* rb_getmin(rbtree_t *t);
void rb_truncate(rbtree_t *t, int size, void (*payload_free)(void* pl)); // public
void rb_rm(rbtree_t *t, uint64_t key); // public

void _left_rotate(rbtree_t *t, rbnode_t *x);
void _right_rotate(rbtree_t *t, rbnode_t *x);
void rb_insert(rbtree_t *t, rbnode_t *z);
rbnode_t* rb_delete(rbtree_t *t, rbnode_t *z);
rbnode_t* rb_search(rbtree_t *t, uint64_t key);
//void rb_increment(rbtree_t *t, uint64_t key);
//void rb_print2(rbtree_t *t);
int rb_check(rbtree_t *t);

traversal_t* traversal_alloc(rbtree_t* t);
void traversal_start(traversal_t* tr);
void* traversal_next(traversal_t* tr);
void traversal_end(traversal_t* tr);
void traversal_free(traversal_t* tr);

#endif

