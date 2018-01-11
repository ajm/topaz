#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <inttypes.h>

#include "rbtree.h"
#include "vector.h"


void _left_rotate(rbtree_t *t, rbnode_t *x) {
    rbnode_t *y;

    y = x->right;
    x->right = y->left;

    if(y->left != t->nil)
        y->left->parent = x;

    y->parent = x->parent;

    if(x->parent == t->nil) {
        t->root = y;
    }
    else {
        if(x->parent->left == x) {
            x->parent->left = y;
        }
        else {
            x->parent->right = y;
        }
    }

    y->left = x;
    x->parent = y;
}

void _right_rotate(rbtree_t *t, rbnode_t *x) {
    rbnode_t *y;

    y = x->left;
    x->left = y->right;

    if(y->right != t->nil)
        y->right->parent = x;

    y->parent = x->parent;

    if(x->parent == t->nil) {
        t->root = y;
    }   
    else {   
        if(x->parent->right == x) {
            x->parent->right = y; 
        }
        else {
            x->parent->left = y;
        }
    }
    
    y->right = x;
    x->parent = y;
}

void rb_insert(rbtree_t *t, rbnode_t *z) {
    rbnode_t *x, *y;

    y = t->nil;
    x = t->root;

    while(x != t->nil) {
        y = x;
        if(z->key < x->key) {
            x = x->left;
        }
        else {
            x = x->right;
        }
    }

    z->parent = y;

    if(y == t->nil) {
        t->root = z;
    }
    else {
        if(z->key < y->key) {
            y->left = z;
        }
        else {
            y->right = z;
        }
    }

    z->left = t->nil;
    z->right = t->nil;
    z->colour = RED;

    // rb_insert_fixup
    while(z->parent->colour == RED) {
        if(z->parent == z->parent->parent->left) {
            y = z->parent->parent->right;
            if(y->colour == RED) {
                z->parent->colour = BLACK;
                y->colour = BLACK;
                z->parent->parent->colour = RED;
                z = z->parent->parent;
            }
            else {
                if(z == z->parent->right) {
                    z = z->parent;
                    _left_rotate(t, z);
                }

                z->parent->colour = BLACK;
                z->parent->parent->colour = RED;
                
                _right_rotate(t, z->parent->parent);
            }
        }
        else {
            // same as before, but right <--> left
            y = z->parent->parent->left;
            if(y->colour == RED) {
                z->parent->colour = BLACK;
                y->colour = BLACK;
                z->parent->parent->colour = RED;
                z = z->parent->parent;
            }
            else {
                if(z == z->parent->left) {
                    z = z->parent;
                    _right_rotate(t, z);
                }

                z->parent->colour = BLACK;
                z->parent->parent->colour = RED;
                
                _left_rotate(t, z->parent->parent);
            }
        }
    }

    t->root->colour = BLACK;
    t->count++;
}

rbnode_t* rb_minimum(rbtree_t* t, rbnode_t* z) {
    rbnode_t *x = z;

    while(x->left != t->nil) {
        x = x->left;
    }

    return x;
}

rbnode_t* rb_maximum(rbtree_t* t, rbnode_t* z) {
    rbnode_t *x = z;

    while(x->right != t->nil) {
        x = x->right;
    }

    return x;
}

rbnode_t* rb_successor(rbtree_t* t, rbnode_t *z) {
    rbnode_t *x, *y;

    x = z;

    if(x->right != t->nil) {
        return rb_minimum(t, x->right);
    }

    y = x->parent;

    while(y != t->nil && x == y->right) {
        x = y;
        y = y->parent;
    }

    return y;
}

void rb_transplant(rbtree_t *t, rbnode_t *u, rbnode_t *v) {
    if(u->parent == t->nil) {
        t->root = v;
    }
    else {
        if(u == u->parent->left) {
            u->parent->left = v;
        }
        else {
            u->parent->right = v;
        }
    }
    v->parent = u->parent;
}

rbnode_t* rb_delete(rbtree_t *t, rbnode_t *z) {
    rbnode_t *w, *x, *y;
    enum colour y_original_colour;

    y = z;
    y_original_colour = y->colour;

    if(z->left == t->nil) {
        //fprintf(stderr, "case 1\n");
        x = z->right;
        rb_transplant(t, z, z->right);
    }
    else if(z->right == t->nil) {
        //fprintf(stderr, "case 2\n");
        x = z->left;
        rb_transplant(t, z, z->left);
    }
    else {
        //fprintf(stderr, "case 3\n");
        y = rb_minimum(t, z->right);
        y_original_colour = y->colour;
        x = y->right;
        
        if(y->parent == z) {
            x->parent = y;
        }
        else {
            rb_transplant(t, y, y->right);
            y->right = z->right;
            y->right->parent = y;
        }

        rb_transplant(t, z, y);
        y->left = z->left;
        y->left->parent = y;
        y->colour = z->colour;
    }

    if(y_original_colour == BLACK) {

        //fprintf(stderr, "FIX THE TREE\n");

        while(x != t->root && x->colour == BLACK) {
            if(x == x->parent->left) {
                w = x->parent->right;

                assert(w != t->nil);

                if(w->colour == RED) {
                    w->colour = BLACK;
                    x->parent->colour = RED;
                    _left_rotate(t, x->parent);
                    w = x->parent->right;
                }

                assert(w != t->nil);

                if(w->left->colour == BLACK && w->right->colour == BLACK) {
                    w->colour = RED;
                    x = x->parent;
                }
                else {
                    if(w->right->colour == BLACK) {
                        w->left->colour = BLACK;
                        w->colour = RED;
                        _right_rotate(t, w);
                        w = x->parent->right;
                    }

                    w->colour = x->parent->colour;
                    x->parent->colour = BLACK;
                    w->right->colour = BLACK;
                    _left_rotate(t, x->parent);
                    x = t->root;
                }
            }
            else {
                w = x->parent->left;

                assert(w != t->nil);

                if(w->colour == RED) {
                    w->colour = BLACK;
                    x->parent->colour = RED;
                    _right_rotate(t, x->parent);
                    w = x->parent->left;
                }

                assert(w != t->nil);

                if(w->left->colour == BLACK && w->right->colour == BLACK) {
                    w->colour = RED;
                    x = x->parent;
                }
                else {
                    if(w->left->colour == BLACK) {
                        w->right->colour = BLACK;
                        w->colour = RED;
                        _left_rotate(t, w);
                        w = x->parent->left;
                    }

                    w->colour = x->parent->colour;
                    x->parent->colour = BLACK;
                    w->left->colour = BLACK;
                    _right_rotate(t, x->parent);
                    x = t->root;
                }
            }
        }

        x->colour = BLACK;
    }

    t->count--;

    // ensure it can be safely freed
    z->parent = t->nil;
    z->left = t->nil;
    z->right = t->nil;

    return z;
}

rbnode_t *rb_search(rbtree_t *t, uint64_t key) {
    rbnode_t *x;

    x = t->root;

    while(x != t->nil && key != x->key) {
        if(key < x->key) {
            x = x->left;
        }
        else {
            x = x->right;
        }
    }

    return x;
}

rbtree_t *rb_alloc() {
    rbtree_t *t;

    t = (rbtree_t*) malloc(sizeof(rbtree_t));

    t->count = 0;
    t->nil = (rbnode_t*) malloc(sizeof(rbnode_t));

    t->nil->colour = BLACK;
    t->nil->parent = NULL;
    t->nil->left = NULL;
    t->nil->right = NULL;
    t->nil->payload = NULL;

    t->root = t->nil;

    return t;
}

void rb_node_free(rbtree_t *t, rbnode_t *x, void (*payload_free)(void* pl)) {

    if(x != t->nil) {

        if(!x->left) fprintf(stderr, "left is NULL\n");

        if(x->left != t->nil) {
            rb_node_free(t, x->left, payload_free);
        }

        if(!x->right) fprintf(stderr, "right is NULL\n");

        if(x->right != t->nil) {
            rb_node_free(t, x->right, payload_free);
        }

        if((payload_free != NULL) && x->payload) {
            payload_free(x->payload);
        }

        free(x);
    }
}

void rb_free(rbtree_t *t, void (*payload_free)(void* pl)) {

    if(t) {
        assert(((t->count == 0) && (t->root == t->nil)) \
            || ((t->count != 0) && (t->root != t->nil)));

        if(t->root != t->nil) {
            rb_node_free(t, t->root, payload_free);
        }

        free(t->nil);
        free(t);
    }
}

void* rb_get(rbtree_t *t, uint64_t key) {
    rbnode_t *tmp;
    
    if((tmp = rb_search(t, key)) != t->nil) {
        return tmp->payload;
    }

    return NULL;
}

void* rb_popmax(rbtree_t *t) {
    rbnode_t *tmp;
    void *pl;
    int before = rb_size(t);

    if(rb_empty(t))
        return NULL;

    tmp = rb_maximum(t, t->root);   // find biggest key
    tmp = rb_delete(t, tmp);        // delete rb node
    pl = tmp->payload;              // copy pointer to the payload
    tmp->payload = NULL;            // make the original NULL so it cannot be freed
    rb_node_free(t, tmp, NULL);     // free the rb_node, but not the payload

    assert(before - 1 == rb_size(t));

    return pl;
}

void* rb_popmin(rbtree_t *t) {
    rbnode_t *tmp;
    void *pl;
    int before = rb_size(t);

    if(rb_empty(t))
        return NULL;

    tmp = rb_minimum(t, t->root);   // find smallest key
    tmp = rb_delete(t, tmp);        // delete rb node
    pl = tmp->payload;              // copy pointer to the payload
    tmp->payload = NULL;            // make the original NULL so it cannot be freed
    rb_node_free(t, tmp, NULL);     // free the rb_node, but not the payload

    assert(before - 1 == rb_size(t));

    return pl;
}

void* rb_getmin(rbtree_t *t) {
    rbnode_t *tmp;

    if(rb_empty(t))
        return NULL;

    tmp = rb_minimum(t, t->root);

    return tmp->payload;
}

void rb_truncate(rbtree_t *t, int size, void (*payload_free)(void* pl)) {
    rbnode_t *tmp;

    while(rb_size(t) > size) {
        tmp = rb_minimum(t, t->root);       // find smallest key
        tmp = rb_delete(t, tmp);            // delete rb node from the tree
        rb_node_free(t, tmp, payload_free); // free the rb_node
    }
}

// delete by key, ignore payload, assumes caller already has pointer
void rb_rm(rbtree_t *t, uint64_t key) {
    rbnode_t *tmp;
   
    if((tmp = rb_search(t, key)) != t->nil) {
        tmp = rb_delete(t, tmp);
        rb_node_free(t, tmp, NULL);
    }
}

rbnode_t* rb_put(rbtree_t *t, uint64_t key, void *payload) {
    rbnode_t *tmp;

    tmp = (rbnode_t *) malloc(sizeof(rbnode_t));
        
    tmp->key = key;
    tmp->payload = payload;

    tmp->left = t->nil;
    tmp->right = t->nil;
    tmp->parent = t->nil;

    rb_insert(t, tmp);
    
    return tmp;
}

rbnode_t* rb_reput(rbtree_t *t, uint64_t key, rbnode_t *tmp) {

    tmp = rb_delete(t, tmp);

    tmp->key = key;

    rb_insert(t, tmp);

    return tmp;
}

int _compare_int(const void *l, const void *r) {
    int li = *((int*) l);
    int ri = *((int*) r);

    return li - ri;
}

int rb_nth_score(rbtree_t *t, int num_items, int (*get_value)(void*)) {
    int *scores;
    int scores_ptr;
    vector_t *stack;
    int nth_value;
    rbnode_t *tmp;

    scores = (int*) malloc(sizeof(int) * t->count);
    scores_ptr = 0;

    stack = vector_alloc();
    vector_push(stack, t->root);


    while((tmp = (rbnode_t*) vector_pop(stack)) != NULL) {
        scores[scores_ptr++] = get_value(tmp->payload);

        if(tmp->left != t->nil)
            vector_push(stack, tmp->left);

        if(tmp->right != t->nil)
            vector_push(stack, tmp->right);
    }

    qsort(scores, t->count, sizeof(int), _compare_int);
    
    nth_value = scores[t->count - num_items];

    vector_free(stack);
    free(scores);

    return nth_value;
}

void rb_prune(rbtree_t *t, int threshold, int (*prune_test)(void*, int), void (*payload_free)(void*)) {
    rbnode_t *tmp; //, *copy;
    vector_t *stack;
    vector_t *rmlist;
    
    rmlist = vector_alloc();
    stack = vector_alloc();
    vector_push(stack, t->root);


    while((tmp = (rbnode_t*) vector_pop(stack)) != NULL) {

        if(tmp->left != t->nil)
            vector_push(stack, tmp->left);

        if(tmp->right != t->nil)
            vector_push(stack, tmp->right);

        if(prune_test(tmp->payload, threshold))
            vector_push(rmlist, tmp);
    }

    while((tmp = vector_pop(rmlist)) != NULL) {
        //copy = tmp;        

        tmp = rb_delete(t, tmp);

        rb_node_free(t, tmp, payload_free); // XXX
    }

    vector_free(stack);
    vector_free(rmlist);
}

int rb_size(rbtree_t* t) {
    return t->count;
}

// used for debugging
int rb_size_traverse(rbtree_t* t) {
    traversal_t* tr;
    int count = 0;

    tr = traversal_alloc(t);

    while(traversal_next(tr) != NULL) {
        ++count;
    }

    traversal_free(tr);
    return count;
}

int rb_empty(rbtree_t* t) {
    return t->count == 0;
}

traversal_t* traversal_alloc(rbtree_t* t) {
    traversal_t* tr;

    tr = (traversal_t*) malloc(sizeof(traversal_t));
    tr->stack = vector_alloc();
    tr->tree = t;
    tr->node = t->nil;

    //if(tr->tree->root != tr->tree->nil)
    //    vector_push(tr->stack, tr->tree->root);
    
    return tr;
}

/*

// post-order traversal
void traversal_start(traversal_t* tr) {
    vector_reset(tr->stack);

    if(tr->tree->root != tr->tree->nil)
        vector_push(tr->stack, tr->tree->root);
}

void* traversal_next(traversal_t* tr) {
    rbnode_t* tmp = (rbnode_t*) vector_pop(tr->stack);

    if(!tmp)
        return NULL;

    if(tmp->left != tr->tree->nil)
        vector_push(tr->stack, tmp->left);

    if(tmp->right != tr->tree->nil)
        vector_push(tr->stack, tmp->right);

    return tmp->payload;
}
*/

// inorder traversal
void traversal_start(traversal_t* tr) {
    vector_reset(tr->stack);

    if(tr->tree->root != tr->tree->nil)
        tr->node = tr->tree->root;
}

void* traversal_next(traversal_t* tr) {
    rbnode_t *tmp;
    
    while(!vector_empty(tr->stack) || (tr->node != tr->tree->nil)) {
        if(tr->node != tr->tree->nil) {
            vector_push(tr->stack, tr->node);
            tr->node = tr->node->right;
        }
        else {
            tmp = (rbnode_t*) vector_pop(tr->stack);
            tr->node = tmp->left;
            return tmp->payload;
        }
    }

    return NULL;
}

void traversal_end(traversal_t* tr) {
    vector_reset(tr->stack);
}

void traversal_free(traversal_t* tr) {
    if(tr) {
        vector_free(tr->stack);
        free(tr);
    }
}

void rb_print(rbtree_t *t) {
    rbnode_t* tmp;
    vector_t* stack;

    if(t->root == t->nil) {
        return;
    }

    stack = vector_alloc();
    vector_push(stack, t->root);

    while((tmp = (rbnode_t*) vector_pop(stack)) != NULL) {
        fprintf(stderr, "(%" PRIu64 ", %c) -> (%" PRIu64 ", %c) (%" PRIu64 ", %c), P=%c, root=%d\n", 
                                                tmp->key, tmp->colour == BLACK ? 'B' : 'R',
                                                tmp->left->key, tmp->left->colour == BLACK ? 'B' : 'R',
                                                tmp->right->key, tmp->right->colour == BLACK ? 'B' : 'R',
                                                tmp->parent->colour == BLACK ? 'B' : 'R',
                                                tmp == t->root );
        if(tmp->left != t->nil)
            vector_push(stack, tmp->left);

        if(tmp->right != t->nil)
            vector_push(stack, tmp->right);
    }

    vector_free(stack);
}

int rb_check(rbtree_t *t) {
    rbnode_t* tmp;
    vector_t* stack;
    int ret = 1;

    if(t->root == t->nil)
        return 1;

    if(t->root->colour != BLACK)
        return 0;

    stack = vector_alloc();
    vector_push(stack, t->root);

    while((tmp = (rbnode_t*) vector_pop(stack)) != NULL) {
        if(tmp->left == NULL) {
            fprintf(stderr, "(%" PRIu64 ") left is NULL!\n", tmp->key);
            return 0;
        }

        if(tmp->right == NULL) {
            fprintf(stderr, "(%" PRIu64 ") right is NULL!\n", tmp->key);
            return 0;
        }

        if(tmp->parent == NULL) {
            fprintf(stderr, "(%" PRIu64 ") parent is NULL!\n", tmp->key);
            return 0;
        }

        if(tmp->colour == RED) {
            ret &= tmp->left->colour == BLACK;
            ret &= tmp->right->colour == BLACK;
            ret &= tmp->parent->colour == BLACK;
        }

        if(!ret) {
            fprintf(stderr, "[BAD] (%" PRIu64 ") %c -> (%c, %c), p=%c, root=%d\n", tmp->key, tmp->colour == BLACK ? 'B' : 'R',
                                                tmp->left->colour == BLACK ? 'B' : 'R',
                                                tmp->right->colour == BLACK ? 'B' : 'R',
                                                tmp->parent->colour == BLACK ? 'B' : 'R',
                                                tmp == t->root );
            break;
        }

        if(tmp->left != t->nil)
            vector_push(stack, tmp->left);

        if(tmp->right != t->nil)
            vector_push(stack, tmp->right);
    }

    vector_free(stack);
    return ret;
}

