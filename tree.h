
#ifndef TREE_H
#define TREE_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

typedef struct tree_node {
    void *data;
    struct tree_node *child1;
    struct tree_node *child2;
    unsigned int depth;
} tree_node;

typedef struct {
    tree_node *root;
    unsigned int tree_depth;
} tree; 

void init_tree(tree *T , void *data_root);

void destroy_tree(tree *T);

void destroy_subtree(tree_node *root);

void add_node_tree(tree *T , tree_node *parent, void *data_node);

void print_tree_node(tree_node *node);

void print_tree(tree *tree);

void print_subtree(tree_node *root);



#endif