
#include "tree.h"

void init_tree(tree *T , void *data_root) {

    tree_node *root = (tree_node*)malloc(sizeof(tree_node));

    root->data = data_root;
    root->child1 = NULL;
    root->child2 = NULL;
    root->depth = 0;

    T->root = root;
    T->tree_depth = 0;
    
}

/* ********************************************************************************************************** */

void destroy_tree(tree *tree) {
    destroy_subtree(tree->root);
}

/* ********************************************************************************************************** */

void destroy_subtree(tree_node *root) {
    if(root->child1 != NULL) destroy_subtree(root->child1);
    if(root->child2 != NULL) destroy_subtree(root->child2);
    free(root->data);
    free(root);
    
}

/* ********************************************************************************************************** */

void add_node_tree(tree *tree , tree_node *parent, void *data_node) {
    if(parent->child1 != NULL && parent->child2 != NULL) {
        perror("Try to add a child but no space available for this parent");
        return;
    }
    tree_node *new_child;
    new_child = (tree_node*)malloc(sizeof(tree_node));

    new_child->data = data_node;
    new_child->child1 = NULL;
    new_child->child2 = NULL;
    new_child->depth = parent->depth + 1;

    if(new_child->depth > tree->tree_depth) tree->tree_depth = new_child->depth;

    if(parent->child1 != NULL) parent->child1 = new_child;
    else parent->child2 = new_child;
}

/* ********************************************************************************************************** */

void print_tree_node(tree_node *node) {
    if(node->depth == 0) printf("ROOT\n");
    printf("Address data : %p\n",node->data);
    printf("Nb childs : %d\n", (node->child1 != NULL) + (node->child2 != NULL));
}

/* ********************************************************************************************************** */


void print_subtree(tree_node *root) {
    print_tree_node(root);
    if(root->child1 != NULL) print_subtree(root->child1);
    if(root->child2 != NULL) print_subtree(root->child2);

}

/* ********************************************************************************************************** */

void print_tree(tree *tree) {
    print_subtree(tree->root);
}

