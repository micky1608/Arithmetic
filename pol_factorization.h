/*
* Created on 2019-05-13
*/

#ifndef FACTORIZATION_H
#define FACTORIZATION_H

#include "pol.h"

typedef struct triplet {
    pol G;
    long q;
    long m;
} triplet;

typedef struct fact_list_element{
    triplet triplet;
    struct fact_list_element *next;
    struct fact_list_element *prev;
} fact_list_element;

typedef struct fact_list_info {
    int size;
    struct fact_list_element *first;
    struct fact_list_element *last;
} fact_list_info;


void init_fact_list(fact_list_info **info);

void add_triplet_fact_list_end(fact_list_info **info , triplet t);

void print_fact_list(fact_list_info *info);

void clean_fact_list(fact_list_info *info);

/* ************************************************************************ */

void fact_algo1(fact_list_info **L , pol *C , pol F , long P);

#endif