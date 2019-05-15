

#include "pol_factorization.h"



void init_fact_list(fact_list_info **info) {
    (*info) = (fact_list_info*)malloc(sizeof(fact_list_info));
    (*info)->size = 0;
    (*info)->first = (*info)->last = NULL;
}

/* ********************************************************************************************************************** */

void clean_fact_list(fact_list_info *info) {
    if(info->size == 0) return;
    fact_list_element *toRemove, *toRemoveNext;

    toRemove = info->first;

    while(toRemove != NULL) {
        toRemoveNext = toRemove->next;
        destroy_pol(toRemove->triplet.G);
        free(toRemove);
        toRemove = toRemoveNext;
    }
    free(info);
}

/* ********************************************************************************************************************** */

void add_triplet_fact_list_end(fact_list_info **info , triplet t) {
    fact_list_element *new_element = (fact_list_element *)malloc(sizeof(fact_list_element));
    new_element->triplet = t;
    new_element->next = NULL;
    new_element->prev = (*info)->last;

    if((*info)->size) 
        (*info)->last->next = new_element;
    else (*info)->first = new_element;

    (*info)->last = new_element;
    ((*info)->size)++;
}

/* ********************************************************************************************************************** */

void print_fact_list(fact_list_info *info) {
    if(info->size == 0) return;

    fact_list_element *iterator = info->first;
    while(iterator != NULL) {
        printf("Triplet : q = %ld\tm = %ld\t",iterator->triplet.q,iterator->triplet.m);
        print_pol(iterator->triplet.G , "G");
        iterator = iterator->next;
    }
}

/* ********************************************************************************************************************** */

void fact_algo1(fact_list_info **L , pol *C , pol F , long P) {
    
    init_fact_list(L);
   
    int i=1;

    pol F_der , G , trash , trash2, PPol;

    derivate_pol(&F_der , F);
    reduce_pol_ff(&F_der , P);

    gcd_pol_ff(C , F , F_der , P);

    euclide_div_pol_ff(&G , &trash , F , *C , P);


    while(G.degree >= 1) {
        gcd_pol_ff(&PPol , *C , G, P);
        euclide_div_pol_ff(&trash2 , &trash , *C , PPol , P);
        copy_pol(C , trash2);
    
        if(G.degree > PPol.degree) {
            triplet t;

            euclide_div_pol_ff(&t.G , &trash , G , PPol , P);
            t.q = 1;
            t.m = i;

            add_triplet_fact_list_end(L , t);
        }

        copy_pol(&G , PPol);
        i++;
    }
}

/* ********************************************************************************************************************** */

void fact_algo2(fact_list_info **L , pol F , long P) {
    
    pol C, F2;

    fact_algo1(L , &C , F , P);

    if(C.degree == 0 && C.coeffs[0] == 1) return;

    init_pol(&F2 , C.degree/P);
    for(int i = 0 ; i <= F2.degree ; i++) {
        F2.coeffs[i] = C.coeffs[i*P];
    }

    fact_list_info *l;
    fact_algo2(&l , F2 , P);

    fact_list_element *element = l->first;

    while(element != NULL) {
        element->triplet.q *= P;
        add_triplet_fact_list_end(L , element->triplet);
    }

    destroy_pol(C);
    destroy_pol(F2);

}