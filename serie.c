#include "serie.h"


void init_serie(serie *s, unsigned int precision) {
    if(precision > 0) {
        s->precision = precision;
        s->coeffs = (s_coeff_t*)malloc(sizeof(s_coeff_t) * (precision));
        memset(s->coeffs , 0 , precision*sizeof(s_coeff_t));
    }
}

void set_all_coeffs_serie(serie *s, s_coeff_t *coeffs, size_t size_coeffs_array) {
    s->precision = size_coeffs_array;
    for(int i = 0 ; i < s->precision ; i++) {
        s->coeffs[i] = coeffs[i];
    }
}

void destroy_serie(serie s) {
    free(s.coeffs);
}

void change_precision_serie(serie *s, unsigned int new_precision) {
    if(new_precision == s->precision) return;

     for(int i = 0 ; i < s->precision ; i++) {
        printf("%ld ",s->coeffs[i]);
    }
    printf("\n");
    
    s_coeff_t *r = realloc(s->coeffs , new_precision*sizeof(s_coeff_t));

    if(r) {
        s->coeffs = r;
        if(new_precision > s->precision) {
            memset(s->coeffs + s->precision*sizeof(s_coeff_t) , 0 , (new_precision - s->precision)*sizeof(s_coeff_t));
        }
        s->precision = new_precision;
    }

}


void print_serie(serie s , char *name) {
    if(name != NULL) printf("%s (in precision %u): ",name,s.precision);

    s_coeff_t temp;

    int index=0;
    while(index < s.precision && !s.coeffs[index]) ++index;

    if(index == s.precision) {
        printf("0\n");
        return;
    }

    if(index) printf("(");

    #ifdef S_COEFF_LONG
        printf("%ld",s.coeffs[index]);
    #endif

    #ifdef S_COEFF_DOUBLE
        printf("%.2f",s.coeffs[index]);
    #endif

    if(index) printf(" * X^%d)" , index);

    for(int i=index+1 ; i <= s.precision ; i++) {
        if (s.coeffs[i] < 0 )
            printf(" - ");
        else if (s.coeffs[i] > 0 )
            printf(" + ");
        else
            continue;

        temp = (s.coeffs[i] > 0) ? s.coeffs[i] : -1*s.coeffs[i];
        printf("(");

        #ifdef S_COEFF_LONG
            printf("%ld",temp);
        #endif

        #ifdef S_COEFF_DOUBLE
            printf("%.2f",temp);
        #endif
        
        printf(" * X^%d)" , i);

    }
    
    if(name != NULL) printf("\n");
}
