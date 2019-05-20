
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

    for(int i=index+1 ; i < s.precision ; i++) {
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

void reduce_serie_ff(serie s , long p) {
    for(int i = 0 ; i < s.precision ; i++) {
        s.coeffs[i] %= p;
        while(s.coeffs[i] < 0) s.coeffs[i] += p;
    }
}

void add_serie(serie *s , serie a , serie b) {
    change_precision_serie(s , MAX(a.precision,b.precision));
    for(int i = 0 ; i < s->precision ; i++) 
        s->coeffs[i] = (i < a.precision ? a.coeffs[i] : 0) + (i < b.precision ? b.coeffs[i] : 0);
}

void add_serie_ff(serie *s , serie a , serie b , long p) {
    add_serie(s , a , b);
    reduce_serie_ff(*s , p);
}

void mul_serie(serie *s , serie a , serie b) {
    change_precision_serie(s , a.precision + b.precision - 1); // (a.precision - 1)  + (b.precision - 1) + 1
    memset(s->coeffs , 0 , s->precision * sizeof(s_coeff_t));
    
    for(int i = 0 ; i < s->precision ; i++) {
        for(int j = 0 ; j <= i ; j++) {
            s->coeffs[i] += (j < a.precision ? a.coeffs[j] : 0) * (i-j < b.precision ? b.coeffs[i-j] : 0);
        }
    }
}

void mul_serie_ff(serie *s , serie a , serie b , long p) {
    mul_serie(s , a , b);
    reduce_serie_ff(*s , p);
}

void inv_serie_ff(serie *inv , serie *s , unsigned int precision , long p) {
    if(!s->coeffs[0]) return;

    if(precision == 1) {
        change_precision_serie(inv , 1);
        inv->coeffs[0] = modular_inverse(s->coeffs[0] , p);
        return;
    }

    serie t,s2;
    int n = ceil(precision/2);
    init_serie(&t , n);
    init_serie(&s2 , n);
    memcpy(s2.coeffs , s->coeffs , n*sizeof(s_coeff_t));

    
    inv_serie_ff(&t , &s2 , n , p);
    
    serie temp;
    init_serie(&temp , precision);

    // temp = 2 - T*S
    for(int i = 0 ; i < temp.precision ; i++) {
        for(int j = 0 ; j <= i ; j++) {
            temp.coeffs[i] += -1 * (j < t.precision ? t.coeffs[j] : 0) * (i-j < s->precision ? s->coeffs[i-j] : 0);
            if(i == 0) temp.coeffs[i] += 2;
        }
    }
    
    change_precision_serie(inv , precision);

    for(int i = 0 ; i < precision ; i++) {
        inv->coeffs[i] = 0;
        for(int j = 0 ; j <= i ; j++) {
            inv->coeffs[i] += (j < t.precision ? t.coeffs[j] : 0) * (i-j < temp.precision ? temp.coeffs[i-j] : 0);
        }
        
    }
    reduce_serie_ff(*inv , p);

    destroy_serie(t);
    destroy_serie(temp);
    destroy_serie(s2);
    

    
}
