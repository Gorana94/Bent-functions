#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <stdbool.h>
#include <emmintrin.h> 
#include <pthread.h>

#ifdef _WIN32
#include <windows.h>
#elif MACOS
#include <sys/param.h>
#include <sys/sysctl.h>
#else
#include <unistd.h>
#endif

int getNumCores() {
#ifdef WIN32
    SYSTEM_INFO sysinfo;
    GetSystemInfo(&sysinfo);
    return sysinfo.dwNumberOfProcessors;
#elif MACOS
    int nm[2];


    size_t len = 4;
    uint32_t count;

    nm[0] = CTL_HW; nm[1] = HW_AVAILCPU;
    sysctl(nm, 2, &count, &len, NULL, 0);

    if(count < 1) {
        nm[1] = HW_NCPU;
        sysctl(nm, 2, &count, &len, NULL, 0);
        if(count < 1) { count = 1; }
    }
    return count;
#else
    return sysconf(_SC_NPROCESSORS_ONLN);
#endif
}

typedef unsigned uint;
typedef unsigned long uint_32;

/* Global variables */

uint g_bent_nonlinearity = 0;
uint g_higher_correct_TT_wt = 0;
uint g_lower_correct_TT_wt = 0;
uint g_correct_sub_wt = 0;



uint random_at_most(uint max) {
  uint
    // max <= RAND_MAX < ULONG_MAX, so this is okay.
    num_bins = (uint) max + 1,
    num_rand = (uint) RAND_MAX + 1,
    bin_size = num_rand / num_bins,
    defect   = num_rand % num_bins;

  uint x;
  do {
   x = rand();
  }
  // This is carefully written not to overflow
  while (num_rand - defect <= (uint_32)x);

  // Truncated division is intentional
  return x/bin_size;

}

void create_array_of_length(uint ANF_size, uint ** out_anf)
{
    *out_anf = (uint *)malloc(ANF_size * sizeof(uint));
}

uint scalar_product_unit_vector(long long unsigned x, uint unit_vector) {
    return x & unit_vector;
}


void split_along_direction_e_n(uint * f, uint* f0, uint * f1, uint n, uint ei)  /* TODO Split_along_direction */
{
    uint f0_counter = 0, f1_counter = 0;
    uint total_split = (uint)pow(2,n);
    // printf("there are only %u elements in array\n", total_split);
    for (long long unsigned x = 0; x < (uint)pow(2,n); x++) {
        if (scalar_product_unit_vector(x, ei) == 0) {
            // printf("updating f0 %u\n", f0_counter);
            f0[f0_counter] = f[x];
            f0_counter++;
        } else {
            // printf("updating f1 %u\n", f1_counter);
            f1[f1_counter] = f[x];
            f1_counter++;
        }
    }
}

uint wt(uint* TT, uint length) /* TODO WT */
{
    uint result = 0;

    // printf("TT unutar wt\n");
    //         for (int i = 0; i < length; i++) {
    //             printf("%u ", TT[i]);
    //         }
    //         printf("\n");

    for (uint i = 0; i < length; i++) {
        result += TT[i];
    }

    return result;
}

uint ord2(long long unsigned term_index, uint n)
{
    uint bits = sizeof(term_index) * 8;

    uint result = 0;
    for (int i = 0; i < n; i++){
        result += (1 & (term_index >> i));
        // printf("%d", (1 & (term_index >> i)));
    }
    // printf("\n");
    return result;
}

void generate_wht2(uint* TT, uint* WHT, uint * temp, uint length, uint n) {
    for (uint i=0; i<length; i++) {
        WHT[i]=(1 - 2*TT[i]);
    }

    for (uint i=0; i<n; i++) {
        for (uint j=0; j<length; j++)
            temp[j]=WHT[j];
        for (uint j=0; j<length; j++) {
            if (j<((int)pow(2,n-1))) {
                WHT[j]=temp[2*j]+temp[2*j+1];
            } else {
                WHT[j]=temp[2*j-(length)]-temp[2*j-(length)+1];
            }
        }
    }
}

void generate_wht(uint* tt, int* wht, int* temp, uint n)
{
    for (uint i=0; i<(int)pow(2,n); i++) {
        wht[i]=(1 - 2*tt[i]);
    }

    for (uint i=0; i<n; i++) {
        for (uint j=0; j<(int)pow(2,n); j=j+(int)pow(2,i+1)) {
            for (uint k=j; k<j+(int)pow(2,i); k++) {
                int a=wht[k]-wht[k+(int)pow(2,i)];
                wht[k]=wht[k]+wht[k+(int)pow(2,i)];
                wht[k+(int)pow(2,i)]=a;
            }
        }
    }

}

uint skalarni_proizvod(uint x, uint w, uint n)
{
    uint res = 0;
    uint sume = w & x;
    for (int i = 0; i < n; i++) {
        res ^= 0x1 & (sume >> i);
    }

    return res;
}

int rucni_wht_za_jedno_w(uint* polar_tt, uint length, uint n, uint w)
{
    int wht = 0;

    for (int x = 0; x < length; x++) {
        int lin_rez = 1 - 2*skalarni_proizvod(x, w, n);
        wht += polar_tt[x]*lin_rez;
    }

    return wht;
}

void rucni_wht(uint* TT, int* WHT, uint* polar_tt, uint length, uint n)
{
    for (int i = 0; i < length; i++) {
        polar_tt[i] = 1 - 2*TT[i];
    }

    for (int w = 0; w < length; w++){
        WHT[w] = rucni_wht_za_jedno_w(polar_tt, length, n, w);
    }
}

uint nonlinearity(uint* TT, int* WHT, uint* temp, uint length, uint n) {

    rucni_wht(TT, WHT, temp, length, n);
    uint max_wht = 0;
    int the_only_correct_wht_value_for_an_element_of_wht_of_a_bent_function = (int)pow(2, n/2);
    bool good_wht = true;
    for (uint i = 0; i < length; i++) {
        if (abs(WHT[i]) > max_wht) {
            max_wht = abs(WHT[i]);
            if(abs(WHT[i]) != the_only_correct_wht_value_for_an_element_of_wht_of_a_bent_function) {
                good_wht = false;
                // printf("we dont have a correct WHT value: %d at index %d\n", WHT[i], i);
            }
        }
    }

    return (length - max_wht) >> 1;
}


void anft_naive(uint * tt, uint* anf, uint n) {
    for (int i=0; i<(int)pow(2,n); i++) {
        anf[i]=tt[i];
    }
    for (int i=0; i<n; i++) {
        for (int j=0; j<(int)pow(2,n); j=j+(int)pow(2,i+1)) {
            for (int k=j; k<j+(int)pow(2,i); k++) {
                anf[k+(int)pow(2,i)]=(anf[k] + anf[k+(int)pow(2,i)])%2;
            }
        }
    }
}

void anft_naive_math(uint * anf, uint * tt, uint * temp, uint n)
{
    for (int i=0; i<(int)pow(2,n); i++) anf[i]=tt[i];
    int k1=2;
    int k2=1;
    for (int i=0; i<n; i++) {
        for (int j=0; j<(int)pow(2,n); j++)
            temp[j]=anf[j];
        for (int j=0; j<(int)pow(2,n); j++)
            if ((j%k1)<k2) { /*magicni uslov*/
                anf[j]=temp[j];
            } else {
                anf[j]=(temp[j] + temp[j-k2])%2;
            }
        k1=k1*2;
        k2=k2*2;
    }
}

int main(int argc, char** argv)
{
    uint n = 8;
    uint max_ord = 3;
    if (argc == 2) {
        n = atoi(argv[1]);
    } else if (argc == 3) {
        n = atoi(argv[1]);
        max_ord = atoi(argv[2]);
    }

    uint truth_table_length = (uint)(pow(2,n));
    uint *term_thresholds = NULL;
    time_t t;
    uint bent_counter = 0;
    uint ANF_positions = 0;
    uint total_positions = 0;
    uint half_balance_changes = 0;
    uint max_bent_awt = 0;
    FILE *f = NULL;
    f = fopen("bent_functions.txt", "wr+");
    if (f == NULL) {
        printf("ERROR WHILE OPENING FILE");
    }
    fclose(f);

    srand(time(0));


step_1: ; /*blank statement because of how C works*/                                /* 1. Create ANF array of size 2^n and initialize: A_i = 0, 0 <= i <= 2^n - 1*/
    uint *ANF;
    create_array_of_length(truth_table_length, &ANF);
    memset(ANF, 0, truth_table_length*sizeof(uint));
    /* Create all other arrays and ANFs */
    uint *ANF_tmp;
    create_array_of_length(truth_table_length, &ANF_tmp);

    uint* TT;
    create_array_of_length(truth_table_length, &TT);

    uint * f0, *f1;
    create_array_of_length(truth_table_length/2, &f0);
    create_array_of_length(truth_table_length/2, &f1);

    uint * WHT;
    create_array_of_length(truth_table_length, &WHT);

    uint * temp_WHT;
    create_array_of_length(truth_table_length, &temp_WHT);

step_2: /* Create a quadratic bent function */
    for (int i = 0; i <= n/2 - 1; i ++) {                                           /* 2. For i = 0,...,N/2-1      ORIGINALLY 2^(N/2-1) */
        int j = (uint)(pow(2,i)) + (uint)(pow(2, n/2+i));                           /*      (a) Let j = 2^i + 2^(n/2+i) */
        ANF[j] = 1;                                                                 /*      (b) Set A_j = 1 */
    }


step_3: ; /*blank statement because of how C works*/                                /* 3. Select a random r, 0<= r <= (2^n) - 1 */
    uint r = rand();
    r = r % truth_table_length;
    uint j = r;
    printf("Random R is: %u\n", r);
step_4: ;/* and now for the tricky part*/
    for (uint i = 0; i < truth_table_length; i++) {                                 /* 4. For i = 0,...,2^n -1 */
        j = (i + j) % truth_table_length;                                      /*      (a) Let j = (i + r) mod 2^n */
        total_positions++;
        /* (b) If A_j = 0, and 3 <= ord2(A_j)<= maxOrd then: */
        if ((ANF[j] == 0 ) && (3 <= ord2(j, n)) && (ord2(j, n) <= max_ord) ) {      /*      (b) If A_j = 0, and 3 <=ord(A_j)<=maxOrd then: */
            // printf("term index %u passed the test\n", j);
            ANF_positions++;
            memset(ANF_tmp, 0, truth_table_length*sizeof(uint));
            memcpy(ANF_tmp, ANF, truth_table_length*sizeof(uint));
            ANF_tmp[j] = 1;                                                         /*          i. Create A* = A and set A*_j = 1 (the addition of jth term) */
            memset(TT, 0, truth_table_length*sizeof(uint));
            anft_naive(ANF_tmp, TT, n);                                             /*          ii. Find the truth table, f(x), corresponding to A* */
            bool found = true;
            bool TT_weight = true;
            bool f0_weight = true;
            bool f1_weight = true;
            for (int ei = 1; ei < truth_table_length; ei  = ei << 1) {                                        /*          There are n unit vectors with which we can split f(x)        */
                memset(f0, 0, truth_table_length*sizeof(uint)/2);
                memset(f1, 0, truth_table_length*sizeof(uint)/2);
                split_along_direction_e_n(TT, f0, f1, n, ei);
                uint TT_wt = wt(TT, truth_table_length);
                uint f0_wt = wt(f0, truth_table_length/2);
                uint f1_wt = wt(f1, truth_table_length/2);
                TT_weight = (TT_wt == (unsigned)(pow(2,n-1) + pow(2,n/2-1)) ||  (TT_wt == (unsigned)(pow(2,n-1) - pow(2,n/2-1))));
                f0_weight = f0_wt == (unsigned)(pow(2,n-2));
                f1_weight = f1_wt == (unsigned)(pow(2,n-2));
                found = (TT_weight) && (f0_weight || f1_weight);
                if (!found) {
                    break;
                }
            }
            if(found ){
                half_balance_changes++;
                memset(temp_WHT, 0, truth_table_length);
                uint non_lin = nonlinearity(TT, WHT, temp_WHT, truth_table_length, n);          /*              A. Find nonlinearity of f, N_f                  */
                if (non_lin == (unsigned)(pow(2,n-1) - pow(2,n/2-1))) {              /*              B. if N_f = bla bla (i.e. bent) then update A=A*, store f and return to step 2*/
                    memcpy(ANF, ANF_tmp, truth_table_length*sizeof(uint));


                    bent_counter++;
                    uint temp = wt(ANF_tmp, truth_table_length);
                    if (temp > max_bent_awt) {
                        max_bent_awt = temp;
                    }
                    memset(TT, 0, truth_table_length*sizeof(uint));
                    memset(f0, 0, truth_table_length*sizeof(uint)/2);
                    memset(f1, 0, truth_table_length*sizeof(uint)/2);
                    memset(ANF_tmp, 0, truth_table_length*sizeof(uint));
                    goto step_2;
                } else {
                }
            }
        } else {
        }
    }

    free(TT);
    free(f0);
    free(f1);
    free(ANF_tmp);

    printf("Bent counter: %u\n", bent_counter);
    printf("Total positions: %u\n", total_positions);
    printf("ANf positions tested: %u\n", ANF_positions);
    printf("Half balance changes: %u\n", half_balance_changes);
    printf("Maximum awt: %u\n", max_bent_awt);

    return 0;
}


