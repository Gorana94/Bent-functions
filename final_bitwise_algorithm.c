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

typedef struct _split_thread_arg {
    uint *TT;
    uint n;
    long long unsigned length_in_bits;
    uint ei;
    uint split_result;
} split_thread_arg;

typedef struct _anf_thread_arg {
    uint* TT;
    uint* ANF;
    long long unsigned range_start;
    long long unsigned range_end;
    uint n;
} anf_thread_arg;

typedef struct _wht_thread_arg {
    int *WHT;
    uint *TT;
    uint n;
    long long unsigned range_start;
    long long unsigned range_end;
} wht_thread_arg;

typedef struct _inner_wht_thread_arg {
    int *WHT;
    uint stride;
    long long unsigned range_start;
    long long unsigned range_end;
} inner_wht_thread_arg;

/* Global variables */

uint g_bent_nonlinearity = 0;
uint g_higher_correct_TT_wt = 0;
uint g_lower_correct_TT_wt = 0;
uint g_correct_sub_wt = 0;

int A[256];
int W[256][8];

uint count(const char a[], long long unsigned length) 
{ 
    uint n = 0; 
    typedef unsigned long long __attribute__((may_alias)) ull_a; 
    const ull_a *p; 
    // printf("doing a popcount for %ull\n", (long long unsigned)a);
    for (p = (const ull_a *)a; p != (const ull_a *)((long long unsigned)a + length); p++) {
        n += __builtin_popcountll(*p); 
    }
    // printf(" popcount is %u\n", n);
    return n; 
} 

typedef unsigned long long ull;
const ull m1= ~0XAAAAAAAAAAAAAAAA,
          m2= ~0XCCCCCCCCCCCCCCCC,
          m3= ~0XF0F0F0F0F0F0F0F0,
          m4= ~0XFF00FF00FF00FF00,
          m5= ~0XFFFF0000FFFF0000;

void anft_x64 ( ull f[], uint num_of_vars ) {
    int num_of_steps = num_of_vars - 6;
    int num_of_comp_words= 1<<num_of_steps ;
// Step 1 − bitwise ANF on computer words 
    for ( int k= 0; k < num_of_comp_words; k++) {
        ull temp= f[k];
        temp ^= (temp & m1)<<1;
        temp ^= (temp & m2)<<2;
        temp ^= (temp & m3)<<4;
        temp ^= (temp & m4)<<8;
        temp ^= (temp & m5)<<16;
        temp ^= temp<<32;
        f[k]= temp;
    }
 
    // Step 2 
    int blocksize = 1;
    for (int step= 1; step <= num_of_steps ; step++) {
        int source= 0;
        while ( source < num_of_comp_words) {
            int target= source + blocksize ;
            for ( int i= 0; i < blocksize ; i++) {
                f[ target + i ] ^= f[ source + i ] ;
            }
            source += 2 *blocksize;
        }
        blocksize *= 2;
    }
}

void Matrix_W() {

int i, j, p, q, fn=0, w[4][2]={{2,0},{0,-2},{0,2},{-2,0}};
for (i=0; i<4; i++)
    for (j=0; j<4; j++)
        for (p=0; p<4; p++)
            for (q=0; q<4; q++) {
                W[fn][0]=w[q][0]+w[p][0]+w[j][0]+w[i][0];
                W[fn][1]=w[q][1]+w[p][1]+w[j][1]+w[i][1];
                W[fn][2]=w[q][0]-w[p][0]+w[j][0]-w[i][0];
                W[fn][3]=w[q][1]-w[p][1]+w[j][1]-w[i][1];
                W[fn][4]=w[q][0]+w[p][0]-w[j][0]-w[i][0];
                W[fn][5]=w[q][1]+w[p][1]-w[j][1]-w[i][1];
                W[fn][6]=w[q][0]-w[p][0]-w[j][0]+w[i][0];
                W[fn][7]=w[q][1]-w[p][1]-w[j][1]+w[i][1];
                fn++;
            }

}; 

void inner_wht_loop (uint range_start, uint range_end, uint stride, int* wht) {
    for (int k = range_start; k < range_end; k++) {
        int a = wht[k] - wht[k+stride];
        wht[k] = wht[k] + wht[k+stride];
        wht[k+stride] = a;
    }
}

void *WHT_thread(void *arg)
{
    wht_thread_arg* t_arg = (wht_thread_arg*)arg;

    int* wht = t_arg->WHT;
    uint* tt = t_arg->TT;

    long long int i, j, k, a, b, l, N=((long long unsigned)1<<t_arg->n)>>5;
    for (i=t_arg->range_start; i<t_arg->range_end; i++) {
        l=i*32;
        a=tt[i]&0xFF;
        b=(tt[i]&0xFF00)>>8;
        j=(tt[i]&0xFF0000)>>16;
        k=(tt[i]&0xFF000000)>>24;
        wht[l]=W[a][0]+W[b][0]+W[j][0]+W[k][0];
        wht[l+1]=W[a][1]+W[b][1]+W[j][1]+W[k][1];
        wht[l+2]=W[a][2]+W[b][2]+W[j][2]+W[k][2];
        wht[l+3]=W[a][3]+W[b][3]+W[j][3]+W[k][3];
        wht[l+4]=W[a][4]+W[b][4]+W[j][4]+W[k][4];
        wht[l+5]=W[a][5]+W[b][5]+W[j][5]+W[k][5];
        wht[l+6]=W[a][6]+W[b][6]+W[j][6]+W[k][6];
        wht[l+7]=W[a][7]+W[b][7]+W[j][7]+W[k][7];
        wht[l+8]=W[a][0]-W[b][0]+W[j][0]-W[k][0];
        wht[l+9]=W[a][1]-W[b][1]+W[j][1]-W[k][1];
        wht[l+10]=W[a][2]-W[b][2]+W[j][2]-W[k][2];
        wht[l+11]=W[a][3]-W[b][3]+W[j][3]-W[k][3];
        wht[l+12]=W[a][4]-W[b][4]+W[j][4]-W[k][4];
        wht[l+13]=W[a][5]-W[b][5]+W[j][5]-W[k][5];
        wht[l+14]=W[a][6]-W[b][6]+W[j][6]-W[k][6];
        wht[l+15]=W[a][7]-W[b][7]+W[j][7]-W[k][7];
        wht[l+16]=W[a][0]+W[b][0]-W[j][0]-W[k][0];
        wht[l+17]=W[a][1]+W[b][1]-W[j][1]-W[k][1];
        wht[l+18]=W[a][2]+W[b][2]-W[j][2]-W[k][2];
        wht[l+19]=W[a][3]+W[b][3]-W[j][3]-W[k][3];
        wht[l+20]=W[a][4]+W[b][4]-W[j][4]-W[k][4];
        wht[l+21]=W[a][5]+W[b][5]-W[j][5]-W[k][5];
        wht[l+22]=W[a][6]+W[b][6]-W[j][6]-W[k][6];
        wht[l+23]=W[a][7]+W[b][7]-W[j][7]-W[k][7];
        wht[l+24]=W[a][0]-W[b][0]-W[j][0]+W[k][0];
        wht[l+25]=W[a][1]-W[b][1]-W[j][1]+W[k][1];
        wht[l+26]=W[a][2]-W[b][2]-W[j][2]+W[k][2];
        wht[l+27]=W[a][3]-W[b][3]-W[j][3]+W[k][3];
        wht[l+28]=W[a][4]-W[b][4]-W[j][4]+W[k][4];
        wht[l+29]=W[a][5]-W[b][5]-W[j][5]+W[k][5];
        wht[l+30]=W[a][6]-W[b][6]-W[j][6]+W[k][6];
        wht[l+31]=W[a][7]-W[b][7]-W[j][7]+W[k][7];
    }

}

void* inner_wht_thread(void* arg)
{
    inner_wht_thread_arg* t_arg = (inner_wht_thread_arg*)arg;
    int *wht = t_arg->WHT;
    int stride = t_arg->stride;
    for (int k = t_arg->range_start; k < t_arg->range_end; k++){
        int a = wht[k] - wht[k+stride];
        wht[k] = wht[k] + wht[k+stride];
        wht[k+stride] = a;
    }
}

void fully_threaded_WHT(int *tt, int *wht, int n) {

    uint num_cores = getNumCores();
    if (num_cores > 16) {
        num_cores = 16;
    }
    long long unsigned N=1<<(n-5);
    uint range = 0;
    uint modulo_range = 0;
    if (N < num_cores) {
        range = N;
        num_cores = 1;
        modulo_range = 0;
    } else {
        range = N / num_cores;
        modulo_range = N % range;
    }

    /* Create thread objects */
    pthread_t* threads = malloc(num_cores * sizeof(pthread_t));
    wht_thread_arg* wht_thread_args = malloc(num_cores * sizeof(wht_thread_arg));
    void **proforme_args;
    proforme_args = (void**) malloc(num_cores * sizeof(void*));

    pthread_t modulo_thread;
    wht_thread_arg modulo_arg;
    void *modulo_proforme_arg;
    for (int i = 0; i < num_cores; i++) {
        wht_thread_args[i].WHT = wht;
        wht_thread_args[i].TT = tt;
        wht_thread_args[i].n = n;
        wht_thread_args[i].range_start = i*range;
        wht_thread_args[i].range_end = (i+1)*range;
        pthread_create((threads + i), NULL, &WHT_thread, (wht_thread_args+i));
    }
    if (modulo_range > 0) {
        modulo_arg.WHT = wht;
        modulo_arg.TT = tt;
        modulo_arg.n = n;
        modulo_arg.range_start = num_cores * range;
        modulo_arg.range_end = N;
        pthread_create(&modulo_thread, NULL, &WHT_thread, &modulo_arg);
    }
    for(int i = 0; i < num_cores; i++){
        pthread_join(threads[i], proforme_args+i);
    }
    if (modulo_range > 0) {
        pthread_join(modulo_thread, &modulo_proforme_arg);
    }
    free(threads);
    free(wht_thread_args);
    free(proforme_args);

    pthread_t* inner_pthreads = malloc(num_cores * sizeof(pthread_t));
    inner_wht_thread_arg* inner_args = malloc(num_cores * sizeof(inner_wht_thread_arg));
    void **inner_proforme_args = malloc(num_cores * sizeof(void*));

    long long int i, j, k, a, b, l, size=((long long unsigned)1<<n);

    for (i=32; i<size; i=i<<1) {
        /*we are gonna paralelize the inner two loops of WHT, they'll be called outer and inner loop from now*/
        /* whichever loop has more steps shall be parallelized by available cores (16) */
        uint total_outer_steps = size / (i<<1);
        uint total_inner_steps = i;
        uint num_outer_threads = num_cores; /* To be either 1 or num_cores */
        uint inner_threads = 1; /* To be either 1 or num_cores */
        uint steps_per_outer_thread = total_outer_steps / num_outer_threads;

        for (j = 0; j < size; j+=(i<<1)){
            if (i > 1024*1024*2) { 
                inner_threads = 16;
                uint steps_per_inner_thread = total_inner_steps / inner_threads;
                for (int t= 0;  t < inner_threads; t++) {
                    inner_args[t].stride = i;
                    inner_args[t].WHT = wht;
                    inner_args[t].range_start = j + t * steps_per_inner_thread;
                    inner_args[t].range_end = j + (t+1)*steps_per_inner_thread;
                    pthread_create((inner_pthreads+t), NULL, &inner_wht_thread, (inner_args+t));
                }
                for (int t = 0; t < inner_threads; t++) {
                    pthread_join(inner_pthreads[t], (inner_proforme_args+t));
                }
            } else { 
                for (k=j; k < j + i; k++) {
                    a = wht[k] - wht[k+i];
                    wht[k] = wht[k] + wht[k+i];
                    wht[k+i] = a;
                }
            }
        }
    }

    free(inner_pthreads);
    free(inner_args);
    free(inner_proforme_args);

};


void WHT(int *tt, int *wht, int n) {

    long long int i, j, k, a, b, l, size=((long long unsigned)1<<n), N=size>>5;
    for (i=0; i<N; i++) {
        l=i*32;
        a=tt[i]&0xFF;
        b=(tt[i]&0xFF00)>>8;
        j=(tt[i]&0xFF0000)>>16;
        k=(tt[i]&0xFF000000)>>24;
        wht[l]=W[a][0]+W[b][0]+W[j][0]+W[k][0];
        wht[l+1]=W[a][1]+W[b][1]+W[j][1]+W[k][1];
        wht[l+2]=W[a][2]+W[b][2]+W[j][2]+W[k][2];
        wht[l+3]=W[a][3]+W[b][3]+W[j][3]+W[k][3];
        wht[l+4]=W[a][4]+W[b][4]+W[j][4]+W[k][4];
        wht[l+5]=W[a][5]+W[b][5]+W[j][5]+W[k][5];
        wht[l+6]=W[a][6]+W[b][6]+W[j][6]+W[k][6];
        wht[l+7]=W[a][7]+W[b][7]+W[j][7]+W[k][7];
        wht[l+8]=W[a][0]-W[b][0]+W[j][0]-W[k][0];
        wht[l+9]=W[a][1]-W[b][1]+W[j][1]-W[k][1];
        wht[l+10]=W[a][2]-W[b][2]+W[j][2]-W[k][2];
        wht[l+11]=W[a][3]-W[b][3]+W[j][3]-W[k][3];
        wht[l+12]=W[a][4]-W[b][4]+W[j][4]-W[k][4];
        wht[l+13]=W[a][5]-W[b][5]+W[j][5]-W[k][5];
        wht[l+14]=W[a][6]-W[b][6]+W[j][6]-W[k][6];
        wht[l+15]=W[a][7]-W[b][7]+W[j][7]-W[k][7];
        wht[l+16]=W[a][0]+W[b][0]-W[j][0]-W[k][0];
        wht[l+17]=W[a][1]+W[b][1]-W[j][1]-W[k][1];
        wht[l+18]=W[a][2]+W[b][2]-W[j][2]-W[k][2];
        wht[l+19]=W[a][3]+W[b][3]-W[j][3]-W[k][3];
        wht[l+20]=W[a][4]+W[b][4]-W[j][4]-W[k][4];
        wht[l+21]=W[a][5]+W[b][5]-W[j][5]-W[k][5];
        wht[l+22]=W[a][6]+W[b][6]-W[j][6]-W[k][6];
        wht[l+23]=W[a][7]+W[b][7]-W[j][7]-W[k][7];
        wht[l+24]=W[a][0]-W[b][0]-W[j][0]+W[k][0];
        wht[l+25]=W[a][1]-W[b][1]-W[j][1]+W[k][1];
        wht[l+26]=W[a][2]-W[b][2]-W[j][2]+W[k][2];
        wht[l+27]=W[a][3]-W[b][3]-W[j][3]+W[k][3];
        wht[l+28]=W[a][4]-W[b][4]-W[j][4]+W[k][4];
        wht[l+29]=W[a][5]-W[b][5]-W[j][5]+W[k][5];
        wht[l+30]=W[a][6]-W[b][6]-W[j][6]+W[k][6];
        wht[l+31]=W[a][7]-W[b][7]-W[j][7]+W[k][7];
    }

    for (i=32; i<size; i=i<<1)
        for (j=0; j<size; j=j+(i<<1))
            for (k=j; k<j+i; k++) {
                a=wht[k]-wht[k+i];
                wht[k]=wht[k]+wht[k+i];
                wht[k+i]=a;
            }
};

void Matrix_A() 
{
    int i, j, fn = 0, a[16] = {0,15,10,5,12,3,6,9,8,7,2,13,4,11,14,1};
    for (i = 0; i < 16; i++) {
        for (j = 0; j < 16; j++){
            A[fn++] = a[j]^((a[i]^a[j]) << 4);
        }
    }
}

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

static inline uint scalar_product_unit_vector(long long unsigned x, uint unit_vector) {
    return x & unit_vector;
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

void populate_A() {
    int fn=0, a[16] = {0,15,10,5,12,3,6,9,8,7,2,13,4,11,14,1};
    for (int i = 0; i < 16; i++) {
        for (int j = 0; j < 16; j++) {
            A[fn++]=a[j]^((a[i]^a[j]) << 4);
        }
    }
}

void ANFT(uint *tt, uint *anf, int n) {
    long long unsigned i, j, k, l, m, N=1<<(n-5);
    for (i=0; i<N; i++) {
        j=tt[i]&0x000000FF;
        k=(tt[i]&0x0000FF00)>>8;
        l=(tt[i]&0x00FF0000)>>16;
        m=(tt[i]&0xFF000000)>>24;
        anf[i]=A[j] | ((A[j] ^ A[k]) << 8) | ((A[j] ^ A[l]) << 16)
            | ((A[j] ^ A[k] ^ A[l] ^ A[m]) << 24);
    }
    for (i=1; i<N; i=i<<1){
        for (j=0; j<N; j=j+(i<<1)){
            for (k=j; k<j+i; k++) {
                anf[k+i]^=anf[k];
            }
        }
    }
}

void *ANF_thread(void* arg) {
    anf_thread_arg* t_arg = (anf_thread_arg*)arg;
    uint *tt = t_arg->TT;
    uint *anf = t_arg->ANF;
    long long unsigned i, j, k, l, m;
    for (i=t_arg->range_start; i<t_arg->range_end; i++) {
        j=tt[i]&0x000000FF;
        k=(tt[i]&0x0000FF00)>>8;
        l=(tt[i]&0x00FF0000)>>16;
        m=(tt[i]&0xFF000000)>>24;
        anf[i]=A[j] | ((A[j] ^ A[k]) << 8) | ((A[j] ^ A[l]) << 16)
            | ((A[j] ^ A[k] ^ A[l] ^ A[m]) << 24);
    }
}

void threaded_ANFT(uint *tt, uint *anf, int n) {

    uint num_cores = getNumCores();
    if (num_cores > 16) {
        num_cores = 16;
    }
    long long unsigned N=1<<(n-5);
    uint range = 0;
    uint modulo_range = 0;
    if (N < num_cores) {
        range = N;
        num_cores = 1;
        modulo_range = 0;
    } else {
        range = N / num_cores;
        modulo_range = N % range;
    }

    /* Create thread objects */
    pthread_t* threads = malloc(num_cores * sizeof(pthread_t));
    anf_thread_arg* anf_thread_args = malloc(num_cores * sizeof(anf_thread_arg));
    void **proforme_args;
    proforme_args = (void**) malloc(num_cores * sizeof(void*));

    pthread_t modulo_thread;
    anf_thread_arg modulo_arg;
    void *modulo_proforme_arg;
    for (int i = 0; i < num_cores; i++) {
        anf_thread_args[i].TT = tt;
        anf_thread_args[i].ANF = anf;
        anf_thread_args[i].range_start = i*range;
        anf_thread_args[i].range_end = (i+1)*range;
        pthread_create((threads + i), NULL, &ANF_thread, (anf_thread_args+i));
    }
    if (modulo_range > 0) {
        modulo_arg.TT = tt;
        modulo_arg.ANF = anf;
        modulo_arg.range_start = num_cores * range;
        modulo_arg.range_end = N;
        pthread_create(&modulo_thread, NULL, &ANF_thread, &modulo_arg);
    }
    for(int i = 0; i < num_cores; i++){
        pthread_join(threads[i], proforme_args+i);
    }
    if (modulo_range > 0) {
        pthread_join(modulo_thread, &modulo_proforme_arg);
    }
    free(threads);
    free(anf_thread_args);
    free(proforme_args);
    long long unsigned i, j, k, l, m;
    
    for (i=1; i<N; i=i<<1){
        for (j=0; j<N; j=j+(i<<1)){
            for (k=j; k<j+i; k++) {
                anf[k+i]^=anf[k];
            }
        }
    }
}

void print_packed_array(uint* packed_array, uint packed_array_length)
{
    int total = 0;

    for (long long unsigned i = packed_array_length-1; i >= 0; i--) {
        for (int j = (sizeof(uint)*8)-1; j >= 0; j--) {
            printf("%u", (packed_array[i] >> j) & 1);
            total++;
        }
    }
            return;
}

static inline void set_bit(uint* packed_array, long long unsigned bit_index, uint value)
{
    if (value) {
        packed_array[bit_index >> 5] |= (1<<(bit_index&31));
    } else {
        packed_array[bit_index >> 5] &= ~(1<<(bit_index&31));
    }
}

static inline uint get_bit(uint* packed_array, long long unsigned bit_index)
{
    return !!((1<<(bit_index&31)) & (packed_array[bit_index >> 5]));
}

void bit_ANF(uint* tt, uint* anf, uint n)
{
    int i, j, k, l, m, N = 1 << (n-5);
    for (i = 0; i < N; i++) {
        j = tt[i] & 0x000000FF;
        k = (tt[i] & 0x0000FF00) >> 8;
        l = (tt[i] & 0x00FF0000) >> 16;
        m = (tt[i] & 0xFF000000) >> 24;
        anf[i] = A[j] | ((A[j] ^ A[k]) << 8) | ((A[j] ^ A[l]) << 16) 
                | ((A[j] ^ A[k] ^ A[l] ^ A[m]) << 24);
    }

    for (i = 1; i < N; i = i << 1) {
        for (j = 0; j < N; j = j + (i << 1)) {
            for (k = j; k < j+i; k++) {
                anf[k+i] ^= anf[k];
            }
        }
    }
}

uint bit_nonlinearity(uint* tt, int* wht, uint n, long long unsigned length_in_bits) {
    if (n < 20) {
        WHT(tt, wht, n);
    } else {
        fully_threaded_WHT(tt, wht, n);
    }
    uint max_wht = 0;
    long long unsigned the_only_correct_wht_value_for_an_element_of_wht_of_a_bent_function = (long long unsigned)pow(2, n/2);
    for (long long unsigned i = 0; i < length_in_bits; i++) {
        if (abs(wht[i]) > max_wht) {
            max_wht = abs(wht[i]);
            if(abs(wht[i]) != the_only_correct_wht_value_for_an_element_of_wht_of_a_bent_function) {
                return 0; /*false - not a bent function*/
            }
        }
    }

    return 1; /* true - its a bent function */
}

uint popcount_split(uint* TT, long long unsigned ei, uint n, long long unsigned length_in_bits)
{
    uint subfunction_condition = 1; /*true*/
    uint subs[2] = {0,0};
    uint num_bytes = length_in_bits / 8 ;
    uint num_chunks = length_in_bits / ei;
    uint length_of_chunk = num_bytes / num_chunks;

    for (int i = 0; i < num_chunks; i++) {
        subs[i%2] += count((char *)TT + i * length_of_chunk, length_of_chunk);
    }

    subfunction_condition = (subs[0] == g_correct_sub_wt) || (subs[1] == g_correct_sub_wt);

    return subfunction_condition;
}


void* popcount_split2(void *arg)
{
    split_thread_arg* split_args = (split_thread_arg*) arg;
    uint subs[2] = {0,0};
    uint num_bytes = split_args->length_in_bits / 8 ;
    uint num_chunks = split_args->length_in_bits / split_args->ei;
    uint length_of_chunk = num_bytes / num_chunks;

    for (int i = 0; i < num_chunks; i++) {
        subs[i%2] += count((char *)(split_args->TT) + i * length_of_chunk, length_of_chunk);
    }

    split_args->split_result = (subs[0] == g_correct_sub_wt) || (subs[1] == g_correct_sub_wt);
    return NULL;
}

uint smart_split_and_count(uint* TT, long long unsigned ei, uint n, long long unsigned length_in_bits)
{   
    uint sub_weights[2] = {0, 0}; // sub_weights[0] is f0 weight, sub_weights[1] is f1 weight
    uint bit = 0;
    uint scalar_product = 0;
    uint subfunction_condition = 0;
    for (long long unsigned x = 0; x < length_in_bits; x++) {
        bit = get_bit(TT, x);
        scalar_product = scalar_product_unit_vector(x, ei);
        sub_weights[!!scalar_product] += bit;
    }

    subfunction_condition = (sub_weights[0] == g_correct_sub_wt) || (sub_weights[1] == g_correct_sub_wt);
    return subfunction_condition;
}

void* iterative_split(void* arg)
{   
    split_thread_arg* split_args = (split_thread_arg*) arg;
    uint sub_weights[2] = {0, 0}; // sub_weights[0] is f0 weight, sub_weights[1] is f1 weight
    uint bit = 0;
    uint scalar_product = 0;
    for (long long unsigned x = 0; x < split_args->length_in_bits; x++) {
        bit = get_bit(split_args->TT, x);
        scalar_product = scalar_product_unit_vector(x, split_args->ei);
        sub_weights[!!scalar_product] += bit;
    }

    split_args->split_result = (sub_weights[0] == g_correct_sub_wt) || (sub_weights[1] == g_correct_sub_wt);

    return NULL;
}

uint combined_split(uint* TT, uint n, long long unsigned length_in_bits)
{
    uint TT_wt = count((char *)TT, length_in_bits/8);
    uint TT_wt_condition = (TT_wt == g_higher_correct_TT_wt) || (TT_wt == g_lower_correct_TT_wt);
    if (!TT_wt_condition) {
        return 0; /*false*/
    }

    uint result = 1; /*true*/
    for (long long unsigned ei = 1; ei < length_in_bits; ei = ei << 1){                                             /*          iii. if f=[f0|f1] has wt(f) = bla bla and wt(f0) or wt(f1) == bla bla FOR ALL ei*/
        if (ei < sizeof(long long unsigned) * 8) {
            result = smart_split_and_count(TT, ei, n, length_in_bits);
        } else {
            result = popcount_split(TT, ei, n, length_in_bits);
        }
        if (!result) {
            return result;
        }
    }
    return result;
}

uint combined_split2(uint* TT, uint n, long long unsigned length_in_bits)
{
    uint TT_wt = count((char *)TT, length_in_bits/8);
    uint TT_wt_condition = (TT_wt == g_higher_correct_TT_wt) || (TT_wt == g_lower_correct_TT_wt);
    if (!TT_wt_condition) {
        return 0; /*false*/
    }

    /* A batch of threads is defined by number of logical cores, and
     * total number of threads needed is equal to n, so we divide n by number of cores
     * to get the number of batches where all cores will be used at once
     * Some threads might remain otside of a batch and they are executed afterwards - modulo batch*/

    uint complete_batch_size = getNumCores();
    if (complete_batch_size > 16) {
        complete_batch_size = 16;
    }
    uint num_batches = n / complete_batch_size;
    uint modulo_batch_size = n % complete_batch_size;
    void **proforme_arg;

    /* unit vector ei determines how long are chunks of the truth table thatwill
     * belong to subfunciont f0 or f1. If chunk size is smaller than 64 bits we cannot
     * use the popcnt function*/
    uint smallest_chunk_size = sizeof(long long unsigned) * 8;
    pthread_t* threads = (pthread_t*)malloc(complete_batch_size * sizeof(pthread_t));
    split_thread_arg* split_args = (split_thread_arg*)malloc(complete_batch_size * sizeof(split_thread_arg));
    proforme_arg = (void**) malloc(complete_batch_size * sizeof(void*));
    long long unsigned ei = 1;

    for (int i = 0; i < num_batches; i++) {
        for (int j = 0; j < complete_batch_size; j++){
            split_args[j].ei = ei;
            split_args[j].length_in_bits = length_in_bits;
            split_args[j].n = n;
            split_args[j].TT = TT;
            if (ei < smallest_chunk_size) {
                pthread_create((threads + j), NULL, &iterative_split, (split_args + j));
            } else {
                pthread_create((threads + j), NULL, &popcount_split2, (split_args + j));
            }
            ei = ei << 1;
        }

        for (int j = 0; j < complete_batch_size; j++) {
            pthread_join(threads[j], (proforme_arg + j));
        }
        for (int j = 0; j < complete_batch_size; j++) {
            if (!split_args[j].split_result) {
                free(split_args);
                free(threads);
                free(proforme_arg);
                return 0; /*A bad split*/
            }
        }
    }
    /*Now handle rest of threads (if any) that didnt fill out a complete batch*/
    for (int i = 0; i < modulo_batch_size; i++) {
        split_args[i].ei = ei;
        split_args[i].length_in_bits = length_in_bits;
        split_args[i].n = n;
        split_args[i].TT = TT;
        if (ei < smallest_chunk_size) {
            pthread_create((threads + i), NULL, &iterative_split, (split_args + i));
        } else {
            pthread_create((threads + i), NULL, &popcount_split2, (split_args + i));
        }
        ei = ei << 1;
    }
    for (int i = 0; i < modulo_batch_size; i++) {
        pthread_join(threads[i], (proforme_arg + i));
    }

    for (int i = 0; i < modulo_batch_size; i++) {
        if (!split_args[i].split_result) {
            free(split_args);
            free(threads);
            free(proforme_arg);
            return 0; /* A bad split */
        }
    }

    free(split_args);
    free(threads);
    free(proforme_arg);

    return 1; /*A good split*/
}

int main(int argc, char** argv)
{
    uint n = 8;
    uint maxOrd = 3;
    if (argc == 2) {
        n = atoi(argv[1]);
    } else if (argc == 3) {
        n = atoi(argv[1]);
        maxOrd = atoi(argv[2]);
    }

    Matrix_A();
    Matrix_W();
    printf("x64 ANFT, tredovan split, Biramo obican WHT ili WHT tredovan i prvi i drugi deo u zavisnosti od n\n");
    const long long unsigned length_in_bits = (long long unsigned)pow(2,n);
    const uint length_in_bytes = length_in_bits/8;
    const uint length_in_uints = length_in_bytes / sizeof(uint);

    uint *ANF;
    uint *TT;
    int *WHT;
    ANF = (uint *)malloc(length_in_bytes);
    TT = (uint *)malloc(length_in_bytes);
    WHT = (int *)malloc(length_in_bits*sizeof(int));

    g_bent_nonlinearity = (uint)pow(2, n-1) - (uint)pow(2, n/2 - 1);
    g_higher_correct_TT_wt = (uint)pow(2,n-1) + (uint)pow(2,n/2 - 1);
    g_lower_correct_TT_wt = (uint)pow(2,n-1) - (uint)pow(2,n/2 - 1);
    g_correct_sub_wt = (uint)pow(2, n-2);

    srand(time(0));

    struct timespec start_time, end_time, loop_start, loop_end;

    long long unsigned new_term_index = 0;
    uint new_term = 0;
    uint new_term_order = 0;

    long long unsigned ANF_positions_tested = 0;
    long long unsigned half_balance_yielding_changes = 0;
    long long unsigned bent_counter = 0;
    long long unsigned total_iterations = 0;


    memset(ANF, 0, length_in_bytes);
    memset(TT, 0, length_in_bytes);
    memset(WHT, 0, length_in_bits*sizeof(int));

    clock_gettime(CLOCK_REALTIME, &start_time);
    step_1: ;

    /* 1. Create ANF array of size 2^n and initialize: A_i = 0, 0 <= i <= 2^n - 1*/
    memset(ANF, 0, length_in_bytes);

    /* 2. For i = 0,...,N/2-1      ORIGINALLY 2^(N/2-1) */
    step_2: ;
    for (int i = 0; i <= (n/2 - 1); i++) {
        /*      (a) Let j = 2^i + 2^(n/2+i) */
        long long unsigned j = (uint)pow(2,i) + (uint)pow(2, n/2+i);
        /*      (b) Set A_j = 1 */
        set_bit(ANF, j, 1);
    }
    /* 3. Select a random r, 0<= r <= (2^n) - 1 */
    step_3: ;
    long long unsigned r = rand();
    new_term_index = r % length_in_bits;

    /* 4. For i = 0,...,2^n -1 */
    step_4: ;
    for (long long unsigned i = 0; i <= length_in_bits - 1; i++) {
        /*      (a) Let j = (i + r) mod 2^n */
        new_term_index = (i + new_term_index) % length_in_bits;
        new_term = get_bit(ANF, new_term_index);
        new_term_order = ord2(new_term_index, n);
        total_iterations++;
        /*      (b) If A_j = 0, and 3 <=ord(A_j)<=maxOrd then: */
        if ((new_term == 0) && (new_term_order >= 3) && (maxOrd >= new_term_order)) {
            ANF_positions_tested++;
            /*          i. Create A* = A and set A*_j = 1 (the addition of jth term) */
            set_bit(ANF, new_term_index, 1);
            memset(TT, 0, length_in_bytes);
            memcpy(TT, ANF, length_in_bytes);
            /*          ii. Find the truth table, f(x), corresponding to A* */
            anft_x64((ull*)TT, n);
            int good_split = 1;
            /* iii. if f=[f0|f1] has 
             * wt(f) == 2^(n-1) + 2^(n/2-1) or wt(f) == 2^(n-1) - 2^(n/2-1)
             * and
             * wt(f0) or wt(f1) == 2^(n-2) */
            good_split = combined_split2(TT, n , length_in_bits);
            if (good_split) {
                half_balance_yielding_changes++;
                /*              A. Find nonlinearity of f, N_f*/
                uint non_lin = bit_nonlinearity(TT, WHT, n, length_in_bits);
                /*              B. if N_f = 2^(n-1) - 2^(n/2 -1) then update A=A*, store f and return to step 2*/
                if (non_lin) {
                    bent_counter++;
                    goto step_2;
                } 
            }
            /*              C. If N_F is not good nonlinearity, then revert change to A */
            set_bit(ANF, new_term_index, 0);
        } else {
        }
    }
    clock_gettime(CLOCK_REALTIME, &end_time);

    int seconds = end_time.tv_sec - start_time.tv_sec;
    long nanoseconds = 0;
    if ((end_time.tv_nsec - start_time.tv_nsec) < 0) {
        nanoseconds = 1000000000 + end_time.tv_nsec - start_time.tv_nsec;
        seconds--;
    } else {
        nanoseconds = end_time.tv_nsec - start_time.tv_nsec; 
    }

    printf("total iterations: %u ANF positions tested: %u, Half Balance changes: %u, Bent functions: %u\n",total_iterations, ANF_positions_tested, half_balance_yielding_changes, bent_counter);

    printf("Total time was %d seconds, %ldms, %ldus, %ldns\n", seconds, nanoseconds/1000000, (nanoseconds/1000)%1000, nanoseconds%1000);

free(TT);
free(ANF);
free(WHT);

return 0;
}
