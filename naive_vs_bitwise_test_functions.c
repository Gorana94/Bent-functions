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
    long long unsigned ei;
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
    long long unsigned stride;
    long long unsigned range_start;
    long long unsigned range_end;
} inner_wht_thread_arg;

typedef struct _outer_wht_thread_arg {
    int *WHT;
    long long unsigned range_start;
    long long unsigned range_end;
    long long unsigned stride;
    uint inner_threads;
} outer_wht_thread_arg;

/* Global variables */

long long unsigned g_bent_nonlinearity = 0;
long long unsigned g_higher_correct_TT_wt = 0;
long long unsigned g_lower_correct_TT_wt = 0;
long long unsigned g_correct_sub_wt = 0;

long long unsigned A[256];
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
        //   m6= 0XFFFFFFFF00000000
// const int num_of_vars= 28;
 //number of variables 
// const int num_of_steps= num_of_vars - 6;
 //for Step (2) of ANFT 
// const int num_of_comp_words= 1<<num_of_steps ;
 //=4 for 8 vars 
// ull f [num_of_comp_words] ;
 // for representation of the TT( f )
void OPTIMIZED_ANF ( ull f[], uint num_of_vars ) {
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

void WHT(int *tt, int *wht, int n) {

    long long int i, j, k, a, b, l, size=((long long unsigned)1<<n), N=size>>5;
    // printf("u wht pre prve petlje, size je %llu a N je %llu\n", size, N);
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
// printf("u wht posle prve petlje\n");
// printf("u wht pre druge petlje\n");
long long unsigned total_wht_iterator = 0;
    for (i=32; i<size; i=i<<1)
        for (j=0; j<size; j=j+(i<<1))
            for (k=j; k<j+i; k++) {
                total_wht_iterator++;
                a=wht[k]-wht[k+i];
                wht[k]=wht[k]+wht[k+i];
                wht[k+i]=a;
            }
// printf("u wht posle druge petlje iukupno se u noj izvrsi %llu koraka\n", total_wht_iterator);
};

void *WHT_thread(void *arg)
{
    wht_thread_arg* t_arg = (wht_thread_arg*)arg;

    int* wht = t_arg->WHT;
    uint* tt = t_arg->TT;

    long long int i, j, k, a, b, l, N=((long long unsigned)1<<t_arg->n)>>5;
    // printf("u wht pre prve petlje, size je %llu a N je %llu\n", size, N);
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
    long long unsigned stride = t_arg->stride;
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

    // pthread_t* outer_pthreads = malloc(num_cores * sizeof(pthread_t));
    // outer_wht_thread_arg* outer_args = malloc(num_cores * sizeof(outer_wht_thread_arg));
    // void **outer_proforme_args = malloc(num_cores * sizeof(void*));

    pthread_t* inner_pthreads = malloc(num_cores * sizeof(pthread_t));
    inner_wht_thread_arg* inner_args = malloc(num_cores * sizeof(inner_wht_thread_arg));
    void **inner_proforme_args = malloc(num_cores * sizeof(void*));

    long long int i, j, k, a, b, l, size=((long long unsigned)1<<n);

    for (i=32; i<size; i=i<<1) {
        /*we are gonna paralelize the inner two loops of WHT, they'll be called outer and inner loop from now*/
        /* whichever loop has more steps shall be parallelized by available cores (16) */
        uint total_outer_steps = size / (i<<1);
        // printf("total outer steps %u\n", total_outer_steps);
        uint total_inner_steps = i;
        // printf("total inner_steps %u\n", total_inner_steps);
        uint num_outer_threads = num_cores; /* To be either 1 or num_cores */
        uint inner_threads = 1; /* To be either 1 or num_cores */
        uint steps_per_outer_thread = total_outer_steps / num_outer_threads;

        for (j = 0; j < size; j+=(i<<1)){
            if (i > 1024*1024*2) { /*paralelno unutrasnji*/
                inner_threads = 16;
                uint steps_per_inner_thread = total_inner_steps / inner_threads;
                for (int t= 0;  t < inner_threads; t++) {
                    inner_args[t].stride = i;
                    inner_args[t].WHT = wht;
                    inner_args[t].range_start = j + t * steps_per_inner_thread;
                    inner_args[t].range_end = j + (t+1)*steps_per_inner_thread;
                    pthread_create((inner_pthreads+t), NULL, &inner_wht_thread, (inner_args+t));
                    // inner_wht_loop(inner_args[t].range_start, inner_args[t].range_end, i, wht);
                }
                for (int t = 0; t < inner_threads; t++) {
                    pthread_join(inner_pthreads[t], (inner_proforme_args+t));
                }
            } else { /*radi sve normalno*/
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

// printf("u wht posle druge petlje iukupno se u noj izvrsi %llu koraka\n", total_wht_iterator);
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
        // if(i == 0) {
        //     printf("ovo su Aovi %llu %llu %llu %llu\n", A[j], A[k], A[l], A[m]);
        //     printf("a indeksi j, k, l, m su %llu %llu %llu %llu\n", j, k, l, m);
        // }
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

void ANFT_x64(long long unsigned *tt, long long unsigned *anf, int n) {
    long long unsigned i, j, k, l, m, o, p, q, r, N=1<<(n-6);
    for (i=0; i<N; i++) {
        j=tt[i]&0x000000FF;
        k=(tt[i]&0x0000FF00)>>8;
        l=(tt[i]&0x00FF0000)>>16;
        m=(tt[i]&0xFF000000)>>24;
        o=(tt[i]&0xFF00000000)>>32;
        p=(tt[i]&0xFF0000000000)>>40;
        q=(tt[i]&0xFF000000000000)>>48;
        r=(tt[i]&0xFF00000000000000)>>56;

        anf[i] = A[j]                               |
                ((A[j] ^ A[k]) << 8)                |
                ((A[j] ^ A[l]) << 16)               |
                ((A[j] ^ A[k] ^ A[m] ^ A[l]) << 24) |
                ((A[j] ^ A[o]) << 32)               |
                ((A[j] ^ A[k] ^ A[o] ^ A[p]) << 40) |
                ((A[j] ^ A[l] ^ A[o] ^ A[q]) << 48) |
                ((A[j] ^ A[k] ^ A[l] ^ A[m] ^ A[o] ^ A[p] ^ A[q] ^ A[r]) << 56);

        // anf[i]=A[j] | ((A[j] ^ A[k]) << 8) | ((A[j] ^ A[l]) << 16)
        //     | ((A[j] ^ A[k] ^ A[l] ^ A[m]) << 24);
    }
    for (i=1; i<N; i=i<<1){
        for (j=0; j<N; j=j+(i<<1)){
            for (k=j; k<j+i; k++) {
                anf[k+i]^=anf[k];
            }
        }
    }
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


static inline void joanne_set_bit(uint* packed_array, long long unsigned bit_index, uint value)
{
    if (value) {
        packed_array[bit_index >> 5] |= (1<<(bit_index&31));
    } else {
        packed_array[bit_index >> 5] &= ~(1<<(bit_index&31));
    }
}

static inline uint joanne_get_bit(uint* packed_array, long long unsigned bit_index)
{
    return !!((1<<(bit_index&31)) & (packed_array[bit_index >> 5]));
}


uint scalar_product_unit_vector(long long unsigned x, uint unit_vector) {
    return x & unit_vector;
}

uint flipped_scalar_product_unit_vector(long long unsigned x, uint unit_vector, uint n)
{
    return (x >> (n - unit_vector - 1)) & 1;
}

void split_along_direction_e_n(uint * f, uint n, uint* f_wt, uint ei)  /* TODO Split_along_direction */
{
    uint f0_counter = 0, f1_counter = 0;
    for (long long unsigned x = 0; x < (uint)pow(2,n); x++) {
        if (scalar_product_unit_vector(x, ei) == 0) {
            f0_counter += f[x];
        } else {
            f1_counter += f[x];
        }
        *f_wt += f[x];
    }
    // printf("TT wt je %u, a f0 i f1 su %u i %u\n", f_wt, f0_counter, f1_counter);

}

uint smart_split_and_count(uint* TT, long long unsigned ei, uint n, long long unsigned length_in_bits)
{   
    uint sub_weights[2] = {0, 0}; // sub_weights[0] is f0 weight, sub_weights[1] is f1 weight
    uint bit = 0;
    uint scalar_product = 0;
    uint subfunction_condition = 0;
    // printf("pocinjemo smart split\n");
    for (long long unsigned x = 0; x < length_in_bits; x++) {
        bit = joanne_get_bit(TT, x);
        scalar_product = scalar_product_unit_vector(x, ei);
        sub_weights[!!scalar_product] += bit;
    }

    // printf(" f0, f1, su %u %u\n", sub_weights[0], sub_weights[1]);
    subfunction_condition = (sub_weights[0] == g_correct_sub_wt) || (sub_weights[1] == g_correct_sub_wt);
    return subfunction_condition;
}

uint popcount_split(uint* TT, long long unsigned ei, uint n, long long unsigned length_in_bits)
{
    uint subfunction_condition = 1; /*true*/
    uint subs[2] = {0,0};
    uint num_bytes = length_in_bits / 8 ;
    uint num_chunks = length_in_bits / ei;
    uint length_of_chunk = num_bytes / num_chunks;

    // printf("length in bits is %llu total length in bytes is %u length of chunks is %u\n", length_in_bits, num_bytes, length_of_chunk);

    for (int i = 0; i < num_chunks; i++) {
        subs[i%2] += count((char *)TT + i * length_of_chunk, length_of_chunk);
    }

    // printf("popcount f0 is %u and f1 is %u\n", subs[0], subs[1]);
    subfunction_condition = (subs[0] == g_correct_sub_wt) || (subs[1] == g_correct_sub_wt);

    return subfunction_condition;
}

uint combined_split(uint* TT, uint n, long long unsigned length_in_bits)
{
    uint TT_wt = count((char *)TT, length_in_bits/8);
    // printf("wt cele funkcije je %u\n", TT_wt);
    uint TT_wt_condition = (TT_wt == g_higher_correct_TT_wt) || (TT_wt == g_lower_correct_TT_wt);

    // commented out because we dont want to stop the split early; we are measuring worst case scenario.
    // if (!TT_wt_condition) {
    //     // printf("TT has bad weight\n");
    //     return 0; /*false*/ 
    // }

    uint result = 1; /*true*/
    for (long long unsigned ei = 1; ei < length_in_bits; ei = ei << 1){                                             /*          iii. if f=[f0|f1] has wt(f) = bla bla and wt(f0) or wt(f1) == bla bla FOR ALL ei*/
        if (ei < sizeof(long long unsigned) * 8) {
            // printf("iterative split for ei %u\n", ei);
            result = smart_split_and_count(TT, ei, n, length_in_bits);
        } else {
            // printf("popcount split for ei %u\n", ei);
            result = popcount_split(TT, ei, n, length_in_bits);
        }

        // commented out because we dont want to stop the split early; we are measuring worst case scenario.
        // if (!result) {
        //     return result;
        // }
    }
    return result;
}

void* popcount_split2(void *arg)
{
    // printf("negde u popcount split 2\n");
    split_thread_arg* split_args = (split_thread_arg*) arg;
    uint subs[2] = {0,0};
    uint num_bytes = split_args->length_in_bits / 8 ;
    uint num_chunks = split_args->length_in_bits / split_args->ei;
    uint length_of_chunk = num_bytes / num_chunks;

    // printf("length in bits is %llu total length in bytes is %u length of chunks is %u\n", split_args->length_in_bits, num_bytes, length_of_chunk);

    for (int i = 0; i < num_chunks; i++) {
        subs[i%2] += count((char *)(split_args->TT) + i * length_of_chunk, length_of_chunk);
    }

    // printf("popcount f0 is %u and f1 is %u\n", subs[0], subs[1]);
    split_args->split_result = (subs[0] == g_correct_sub_wt) || (subs[1] == g_correct_sub_wt);
    // printf("sacuvali smo zivu glavu\n");
    return NULL;
}

void* iterative_split(void* arg)
{   
    // printf("nedge uiterative split\n");
    split_thread_arg* split_args = (split_thread_arg*) arg;
    uint sub_weights[2] = {0, 0}; // sub_weights[0] is f0 weight, sub_weights[1] is f1 weight
    uint bit = 0;
    uint scalar_product = 0;
    // printf("pocinjemo smart split\n");
    for (long long unsigned x = 0; x < split_args->length_in_bits; x++) {
        bit = joanne_get_bit(split_args->TT, x);
        scalar_product = scalar_product_unit_vector(x, split_args->ei);
        sub_weights[!!scalar_product] += bit;
    }

    // printf(" f0, f1, su %u %u\n", sub_weights[0], sub_weights[1]);
    split_args->split_result = (sub_weights[0] == g_correct_sub_wt) || (sub_weights[1] == g_correct_sub_wt);
    // printf("preziveli smo iterative split\n");

    return NULL;
}

uint combined_split2(uint* TT, uint n, long long unsigned length_in_bits)
{
    // printf("pre TT wt\n");
    uint TT_wt = count((char *)TT, length_in_bits/8);
    uint TT_wt_condition = (TT_wt == g_higher_correct_TT_wt) || (TT_wt == g_lower_correct_TT_wt);
    if (!TT_wt_condition) {
        // printf("TT has bad weight\n");
        return 0; /*false*/
    }
    // printf("posle TT wt\n");

    /* A batch of threads is defined by number of logical cores, and
     * total number of threads needed is equal to n, so we divide n by number of cores
     * to get the number of batches where all cores will be used at once
     * Some threads might remain otside of a batch and they are executed afterwards - modulo batch*/

    uint complete_batch_size = getNumCores();
    /* I'm not going to choke devsrv so lets limit how many threads I'll try out at once*/
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
// printf("posle racunanja\n");
    pthread_t* threads = (pthread_t*)malloc(complete_batch_size * sizeof(pthread_t));
    split_thread_arg* split_args = (split_thread_arg*)malloc(complete_batch_size * sizeof(split_thread_arg));
    proforme_arg = (void**) malloc(complete_batch_size * sizeof(void*));
// printf("posle alociranja a pre batch operacija\n");
    long long unsigned ei = 1;

    for (int i = 0; i < num_batches; i++) {
        for (int j = 0; j < complete_batch_size; j++){
            split_args[j].ei = ei;
            split_args[j].length_in_bits = length_in_bits;
            split_args[j].n = n;
            split_args[j].TT = TT;
            if (ei < smallest_chunk_size) {
                pthread_create((threads + j), NULL, &iterative_split, (split_args + j));
                // iterative_split((void*)(split_args + j));
            } else {
                pthread_create((threads + j), NULL, &popcount_split2, (split_args + j));
                // popcount_split2((void*)(split_args + j));
            }
            // if (!split_args[j].split_result) {
            //     free(split_args);
            //     free(threads);
            //     free(proforme_arg);
            //     return 0; /*A bad split*/
            // }
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
// printf("posle batch operacija a pre modula operacija\n");
    /*Now handle rest of threads (if any) that didnt fill out a complete batch*/
    for (int i = 0; i < modulo_batch_size; i++) {
        split_args[i].ei = ei;
        split_args[i].length_in_bits = length_in_bits;
        split_args[i].n = n;
        split_args[i].TT = TT;
        if (ei < smallest_chunk_size) {
            pthread_create((threads + i), NULL, &iterative_split, (split_args + i));
            // iterative_split((void*)(split_args+i));
        } else {
            pthread_create((threads + i), NULL, &popcount_split2, (split_args + i));
            // popcount_split2((void*)(split_args+i));
        }
        ei = ei << 1;

        // if (!split_args[i].split_result) {
        //         free(split_args);
        //         free(threads);
        //         free(proforme_arg);
        //         return 0; /*A bad split*/
        // }
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

int main(int argc, char** argv) {

    struct timespec start_time, end_time, loop_start, loop_end;
    int seconds;
    long nanoseconds;

    uint num_execs = 10;

    if (argc > 1) {
        num_execs = atoi(argv[1]);
    }

    Matrix_A();
    Matrix_W();

    /* od 8 do 16 za naivni */
    printf("Izvrsavamo naivni ANFT\n");
    for (int n = 8; n <= 16; n+=2) {

        uint duzina = (uint)pow(2,n);
        uint *anf = malloc(duzina*sizeof(uint));
        uint *tt = malloc(duzina*sizeof(uint));

        memset(anf, 0, duzina/8);
        memset(tt, 0, duzina/8);

        /* formiramo kvadratnu bent funkciju za ulaz i */
        for (int j = 0; j <= (n/2 - 1); j++) {
            /*      (a) Let j = 2^i + 2^(n/2+i) */
            long long unsigned k = (uint)pow(2,j) + (uint)pow(2, n/2+j);
            /*      (b) Set A_j = 1 */
            anf[k] = 1;
        }

        /* pocni da meris vreme*/
        clock_gettime(CLOCK_REALTIME, &start_time);
        /*izvrsavanja*/
        for(int j = 0; j < num_execs; j++) {
            anft_naive(anf, tt, n);
        }
        clock_gettime(CLOCK_REALTIME, &end_time);
        seconds = end_time.tv_sec - start_time.tv_sec;
        nanoseconds = 0;
        if ((end_time.tv_nsec - start_time.tv_nsec) < 0) {
            nanoseconds = 1000000000 + end_time.tv_nsec - start_time.tv_nsec;
            seconds--;
        } else {
            nanoseconds = end_time.tv_nsec - start_time.tv_nsec; 
        }
        printf("%d izsvrsavanja naivnog ANFT za stepen %u je trajalao %d seconds, %ldms, %ldus, %ldns\n", num_execs, n, seconds, nanoseconds/1000000, (nanoseconds/1000)%1000, nanoseconds%1000);
        free(anf);
        free(tt);
    }


/////////////////////////////////////////////////////////////////////////////////


    printf("\nPrelazimo na optimizovani anft\n");
    /* od 8 do 16 za  x64 */
    for (int n = 8; n <= 16; n+=2) {

        long long unsigned duzina = (long long unsigned)pow(2,n);
        uint bits_per_uint = sizeof(uint) * 8;
        uint *anf = malloc(duzina/8);
        uint *tt = malloc(duzina/8);

        memset(anf, 0, duzina/8);
        memset(tt, 0, duzina/8);

        /* formiramo kvadratnu bent funkciju za ulaz i */
        for (int j = 0; j <= (n/2 - 1); j++) {
            /*      (a) Let j = 2^i + 2^(n/2+i) */
            long long unsigned k = (uint)pow(2,j) + (uint)pow(2, n/2+j);
            /*      (b) Set A_j = 1 */
            joanne_set_bit(anf, k, 1);
        }

        /*izvrsavanja*/
        clock_gettime(CLOCK_REALTIME, &start_time);
        for(int j = 0; j < num_execs; j++) {
            OPTIMIZED_ANF((ull*)anf, n);
        }
        clock_gettime(CLOCK_REALTIME, &end_time);

        seconds = end_time.tv_sec - start_time.tv_sec;
        nanoseconds = 0;
        if ((end_time.tv_nsec - start_time.tv_nsec) < 0) {
            nanoseconds = 1000000000 + end_time.tv_nsec - start_time.tv_nsec;
            seconds--;
        } else {
            nanoseconds = end_time.tv_nsec - start_time.tv_nsec; 
        }
        printf("%d izsvrsavanja optimizovanog ANFT za stepen %u je trajalao %d seconds, %ldms, %ldus, %ldns\n", num_execs, n, seconds, nanoseconds/1000000, (nanoseconds/1000)%1000, nanoseconds%1000);
        free(anf);
        free(tt);
    }
//////////////////////////////////////////////////

printf("\nPrelazimo na optimizovani anft od 18 do 32\n");
    /* od 8 do 16 za x64 */
    for (int n = 18; n <= 32; n+=2) {

        long long unsigned duzina = (long long unsigned)pow(2,n);
        uint bits_per_uint = sizeof(uint) * 8;
        uint *anf = malloc(duzina/8);
        uint *tt = malloc(duzina/8);

        memset(anf, 0, duzina/8);
        memset(tt, 0, duzina/8);

        /* formiramo kvadratnu bent funkciju za ulaz i */
        for (int j = 0; j <= (n/2 - 1); j++) {
            /*      (a) Let j = 2^i + 2^(n/2+i) */
            long long unsigned k = (uint)pow(2,j) + (uint)pow(2, n/2+j);
            /*      (b) Set A_j = 1 */
            joanne_set_bit(anf, k, 1);
        }

        /*izvrsavanja*/
        clock_gettime(CLOCK_REALTIME, &start_time);
        for(int j = 0; j < num_execs; j++) {
            OPTIMIZED_ANF((ull*)anf, n);
        }
        clock_gettime(CLOCK_REALTIME, &end_time);

        seconds = end_time.tv_sec - start_time.tv_sec;
        nanoseconds = 0;
        if ((end_time.tv_nsec - start_time.tv_nsec) < 0) {
            nanoseconds = 1000000000 + end_time.tv_nsec - start_time.tv_nsec;
            seconds--;
        } else {
            nanoseconds = end_time.tv_nsec - start_time.tv_nsec; 
        }
        printf("%d izsvrsavanja optimizovanog ANFT za stepen %u je trajalao %d seconds, %ldms, %ldus, %ldns\n", num_execs, n, seconds, nanoseconds/1000000, (nanoseconds/1000)%1000, nanoseconds%1000);
        free(anf);
        free(tt);
    }

/////////////////////////////////////////////////////////////////////////////
    printf("\nNaivni WHT\n");

    /* od 8 do 16 za naivni */
    for (int n = 8; n <= 16; n+=2) {

        uint duzina = (uint)pow(2,n);
        uint *anf = malloc(duzina*sizeof(uint));
        uint *tt = malloc(duzina*sizeof(uint));
        int *wht = malloc(duzina*sizeof(int));

        memset(anf, 0, duzina*sizeof(uint));
        memset(tt, 0, duzina*sizeof(uint));
        memset(wht, 0, duzina*sizeof(int));

        /* formiramo kvadratnu bent funkciju za ulaz i */
        for (int j = 0; j <= (n/2 - 1); j++) {
            /*      (a) Let j = 2^i + 2^(n/2+i) */
            long long unsigned k = (uint)pow(2,j) + (uint)pow(2, n/2+j);
            /*      (b) Set A_j = 1 */
            anf[k] = 1;
        }

        anft_naive(anf, tt, n); // moramo bar jedan ANFT da uradimo da bi mogli na necemu WHT
        /* pocni da meris vreme*/
        clock_gettime(CLOCK_REALTIME, &start_time);
        /*izvrsavanja*/
        for(int j = 0; j < num_execs; j++) {
            generate_wht(anf, tt, NULL, n);
        }
        clock_gettime(CLOCK_REALTIME, &end_time);
        seconds = end_time.tv_sec - start_time.tv_sec;
        nanoseconds = 0;
        if ((end_time.tv_nsec - start_time.tv_nsec) < 0) {
            nanoseconds = 1000000000 + end_time.tv_nsec - start_time.tv_nsec;
            seconds--;
        } else {
            nanoseconds = end_time.tv_nsec - start_time.tv_nsec; 
        }
        printf("%d izsvrsavanja naivnog WHT za stepen %u je trajalao %d seconds, %ldms, %ldus, %ldns\n", num_execs, n, seconds, nanoseconds/1000000, (nanoseconds/1000)%1000, nanoseconds%1000);
        free(anf);
        free(tt);
        free(wht);
    }

/////////////////////////////////////////////////////////////////////////////////////
    printf("\nBrzi WHt\n");
    /* od 8 do 16 za x64 */
    for (int n = 8; n <= 16; n+=2) {

        long long unsigned duzina = (long long unsigned)pow(2,n);
        uint bits_per_uint = sizeof(uint) * 8;
        uint *anf = malloc(duzina/8);
        uint *tt = malloc(duzina/8);
        int *wht = malloc(duzina*sizeof(int));

        memset(anf, 0, duzina/8);
        memset(tt, 0, duzina/8);
        memset(wht, 0, duzina*sizeof(int));

        /* formiramo kvadratnu bent funkciju za ulaz i */
        // printf("formiramo kvadratnu bent stepena %d\n", n);
        for (int j = 0; j <= (n/2 - 1); j++) {
            /*      (a) Let j = 2^i + 2^(n/2+i) */
            long long unsigned k = (uint)pow(2,j) + (uint)pow(2, n/2+j);
            /*      (b) Set A_j = 1 */
            // printf("setujemo %llu-ti bit\n", k);
            joanne_set_bit(anf, k, 1);
        }

        ANFT(anf, tt, n); // moramo jedan anft da uradimo da ima nad cime WHT

        /*izvrsavanja*/
        clock_gettime(CLOCK_REALTIME, &start_time);
        for(int j = 0; j < num_execs; j++) {
            WHT(tt, wht, n);
        }
        clock_gettime(CLOCK_REALTIME, &end_time);

        seconds = end_time.tv_sec - start_time.tv_sec;
        nanoseconds = 0;
        if ((end_time.tv_nsec - start_time.tv_nsec) < 0) {
            nanoseconds = 1000000000 + end_time.tv_nsec - start_time.tv_nsec;
            seconds--;
        } else {
            nanoseconds = end_time.tv_nsec - start_time.tv_nsec; 
        }
        printf("%d izsvrsavanja obicnog brzog ANFT za stepen %u je trajalao %d seconds, %ldms, %ldus, %ldns\n", num_execs, n, seconds, nanoseconds/1000000, (nanoseconds/1000)%1000, nanoseconds%1000);
        free(anf);
        free(tt);
        free(wht);
    }

//////////////////////////////////////////////////////////////////////////////

    printf("\nBrzi tredovan WHt\n");
    /* od 8 do 16 za x64 */
    for (int n = 18; n <= 32; n+=2) {

        long long unsigned duzina = (long long unsigned)pow(2,n);
        uint bits_per_uint = sizeof(uint) * 8;
        uint *anf = malloc(duzina/8);
        uint *tt = malloc(duzina/8);
        int *wht = malloc(duzina*sizeof(int));

        memset(anf, 0, duzina/8);
        memset(tt, 0, duzina/8);
        memset(wht, 0, duzina*sizeof(int));

        /* formiramo kvadratnu bent funkciju za ulaz i */
        // printf("formiramo kvadratnu bent stepena %d\n", n);
        for (int j = 0; j <= (n/2 - 1); j++) {
            /*      (a) Let j = 2^i + 2^(n/2+i) */
            long long unsigned k = (uint)pow(2,j) + (uint)pow(2, n/2+j);
            /*      (b) Set A_j = 1 */
            // printf("setujemo %llu-ti bit\n", k);
            joanne_set_bit(anf, k, 1);
        }

        ANFT(anf, tt, n); // moramo jedan anft da uradimo da ima nad cime WHT

        /*izvrsavanja*/
        clock_gettime(CLOCK_REALTIME, &start_time);
        for(int j = 0; j < num_execs; j++) {
            fully_threaded_WHT(tt, wht, n);
        }
        clock_gettime(CLOCK_REALTIME, &end_time);

        seconds = end_time.tv_sec - start_time.tv_sec;
        nanoseconds = 0;
        if ((end_time.tv_nsec - start_time.tv_nsec) < 0) {
            nanoseconds = 1000000000 + end_time.tv_nsec - start_time.tv_nsec;
            seconds--;
        } else {
            nanoseconds = end_time.tv_nsec - start_time.tv_nsec; 
        }
        printf("%d izsvrsavanja obicnog brzog ANFT za stepen %u je trajalao %d seconds, %ldms, %ldus, %ldns\n", num_execs, n, seconds, nanoseconds/1000000, (nanoseconds/1000)%1000, nanoseconds%1000);
        free(anf);
        free(tt);
        free(wht);
    }

/////////////////////////////////////////////////////////////////////////////////////////////

    printf("\nNaivni split\n");

    /* od 8 do 16 za naivni */
    for (int n = 8; n <= 16; n+=2) {

        uint duzina = (uint)pow(2,n);
        uint *anf = malloc(duzina*sizeof(uint));
        uint *tt = malloc(duzina*sizeof(uint));
        uint f_wt = 0;

        memset(anf, 0, duzina*sizeof(uint));
        memset(tt, 0, duzina*sizeof(uint));

        g_bent_nonlinearity = (uint)pow(2, n-1) - (uint)pow(2, n/2 - 1);
        g_higher_correct_TT_wt = (uint)pow(2,n-1) + (uint)pow(2,n/2 - 1);
        g_lower_correct_TT_wt = (uint)pow(2,n-1) - (uint)pow(2,n/2 - 1);
        g_correct_sub_wt = (uint)pow(2, n-2);

        /* formiramo kvadratnu bent funkciju za ulaz i */
        for (int j = 0; j <= (n/2 - 1); j++) {
            /*      (a) Let j = 2^i + 2^(n/2+i) */
            long long unsigned k = (uint)pow(2,j) + (uint)pow(2, n/2+j);
            /*      (b) Set A_j = 1 */
            anf[k] = 1;
        }

        anft_naive(anf, tt, n);
        /* pocni da meris vreme*/
        clock_gettime(CLOCK_REALTIME, &start_time);
        /*izvrsavanja*/
        for(int j = 0; j < num_execs; j++) {
            for (uint ei = 1; ei < duzina; ei = ei << 1) {
                split_along_direction_e_n(tt, n, &f_wt, ei);
            }
        }
        clock_gettime(CLOCK_REALTIME, &end_time);
        printf("printom protiv optimizacije %u\n", f_wt);
        seconds = end_time.tv_sec - start_time.tv_sec;
        nanoseconds = 0;
        if ((end_time.tv_nsec - start_time.tv_nsec) < 0) {
            nanoseconds = 1000000000 + end_time.tv_nsec - start_time.tv_nsec;
            seconds--;
        } else {
            nanoseconds = end_time.tv_nsec - start_time.tv_nsec; 
        }
        printf("%d izsvrsavanja naivnog splita za stepen %u je trajalao %d seconds, %ldms, %ldus, %ldns\n", num_execs, n, seconds, nanoseconds/1000000, (nanoseconds/1000)%1000, nanoseconds%1000);
        free(anf);
        free(tt);

    }

/////////////////////////////////////////////////////////////////////////////////////
    printf("\nBrzi Split\n");
    /* od 8 do 16 za  x64 */
    for (int n = 8; n <= 16; n+=2) {

        long long unsigned duzina = (long long unsigned)pow(2,n);
        uint bits_per_uint = sizeof(uint) * 8;
        uint sp = 0;
        uint *anf = malloc(duzina/8);
        uint *tt = malloc(duzina/8);

        g_bent_nonlinearity = (uint)pow(2, n-1) - (uint)pow(2, n/2 - 1);
        g_higher_correct_TT_wt = (uint)pow(2,n-1) + (uint)pow(2,n/2 - 1);
        g_lower_correct_TT_wt = (uint)pow(2,n-1) - (uint)pow(2,n/2 - 1);
        g_correct_sub_wt = (uint)pow(2, n-2);

        memset(anf, 0, duzina/8);
        memset(tt, 0, duzina/8);

        /* formiramo kvadratnu bent funkciju za ulaz i */
        for (int j = 0; j <= (n/2 - 1); j++) {
            /*      (a) Let j = 2^i + 2^(n/2+i) */
            long long unsigned k = (uint)pow(2,j) + (uint)pow(2, n/2+j);
            /*      (b) Set A_j = 1 */
            joanne_set_bit(anf, k, 1);
        }

        OPTIMIZED_ANF((ull*)anf, n); // moramo jedan anft da uradimo da ima nad cime split

        /*izvrsavanja*/
        clock_gettime(CLOCK_REALTIME, &start_time);
        for(int j = 0; j < num_execs; j++) {
            sp = combined_split(anf, n, duzina);
        }
        clock_gettime(CLOCK_REALTIME, &end_time);
        printf("ovo je da mozak ne bi izoptimizovao moj split %u\n", sp);

        seconds = end_time.tv_sec - start_time.tv_sec;
        nanoseconds = 0;
        if ((end_time.tv_nsec - start_time.tv_nsec) < 0) {
            nanoseconds = 1000000000 + end_time.tv_nsec - start_time.tv_nsec;
            seconds--;
        } else {
            nanoseconds = end_time.tv_nsec - start_time.tv_nsec; 
        }
        printf("%d izsvrsavanja brzog splita za stepen %u je trajalao %d seconds, %ldms, %ldus, %ldns\n", num_execs, n, seconds, nanoseconds/1000000, (nanoseconds/1000)%1000, nanoseconds%1000);
        free(anf);
        free(tt);
    }

    printf("\nBrzi tredovani Split\n");
    /* od 8 do 16 za  x64 */
    for (int n = 18; n <= 32; n+=2) {

        long long unsigned duzina = (long long unsigned)pow(2,n);
        uint bits_per_uint = sizeof(uint) * 8;
        uint sp = 0;
        uint *anf = malloc(duzina/8);
        uint *tt = malloc(duzina/8);

        g_bent_nonlinearity = (uint)pow(2, n-1) - (uint)pow(2, n/2 - 1);
        g_higher_correct_TT_wt = (uint)pow(2,n-1) + (uint)pow(2,n/2 - 1);
        g_lower_correct_TT_wt = (uint)pow(2,n-1) - (uint)pow(2,n/2 - 1);
        g_correct_sub_wt = (uint)pow(2, n-2);

        memset(anf, 0, duzina/8);
        memset(tt, 0, duzina/8);

        /* formiramo kvadratnu bent funkciju za ulaz i */
        for (int j = 0; j <= (n/2 - 1); j++) {
            /*      (a) Let j = 2^i + 2^(n/2+i) */
            long long unsigned k = (uint)pow(2,j) + (uint)pow(2, n/2+j);
            /*      (b) Set A_j = 1 */
            joanne_set_bit(anf, k, 1);
        }

        OPTIMIZED_ANF((ull*)anf, n); // moramo jedan anft da uradimo da ima nad cime split

        /*izvrsavanja*/
        clock_gettime(CLOCK_REALTIME, &start_time);
        for(int j = 0; j < num_execs; j++) {
            sp = combined_split2(anf, n, duzina);
        }
        clock_gettime(CLOCK_REALTIME, &end_time);
        printf("ovo je da mozak ne bi izoptimizovao moj split %u\n", sp);

        seconds = end_time.tv_sec - start_time.tv_sec;
        nanoseconds = 0;
        if ((end_time.tv_nsec - start_time.tv_nsec) < 0) {
            nanoseconds = 1000000000 + end_time.tv_nsec - start_time.tv_nsec;
            seconds--;
        } else {
            nanoseconds = end_time.tv_nsec - start_time.tv_nsec; 
        }
        printf("%d izsvrsavanja brzog splita za stepen %u je trajalao %d seconds, %ldms, %ldus, %ldns\n", num_execs, n, seconds, nanoseconds/1000000, (nanoseconds/1000)%1000, nanoseconds%1000);
        free(anf);
        free(tt);
    }


    return 0;
}
