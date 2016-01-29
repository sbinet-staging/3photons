// (c) Copyright 2013 Lev Barash, Landau Institute for Theoretical Physics, Russian Academy of Sciences
// This is supplement to the paper:
// L.Yu. Barash, L.N. Shchur, "PRAND: GPU accelerated parallel random number generation library: Using most reliable algorithms and applying parallelism of modern GPUs and CPUs".
// e-mail: barash @ itp.ac.ru (remove space)

#define gm55_g       36028797018961904ULL

typedef unsigned long long lt;

typedef struct{
  lt xN[8] __attribute__ ((aligned(16))), 
     xP[8] __attribute__ ((aligned(16)));
} gm55_state;

typedef gm55_state gm55_sse_state;

__device__ __host__ void gm55_skipahead_(gm55_state* state, lt offset64, lt offset0);

__device__ __host__ void gm55_init_(gm55_state* state);

__device__ __host__ void gm55_init_short_sequence_(gm55_state* state,lt SequenceNumber);

__device__ __host__ void gm55_init_long_sequence_(gm55_state* state,lt SequenceNumber);

__device__ __host__ unsigned int gm55_generate_(gm55_state* state);

__device__ __host__ float gm55_generate_uniform_float_(gm55_state* state);

__host__ unsigned gm55_sse_generate_(gm55_sse_state* state);

__device__ __host__ void gm55_get_sse_state_(gm55_state* state,gm55_sse_state* sse_state);

__host__ void gm55_print_state_(gm55_state* state);

__host__ void gm55_print_sse_state_(gm55_sse_state* state);

__host__ void gm55_generate_array_(gm55_state* state, unsigned int* out, long length);

__host__ void gm55_generate_gpu_array_(gm55_state* state, unsigned int* out, long length);

__host__ void gm55_generate_gpu_array_float_(gm55_state* state, float* out, long length);

__host__ void gm55_generate_gpu_array_double_(gm55_state* state, double* out, long length);
