// (c) Copyright 2013 Lev Barash, Landau Institute for Theoretical Physics, Russian Academy of Sciences
// This is supplement to the paper:
// L.Yu. Barash, L.N. Shchur, "PRAND: GPU accelerated parallel random number generation library: Using most reliable algorithms and applying parallelism of modern GPUs and CPUs".
// e-mail: barash @ itp.ac.ru (remove space)

#define gm19_g 524287

typedef struct{
  unsigned xN[32] __attribute__ ((aligned(16))),
           xP[32] __attribute__ ((aligned(16)));
} gm19_state;

typedef gm19_state gm19_sse_state;

__device__ __host__ void gm19_skipahead_(gm19_state* state, unsigned long long offset);

__device__ __host__ void gm19_init_(gm19_state* state);

__device__ __host__ void gm19_init_sequence_(gm19_state* state,unsigned SequenceNumber);

__device__ __host__ unsigned int gm19_generate_(gm19_state* state);

__device__ __host__ float gm19_generate_uniform_float_(gm19_state* state);

__host__ unsigned int gm19_sse_generate_(gm19_sse_state* state);

__device__ __host__ void gm19_get_sse_state_(gm19_state* state,gm19_sse_state* sse_state);

__host__ void gm19_print_state_(gm19_state* state);

__host__ void gm19_print_sse_state_(gm19_sse_state* state);

__host__ void gm19_generate_array_(gm19_state* state, unsigned int* out, long length);

__host__ void gm19_generate_gpu_array_(gm19_state* state, unsigned int* out, long length);

__host__ void gm19_generate_gpu_array_float_(gm19_state* state, float* out, long length);

__host__ void gm19_generate_gpu_array_double_(gm19_state* state, double* out, long length);
