// (c) Copyright 2013 Lev Barash, Landau Institute for Theoretical Physics, Russian Academy of Sciences
// This is supplement to the paper:
// L.Yu. Barash, L.N. Shchur, "PRAND: GPU accelerated parallel random number generation library: Using most reliable algorithms and applying parallelism of modern GPUs and CPUs".
// e-mail: barash @ itp.ac.ru (remove space)

#define gm29_g 536870909U

typedef struct{
  unsigned xN[32] __attribute__ ((aligned(16))), 
           xP[32] __attribute__ ((aligned(16)));
} gm29_state;

typedef gm29_state gm29_sse_state;

__device__ __host__ void gm29_skipahead_(gm29_state* state, unsigned long long offset);

__device__ __host__ void gm29_init_(gm29_state* state);

__device__ __host__ void gm29_init_short_sequence_(gm29_state* state,unsigned SequenceNumber);

__device__ __host__ void gm29_init_medium_sequence_(gm29_state* state,unsigned SequenceNumber);

__device__ __host__ void gm29_init_long_sequence_(gm29_state* state,unsigned SequenceNumber);

__device__ __host__ unsigned int gm29_generate_(gm29_state* state);

__device__ __host__ float gm29_generate_uniform_float_(gm29_state* state);

__host__ unsigned int gm29_sse_generate_(gm29_sse_state* state);

__device__ __host__ void gm29_get_sse_state_(gm29_state* state,gm29_sse_state* sse_state);

__host__ void gm29_print_state_(gm29_state* state);

__host__ void gm29_print_sse_state_(gm29_sse_state* state);

__host__ void gm29_generate_array_(gm29_state* state, unsigned int* out, long length);

__host__ void gm29_generate_gpu_array_(gm29_state* state, unsigned int* out, long length);

__host__ void gm29_generate_gpu_array_float_(gm29_state* state, float* out, long length);

__host__ void gm29_generate_gpu_array_double_(gm29_state* state, double* out, long length);
