// (c) Copyright 2013 Lev Barash, Landau Institute for Theoretical Physics, Russian Academy of Sciences
// This is supplement to the paper:
// L.Yu. Barash, L.N. Shchur, "PRAND: GPU accelerated parallel random number generation library: Using most reliable algorithms and applying parallelism of modern GPUs and CPUs".
// e-mail: barash @ itp.ac.ru (remove space)

#define gq58x4_g       288230374541099008ULL

typedef unsigned long long lt;

typedef struct{
  lt xN[8] __attribute__ ((aligned(16))),
     xP[8] __attribute__ ((aligned(16)));
} gq58x4_state;

typedef gq58x4_state gq58x4_sse_state;

__device__ __host__ void gq58x4_skipahead_(gq58x4_state* state, lt offset64, lt offset0);

__device__ __host__ void gq58x4_init_(gq58x4_state* state);

__device__ __host__ void gq58x4_init_short_sequence_(gq58x4_state* state,unsigned SequenceNumber);

__device__ __host__ void gq58x4_init_medium_sequence_(gq58x4_state* state,unsigned SequenceNumber);

__device__ __host__ void gq58x4_init_long_sequence_(gq58x4_state* state,unsigned SequenceNumber);

__device__ __host__ unsigned int gq58x4_generate_(gq58x4_state* state);

__device__ __host__ float gq58x4_generate_uniform_float_(gq58x4_state* state);

__host__ unsigned int gq58x4_sse_generate_(gq58x4_sse_state* state);

__device__ __host__ void gq58x4_get_sse_state_(gq58x4_state* state,gq58x4_sse_state* sse_state);

__host__ void gq58x4_print_state_(gq58x4_state* state);

__host__ void gq58x4_print_sse_state_(gq58x4_sse_state* state);

__host__ void gq58x4_generate_array_(gq58x4_state* state, unsigned int* out, long length);

__host__ void gq58x4_generate_gpu_array_(gq58x4_state* state, unsigned int* out, long length);

__host__ void gq58x4_generate_gpu_array_float_(gq58x4_state* state, float* out, long length);

__host__ void gq58x4_generate_gpu_array_double_(gq58x4_state* state, double* out, long length);
