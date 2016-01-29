// (c) Copyright 2013 Lev Barash, Landau Institute for Theoretical Physics, Russian Academy of Sciences
// This is supplement to the paper:
// L.Yu. Barash, L.N. Shchur, "PRAND: GPU accelerated parallel random number generation library: Using most reliable algorithms and applying parallelism of modern GPUs and CPUs".
// e-mail: barash @ itp.ac.ru (remove space)

typedef unsigned long long lt;

typedef struct{
  unsigned z1,z2,z3,z4;
} lfsr113_state;

typedef struct{
  unsigned z[4] __attribute__ ((aligned(16)));
} lfsr113_sse_state;

__device__ __host__ void lfsr113_skipahead_(lfsr113_state* state,unsigned long long offset64,unsigned long long offset0);

__device__ __host__ void lfsr113_init_(lfsr113_state* state);

__device__ __host__ void lfsr113_init_sequence_(lfsr113_state* state,lt SequenceNumber);

__device__ __host__ void lfsr113_init_long_sequence_(lfsr113_state* state,lt SequenceNumber);

__device__ __host__ unsigned int lfsr113_generate_(lfsr113_state* state);

__device__ __host__ float lfsr113_generate_uniform_float_(lfsr113_state* state);

__host__ unsigned int lfsr113_sse_generate_(lfsr113_sse_state* state); // use this function only if CPU supports SSE4

__device__ __host__ void lfsr113_get_sse_state_(lfsr113_state* state,lfsr113_sse_state* sse_state);

__host__ void lfsr113_print_state_(lfsr113_state* state);

__host__ void lfsr113_generate_array_(lfsr113_state* state, unsigned int* out, long length);

__host__ void lfsr113_generate_gpu_array_(lfsr113_state* state, unsigned int* out, long length);

__host__ void lfsr113_generate_gpu_array_float_(lfsr113_state* state, float* out, long length);

__host__ void lfsr113_generate_gpu_array_double_(lfsr113_state* state, double* out, long length);
