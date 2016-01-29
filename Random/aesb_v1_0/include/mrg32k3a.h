// (c) Copyright 2013 Lev Barash, Landau Institute for Theoretical Physics, Russian Academy of Sciences
// This is supplement to the paper:
// L.Yu. Barash, L.N. Shchur, "PRAND: GPU accelerated parallel random number generation library: Using most reliable algorithms and applying parallelism of modern GPUs and CPUs".
// e-mail: barash @ itp.ac.ru (remove space)

typedef struct{
  unsigned x0,x1,x2,y0,y1,y2;
} mrg32k3a_state;

typedef struct{
  unsigned s[20] __attribute__ ((aligned(16))); 
} mrg32k3a_sse_state;

__host__ void mrg32k3a_init_device_consts_();

__host__ void mrg32k3a_free_device_consts_();

__host__ void mrg32k3a_skipahead_(mrg32k3a_state* state,unsigned long long offset128,
                                       unsigned long long offset64,unsigned long long offset0);

__device__ void mrg32k3a_dev_skipahead_(mrg32k3a_state* state,unsigned long long offset128,
                                        unsigned long long offset64,unsigned long long offset0);

__device__ __host__ void mrg32k3a_init_(mrg32k3a_state* state);

__host__ void mrg32k3a_init_sequence_(mrg32k3a_state* state,unsigned long long SequenceNumber);

__device__ void mrg32k3a_dev_init_sequence_(mrg32k3a_state* state,unsigned long long SequenceNumber);

__device__ __host__ unsigned int mrg32k3a_generate_(mrg32k3a_state* state);

__device__ __host__ float mrg32k3a_generate_uniform_float_(mrg32k3a_state* state);

__host__ unsigned int mrg32k3a_sse_generate_(mrg32k3a_sse_state* state);

__device__ __host__ void mrg32k3a_get_sse_state_(mrg32k3a_state* state,mrg32k3a_sse_state* sse_state);

__host__ void mrg32k3a_print_state_(mrg32k3a_state* state);

__host__ void mrg32k3a_generate_array_(mrg32k3a_state* state, unsigned int* out, long length);

__host__ void mrg32k3a_generate_gpu_array_(mrg32k3a_state* state, unsigned int* out, long length);

__host__ void mrg32k3a_generate_gpu_array_float_(mrg32k3a_state* state, float* out, long length);

__host__ void mrg32k3a_generate_gpu_array_double_(mrg32k3a_state* state, double* out, long length);
