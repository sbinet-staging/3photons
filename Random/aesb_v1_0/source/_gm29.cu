// (c) Copyright 2013 Lev Barash, Landau Institute for Theoretical Physics, Russian Academy of Sciences
// This is supplement to the paper:
// L.Yu. Barash, L.N. Shchur, "PRAND: GPU accelerated parallel random number generation library: Using most reliable algorithms and applying parallelism of modern GPUs and CPUs".
// e-mail: barash @ itp.ac.ru (remove space)

#include<stdio.h>

#define gm29_CUDA_CALL(x) do { if((x) != cudaSuccess) { printf("Error: %s at %s:%d\n",cudaGetErrorString(cudaGetLastError()),__FILE__,__LINE__); exit(1);}} while(0)

#define gm29_BLOCKS  512
#define gm29_THREADS 128
#define gm29_ARRAY_SECTIONS (gm29_BLOCKS*gm29_THREADS/32)

#define gm29_k 4
#define gm29_q 2
#define gm29_g 536870909U
#define gm29_halfg 268435456U

typedef struct{
  unsigned xN[32] __attribute__ ((aligned(16))), 
           xP[32] __attribute__ ((aligned(16)));
} gm29_state;

typedef gm29_state gm29_sse_state;

unsigned gm29_sse_Consts[16] __attribute__ ((aligned(16))) = 
{536870911,536870911,536870911,536870911,1073741818,1073741818,1073741818,1073741818,
 536870908,536870908,536870908,536870908,536870909,536870909,536870909,536870909};

__host__ unsigned int gm29_sse_generate_(gm29_sse_state* state){
  unsigned output1; unsigned output2 __attribute__ ((unused));
  asm volatile("movaps (%4),%%xmm7\n" \
      "movaps 16(%4),%%xmm6\n" \
      "movaps 32(%4),%%xmm4\n" \

      "movaps (%2),%%xmm0\n" \
      "movaps (%3),%%xmm5\n" \
      "movaps %%xmm0,(%3)\n" \
      "pslld  $2,%%xmm0\n" \
      "paddd  %%xmm6,%%xmm0\n" \
      "pslld  $1,%%xmm5\n" \
      "psubd  %%xmm5,%%xmm0\n" \
      "movaps %%xmm0,%%xmm5\n" \
      "psrld  $29,%%xmm5\n" \
      "pand   %%xmm7,%%xmm0\n" \
      "paddd  %%xmm5,%%xmm0\n" \
      "paddd  %%xmm5,%%xmm0\n" \
      "paddd  %%xmm5,%%xmm0\n" \
      "movaps %%xmm0,%%xmm5\n" \
      "pcmpgtd  %%xmm4,%%xmm5\n" \
      "pand   48(%4),%%xmm5\n" \
      "psubd  %%xmm5,%%xmm0\n" \
      "movaps %%xmm0,(%2)\n" \

      "movaps 16(%2),%%xmm1\n" \
      "movaps 16(%3),%%xmm5\n" \
      "movaps %%xmm1,16(%3)\n" \
      "pslld  $2,%%xmm1\n" \
      "paddd  %%xmm6,%%xmm1\n" \
      "pslld  $1,%%xmm5\n" \
      "psubd  %%xmm5,%%xmm1\n" \
      "movaps %%xmm1,%%xmm5\n" \
      "psrld  $29,%%xmm5\n" \
      "pand   %%xmm7,%%xmm1\n" \
      "paddd  %%xmm5,%%xmm1\n" \
      "paddd  %%xmm5,%%xmm1\n" \
      "paddd  %%xmm5,%%xmm1\n" \
      "movaps %%xmm1,%%xmm5\n" \
      "pcmpgtd  %%xmm4,%%xmm5\n" \
      "pand   48(%4),%%xmm5\n" \
      "psubd  %%xmm5,%%xmm1\n" \
      "movaps %%xmm1,16(%2)\n" \

      "movaps 32(%2),%%xmm2\n" \
      "movaps 32(%3),%%xmm5\n" \
      "movaps %%xmm2,32(%3)\n" \
      "pslld  $2,%%xmm2\n" \
      "paddd  %%xmm6,%%xmm2\n" \
      "pslld  $1,%%xmm5\n" \
      "psubd  %%xmm5,%%xmm2\n" \
      "movaps %%xmm2,%%xmm5\n" \
      "psrld  $29,%%xmm5\n" \
      "pand   %%xmm7,%%xmm2\n" \
      "paddd  %%xmm5,%%xmm2\n" \
      "paddd  %%xmm5,%%xmm2\n" \
      "paddd  %%xmm5,%%xmm2\n" \
      "movaps %%xmm2,%%xmm5\n" \
      "pcmpgtd  %%xmm4,%%xmm5\n" \
      "pand   48(%4),%%xmm5\n" \
      "psubd  %%xmm5,%%xmm2\n" \
      "movaps %%xmm2,32(%2)\n" \

      "movaps 48(%2),%%xmm3\n" \
      "movaps 48(%3),%%xmm5\n" \
      "movaps %%xmm3,48(%3)\n" \
      "pslld  $2,%%xmm3\n" \
      "paddd  %%xmm6,%%xmm3\n" \
      "pslld  $1,%%xmm5\n" \
      "psubd  %%xmm5,%%xmm3\n" \
      "movaps %%xmm3,%%xmm5\n" \
      "psrld  $29,%%xmm5\n" \
      "pand   %%xmm7,%%xmm3\n" \
      "paddd  %%xmm5,%%xmm3\n" \
      "paddd  %%xmm5,%%xmm3\n" \
      "paddd  %%xmm5,%%xmm3\n" \
      "movaps %%xmm3,%%xmm5\n" \
      "pcmpgtd  %%xmm4,%%xmm5\n" \
      "pand   48(%4),%%xmm5\n" \
      "psubd  %%xmm5,%%xmm3\n" \
      "movaps %%xmm3,48(%2)\n" \


      "psrld  $28,%%xmm0\n" \
      "psrld  $28,%%xmm1\n" \
      "psrld  $28,%%xmm2\n" \
      "psrld  $28,%%xmm3\n" \
      "packssdw %%xmm1,%%xmm0\n" \
      "packssdw %%xmm3,%%xmm2\n" \
      "packsswb %%xmm2,%%xmm0\n" \
      "psllw  $7,%%xmm0\n" \
      "pmovmskb %%xmm0,%0\n" \

      "movaps 64(%2),%%xmm0\n" \
      "movaps 64(%3),%%xmm5\n" \
      "movaps %%xmm0,64(%3)\n" \
      "pslld  $2,%%xmm0\n" \
      "paddd  %%xmm6,%%xmm0\n" \
      "pslld  $1,%%xmm5\n" \
      "psubd  %%xmm5,%%xmm0\n" \
      "movaps %%xmm0,%%xmm5\n" \
      "psrld  $29,%%xmm5\n" \
      "pand   %%xmm7,%%xmm0\n" \
      "paddd  %%xmm5,%%xmm0\n" \
      "paddd  %%xmm5,%%xmm0\n" \
      "paddd  %%xmm5,%%xmm0\n" \
      "movaps %%xmm0,%%xmm5\n" \
      "pcmpgtd  %%xmm4,%%xmm5\n" \
      "pand   48(%4),%%xmm5\n" \
      "psubd  %%xmm5,%%xmm0\n" \
      "movaps %%xmm0,64(%2)\n" \

      "movaps 80(%2),%%xmm1\n" \
      "movaps 80(%3),%%xmm5\n" \
      "movaps %%xmm1,80(%3)\n" \
      "pslld  $2,%%xmm1\n" \
      "paddd  %%xmm6,%%xmm1\n" \
      "pslld  $1,%%xmm5\n" \
      "psubd  %%xmm5,%%xmm1\n" \
      "movaps %%xmm1,%%xmm5\n" \
      "psrld  $29,%%xmm5\n" \
      "pand   %%xmm7,%%xmm1\n" \
      "paddd  %%xmm5,%%xmm1\n" \
      "paddd  %%xmm5,%%xmm1\n" \
      "paddd  %%xmm5,%%xmm1\n" \
      "movaps %%xmm1,%%xmm5\n" \
      "pcmpgtd  %%xmm4,%%xmm5\n" \
      "pand   48(%4),%%xmm5\n" \
      "psubd  %%xmm5,%%xmm1\n" \
      "movaps %%xmm1,80(%2)\n" \

      "movaps 96(%2),%%xmm2\n" \
      "movaps 96(%3),%%xmm5\n" \
      "movaps %%xmm2,96(%3)\n" \
      "pslld  $2,%%xmm2\n" \
      "paddd  %%xmm6,%%xmm2\n" \
      "pslld  $1,%%xmm5\n" \
      "psubd  %%xmm5,%%xmm2\n" \
      "movaps %%xmm2,%%xmm5\n" \
      "psrld  $29,%%xmm5\n" \
      "pand   %%xmm7,%%xmm2\n" \
      "paddd  %%xmm5,%%xmm2\n" \
      "paddd  %%xmm5,%%xmm2\n" \
      "paddd  %%xmm5,%%xmm2\n" \
      "movaps %%xmm2,%%xmm5\n" \
      "pcmpgtd  %%xmm4,%%xmm5\n" \
      "pand   48(%4),%%xmm5\n" \
      "psubd  %%xmm5,%%xmm2\n" \
      "movaps %%xmm2,96(%2)\n" \

      "movaps 112(%2),%%xmm3\n" \
      "movaps 112(%3),%%xmm5\n" \
      "movaps %%xmm3,112(%3)\n" \
      "pslld  $2,%%xmm3\n" \
      "paddd  %%xmm6,%%xmm3\n" \
      "pslld  $1,%%xmm5\n" \
      "psubd  %%xmm5,%%xmm3\n" \
      "movaps %%xmm3,%%xmm5\n" \
      "psrld  $29,%%xmm5\n" \
      "pand   %%xmm7,%%xmm3\n" \
      "paddd  %%xmm5,%%xmm3\n" \
      "paddd  %%xmm5,%%xmm3\n" \
      "paddd  %%xmm5,%%xmm3\n" \
      "movaps %%xmm3,%%xmm5\n" \
      "pcmpgtd  %%xmm4,%%xmm5\n" \
      "pand   48(%4),%%xmm5\n" \
      "psubd  %%xmm5,%%xmm3\n" \
      "movaps %%xmm3,112(%2)\n" \

      "psrld  $28,%%xmm0\n" \
      "psrld  $28,%%xmm1\n" \
      "psrld  $28,%%xmm2\n" \
      "psrld  $28,%%xmm3\n" \
      "packssdw %%xmm1,%%xmm0\n" \
      "packssdw %%xmm3,%%xmm2\n" \
      "packsswb %%xmm2,%%xmm0\n" \
      "psllw  $7,%%xmm0\n" \
      "pmovmskb %%xmm0,%1\n" \
      "shll $16,%1\n" \
      "addl %1,%0\n" \
      "":"=&r"(output1),"=&r"(output2):"r"(state->xN),"r"(state->xP),"r"(gm29_sse_Consts));
  return output1;
}

__device__ __host__ void gm29_get_sse_state_(gm29_state* state,gm29_sse_state* sse_state){
  int i; for(i=0;i<32;i++) {sse_state->xN[i]=state->xN[i]; sse_state->xP[i]=state->xP[i];}
}

__device__ __host__ unsigned gm29_CNext(unsigned N,unsigned P){
  return (gm29_k*N+gm29_q*(gm29_g-P))%gm29_g;
}

__device__ __host__ unsigned gm29_CNext2(unsigned N,unsigned P,unsigned myk,unsigned myq){
  unsigned long long NNN,PP,kk,qq,gg,rr;                    // returns (myk*N-myq*P) (mod gm29_g)
  NNN=N; PP=P; kk=myk; qq=myq; gg=gm29_g;
  rr=(kk*NNN+qq*(gg-PP));
  NNN=rr>>29;
  PP=rr-(NNN*gg);
  PP-=((PP>>29)*gg);
  return (unsigned)PP;
}

__device__ __host__ unsigned gm29_GetNextN(unsigned x0,unsigned x1,unsigned n){ //returns x_{2^n}
  unsigned myk=gm29_k,myq=gm29_q,i,x=x1;
  for(i=0;i<n;i++){
    x=gm29_CNext2(x,x0,myk,myq);
    myk=gm29_CNext2(myk,2,myk,myq);
    myq=gm29_CNext2(myq,0,myq,0);
  }
  return x;
}

__device__ __host__ unsigned gm29_GetNextAny(unsigned x0,unsigned x1,unsigned long long N){ // returns x_N
  unsigned long long i; unsigned xp=x0,xn=x1,xpnew,xnnew,shift=0;
  i=N; while(i>0){
    if(i%2==1){                        // xp,xn ----> 2^shift
      xpnew=gm29_GetNextN(xp,xn,shift);
      xnnew=gm29_GetNextN(xn,gm29_CNext(xn,xp),shift);
      xp=xpnew; xn=xnnew;
    }
    i/=2; shift++;
  }
  return xp;
}

__device__ __host__ void gm29_skipahead_(gm29_state* state, unsigned long long offset){
  unsigned xn,xp,j; 
  for(j=0;j<32;j++){
    xp=gm29_GetNextAny(state->xP[j],state->xN[j],offset);
    xn=gm29_GetNextAny(state->xP[j],state->xN[j],offset+1);
    state->xP[j]=xp; state->xN[j]=xn;
  }
}

__device__ __host__ void gm29_init_(gm29_state* state){
   unsigned x0=514932,x1=127293,xp,xn,j;
   for(j=0;j<32;j++){
     xp=gm29_GetNextAny(x0,x1,9007198285571818UL);
     xn=gm29_GetNextAny(x0,x1,9007198285571819UL);
     state->xP[j]=xp; state->xN[j]=xn; x0=xp; x1=xn;
   }
}

__device__ __host__ void gm29_init_short_sequence_(gm29_state* state,unsigned SequenceNumber){
  gm29_init_(state);                     // 0 <= SequenceNumber < 10^8;   length of each sequence <= 8*10^7
  gm29_skipahead_(state,82927047ULL*(unsigned long long)SequenceNumber);
}

__device__ __host__ void gm29_init_medium_sequence_(gm29_state* state,unsigned SequenceNumber){
  gm29_init_(state);                     // 0 <= SequenceNumber < 10^6;   length of each sequence <= 8*10^9
  gm29_skipahead_(state,8799201913ULL*(unsigned long long)SequenceNumber);
}

__device__ __host__ void gm29_init_long_sequence_(gm29_state* state,unsigned SequenceNumber){
  gm29_init_(state);                     // 0 <= SequenceNumber < 10^4;   length of each sequence <= 8*10^11
  gm29_skipahead_(state,828317697521ULL*(unsigned long long)SequenceNumber);
}

__device__ __host__ unsigned int gm29_generate_(gm29_state* state){
  unsigned sum=0, i, temp, bit=1;
  for(i=0;i<32;i++){ 
    temp=(gm29_k*state->xN[i]+gm29_q*(gm29_g-state->xP[i]))%gm29_g;
    state->xP[i]=state->xN[i]; state->xN[i]=temp;
    sum+= ((temp<gm29_halfg)?0:bit); bit*=2;
  }
  return sum;
}

__device__ __host__ float gm29_generate_uniform_float_(gm29_state* state){
  unsigned sum=0, i, temp,bit=1;
  for(i=0;i<32;i++){ 
    temp=(gm29_k*state->xN[i]+gm29_q*(gm29_g-state->xP[i]))%gm29_g;
    state->xP[i]=state->xN[i]; state->xN[i]=temp;
    sum+= ((temp<gm29_halfg)?0:bit); bit*=2;
  }
  return ((float) sum) * 2.3283064365386963e-10;
}

__host__ void gm29_print_state_(gm29_state* state){int i;
    printf("Generator State:\nxN={");
    for(i=0;i<32;i++) {printf("%u",state->xN[i]%gm29_g); printf((i<31)?",":"}\nxP={");}
    for(i=0;i<32;i++) {printf("%u",state->xP[i]%gm29_g); printf((i<31)?",":"}\n\n");}
}

__host__ void gm29_print_sse_state_(gm29_sse_state* state){int i;
    printf("Generator State:\nxN={");
    for(i=0;i<32;i++) {printf("%u",state->xN[i]%gm29_g); printf((i<31)?",":"}\nxP={");}
    for(i=0;i<32;i++) {printf("%u",state->xP[i]%gm29_g); printf((i<31)?",":"}\n\n");}
}

__global__ void gm29_kernel_generate_array(gm29_state* state, unsigned int* out, long* length) {
    unsigned temp,sum,i,orbit,seqNum; long offset;

    __shared__ unsigned xP[gm29_THREADS];  // one generator per s=32 threads, i.e. one orbit
    __shared__ unsigned xN[gm29_THREADS];  // per thread, i.e. blockDim.x orbits per block
    __shared__ unsigned  a[gm29_THREADS];  // array "a" contains corresponding parts of output

    orbit   = threadIdx.x % 32;
    seqNum  = (threadIdx.x + blockIdx.x * blockDim.x)>>5;  // RNG_sequence index
    offset  = seqNum*(*length);                            // start of the section in the output array

    xP[threadIdx.x]=gm29_GetNextAny(state->xP[orbit],state->xN[orbit],offset);
    xN[threadIdx.x]=gm29_GetNextAny(state->xP[orbit],state->xN[orbit],offset+1);

    for(i=0;i<(*length);i++){

      temp = gm29_CNext( xN[threadIdx.x], xP[threadIdx.x] );
      xP[threadIdx.x] = xN[threadIdx.x]; xN[threadIdx.x] = temp; 
      a[threadIdx.x]  = (temp < gm29_halfg ? 0 : (1<<orbit) );

      __syncthreads();              // each s=32 threads result in "length" values in the output array
      if((orbit&3)==0) a[threadIdx.x] = a[threadIdx.x]+a[threadIdx.x+1]+a[threadIdx.x+2]+a[threadIdx.x+3];
      __syncthreads();
      if((orbit&15)==0) a[threadIdx.x] = a[threadIdx.x]+a[threadIdx.x+4]+a[threadIdx.x+8]+a[threadIdx.x+12];
      __syncthreads();
      if(orbit==0){ sum=a[threadIdx.x]+a[threadIdx.x+16];  out[offset+i]=sum; }

    }
}

__host__ void gm29_generate_gpu_array_(gm29_state* state, unsigned int* dev_out, long length){

   long          mylength = length/gm29_ARRAY_SECTIONS;
   gm29_state*   dev_state;
   long*         dev_length;

   if((mylength*gm29_ARRAY_SECTIONS)<length) mylength++;

   gm29_CUDA_CALL(cudaMalloc((void**)&dev_state,sizeof(gm29_state)));
   gm29_CUDA_CALL(cudaMalloc((void**)&dev_length,sizeof(long)));
   gm29_CUDA_CALL(cudaMemcpy(dev_state,state,sizeof(gm29_state),cudaMemcpyHostToDevice));
   gm29_CUDA_CALL(cudaMemcpy(dev_length,&mylength,sizeof(long),cudaMemcpyHostToDevice));

   gm29_kernel_generate_array<<<gm29_BLOCKS,gm29_THREADS>>>(dev_state,dev_out,dev_length);
   gm29_CUDA_CALL(cudaGetLastError());
   
   gm29_CUDA_CALL(cudaFree(dev_state)); gm29_CUDA_CALL(cudaFree(dev_length));

}

__global__ void gm29_kernel_generate_array_float(gm29_state* state, float* out, long* length) {
    unsigned temp,sum,i,orbit,seqNum; long offset;

    __shared__ unsigned xP[gm29_THREADS];  // one generator per s=32 threads, i.e. one orbit
    __shared__ unsigned xN[gm29_THREADS];  // per thread, i.e. blockDim.x orbits per block
    __shared__ unsigned  a[gm29_THREADS];  // array "a" contains corresponding parts of output

    orbit   = threadIdx.x % 32;
    seqNum  = (threadIdx.x + blockIdx.x * blockDim.x)>>5;  // RNG_sequence index
    offset  = seqNum*(*length);                            // start of the section in the output array

    xP[threadIdx.x]=gm29_GetNextAny(state->xP[orbit],state->xN[orbit],offset);
    xN[threadIdx.x]=gm29_GetNextAny(state->xP[orbit],state->xN[orbit],offset+1);

    for(i=0;i<(*length);i++){      // each s=32 threads result in "length" values in the output array

      temp = gm29_CNext( xN[threadIdx.x], xP[threadIdx.x] );
      xP[threadIdx.x] = xN[threadIdx.x]; xN[threadIdx.x] = temp; 
      a[threadIdx.x]  = (temp < gm29_halfg ? 0 : (1<<orbit) );

      __syncthreads();
      if((orbit&3)==0) a[threadIdx.x] = a[threadIdx.x]+a[threadIdx.x+1]+a[threadIdx.x+2]+a[threadIdx.x+3];
      __syncthreads();
      if((orbit&15)==0) a[threadIdx.x] = a[threadIdx.x]+a[threadIdx.x+4]+a[threadIdx.x+8]+a[threadIdx.x+12];
      __syncthreads();
      if(orbit==0){ sum=a[threadIdx.x]+a[threadIdx.x+16];  out[offset+i]=((float)sum) * 2.3283064365386963e-10; }

    }
}

__host__ void gm29_generate_gpu_array_float_(gm29_state* state, float* dev_out, long length){

   long          mylength = length/gm29_ARRAY_SECTIONS;
   gm29_state*   dev_state;
   long*         dev_length;

   if((mylength*gm29_ARRAY_SECTIONS)<length) mylength++;

   gm29_CUDA_CALL(cudaMalloc((void**)&dev_state,sizeof(gm29_state)));
   gm29_CUDA_CALL(cudaMalloc((void**)&dev_length,sizeof(long)));
   gm29_CUDA_CALL(cudaMemcpy(dev_state,state,sizeof(gm29_state),cudaMemcpyHostToDevice));
   gm29_CUDA_CALL(cudaMemcpy(dev_length,&mylength,sizeof(long),cudaMemcpyHostToDevice));

   gm29_kernel_generate_array_float<<<gm29_BLOCKS,gm29_THREADS>>>(dev_state,dev_out,dev_length);
   gm29_CUDA_CALL(cudaGetLastError());
   
   gm29_CUDA_CALL(cudaFree(dev_state)); gm29_CUDA_CALL(cudaFree(dev_length));

}

__global__ void gm29_kernel_generate_array_double(gm29_state* state, double* out, long* length) {
    unsigned temp,sum,i,orbit,seqNum; long offset;

    __shared__ unsigned xP[gm29_THREADS];  // one generator per s=32 threads, i.e. one orbit
    __shared__ unsigned xN[gm29_THREADS];  // per thread, i.e. blockDim.x orbits per block
    __shared__ unsigned  a[gm29_THREADS];  // array "a" contains corresponding parts of output

    orbit   = threadIdx.x % 32;
    seqNum  = (threadIdx.x + blockIdx.x * blockDim.x)>>5;  // RNG_sequence index
    offset  = seqNum*(*length);                            // start of the section in the output array

    xP[threadIdx.x]=gm29_GetNextAny(state->xP[orbit],state->xN[orbit],offset);
    xN[threadIdx.x]=gm29_GetNextAny(state->xP[orbit],state->xN[orbit],offset+1);

    for(i=0;i<(*length);i++){       // each s=32 threads result in "length" values in the output array

      temp = gm29_CNext( xN[threadIdx.x], xP[threadIdx.x] );
      xP[threadIdx.x] = xN[threadIdx.x]; xN[threadIdx.x] = temp; 
      a[threadIdx.x]  = (temp < gm29_halfg ? 0 : (1<<orbit) );

      __syncthreads();
      if((orbit&3)==0) a[threadIdx.x] = a[threadIdx.x]+a[threadIdx.x+1]+a[threadIdx.x+2]+a[threadIdx.x+3];
      __syncthreads();
      if((orbit&15)==0) a[threadIdx.x] = a[threadIdx.x]+a[threadIdx.x+4]+a[threadIdx.x+8]+a[threadIdx.x+12];
      __syncthreads();
      if(orbit==0){ sum=a[threadIdx.x]+a[threadIdx.x+16];  out[offset+i]=((double)sum) * 2.3283064365386963e-10; }

    }
}

__host__ void gm29_generate_gpu_array_double_(gm29_state* state, double* dev_out, long length){

   long          mylength = length/gm29_ARRAY_SECTIONS;
   gm29_state*   dev_state;
   long*         dev_length;

   if((mylength*gm29_ARRAY_SECTIONS)<length) mylength++;

   gm29_CUDA_CALL(cudaMalloc((void**)&dev_state,sizeof(gm29_state)));
   gm29_CUDA_CALL(cudaMalloc((void**)&dev_length,sizeof(long)));
   gm29_CUDA_CALL(cudaMemcpy(dev_state,state,sizeof(gm29_state),cudaMemcpyHostToDevice));
   gm29_CUDA_CALL(cudaMemcpy(dev_length,&mylength,sizeof(long),cudaMemcpyHostToDevice));

   gm29_kernel_generate_array_double<<<gm29_BLOCKS,gm29_THREADS>>>(dev_state,dev_out,dev_length);
   gm29_CUDA_CALL(cudaGetLastError());
   
   gm29_CUDA_CALL(cudaFree(dev_state)); gm29_CUDA_CALL(cudaFree(dev_length));

}

__host__ void gm29_generate_array_(gm29_state* state, unsigned int* out, long length){

   long          mylength = length/gm29_ARRAY_SECTIONS;
   gm29_state*   dev_state;
   unsigned int* dev_out;
   long*         dev_length;

   if((mylength*gm29_ARRAY_SECTIONS)<length) mylength++;

   gm29_CUDA_CALL(cudaMalloc((void**)&dev_state,sizeof(gm29_state)));
   gm29_CUDA_CALL(cudaMalloc((void**)&dev_out,mylength*gm29_ARRAY_SECTIONS*sizeof(unsigned int)));
   gm29_CUDA_CALL(cudaMalloc((void**)&dev_length,sizeof(long)));
   gm29_CUDA_CALL(cudaMemcpy(dev_state,state,sizeof(gm29_state),cudaMemcpyHostToDevice));
   gm29_CUDA_CALL(cudaMemcpy(dev_length,&mylength,sizeof(long),cudaMemcpyHostToDevice));

   gm29_kernel_generate_array<<<gm29_BLOCKS,gm29_THREADS>>>(dev_state,dev_out,dev_length);
   gm29_CUDA_CALL(cudaGetLastError());
   
   gm29_CUDA_CALL(cudaMemcpy(out,dev_out,length*sizeof(unsigned int),cudaMemcpyDeviceToHost));
   gm29_CUDA_CALL(cudaFree(dev_state)); gm29_CUDA_CALL(cudaFree(dev_out)); 
   gm29_CUDA_CALL(cudaFree(dev_length));

}
