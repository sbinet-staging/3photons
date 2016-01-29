// (c) Copyright 2013 Lev Barash, Landau Institute for Theoretical Physics, Russian Academy of Sciences
// This is supplement to the paper:
// L.Yu. Barash, L.N. Shchur, "PRAND: GPU accelerated parallel random number generation library: Using most reliable algorithms and applying parallelism of modern GPUs and CPUs".
// e-mail: barash @ itp.ac.ru (remove space)

#include<stdio.h>

#define gm31_CUDA_CALL(x) do { if((x) != cudaSuccess) { printf("Error: %s at %s:%d\n",cudaGetErrorString(cudaGetLastError()),__FILE__,__LINE__); exit(1);}} while(0)

#define gm31_BLOCKS  512
#define gm31_THREADS 128
#define gm31_ARRAY_SECTIONS  (gm31_BLOCKS*gm31_THREADS/32)

#define gm31_qg    30064771058ULL
#define gm31_g     2147483647
#define gm31_halfg 1073741824
#define gm31_k     11
#define gm31_q     14

typedef struct{
  unsigned xN[32],xP[32];
} gm31_state;

typedef struct{
  unsigned xN[64] __attribute__ ((aligned(16))),
           xP[64] __attribute__ ((aligned(16)));
} gm31_sse_state;

unsigned gm31_Consts[16] __attribute__ ((aligned(16))) =
{4294967222UL,36UL,4294967222UL,36UL,gm31_k,0,gm31_k,0,gm31_q,0,gm31_q,0,gm31_g,0,gm31_g,0};

extern "C" __host__ unsigned int gm31_sse_generate_(gm31_sse_state* state){
  unsigned output1,output2;
  asm volatile("\n" \
      "movaps 48(%3),%%xmm5\n" \
      "\n" \
      "movaps (%1),%%xmm0\n" \
      "movaps %%xmm0,%%xmm7\n" \
      "movaps (%2),%%xmm6\n" \
      "pmuludq 16(%3),%%xmm0\n" \
      "pmuludq 32(%3),%%xmm6\n" \
      "paddq (%3),%%xmm0\n" \
      "psubq %%xmm6,%%xmm0\n" \
      "movaps %%xmm0,%%xmm6\n" \
      "psrlq $31,%%xmm6\n" \
      "andps %%xmm5,%%xmm0\n" \
      "paddq %%xmm6,%%xmm0\n" \
      "movaps %%xmm0,(%1)\n" \
      "movaps %%xmm7,(%2)\n" \
      "\n" \
      "movaps 16(%1),%%xmm1\n" \
      "movaps %%xmm1,%%xmm7\n" \
      "movaps 16(%2),%%xmm6\n" \
      "pmuludq 16(%3),%%xmm1\n" \
      "pmuludq 32(%3),%%xmm6\n" \
      "paddq (%3),%%xmm1\n" \
      "psubq %%xmm6,%%xmm1\n" \
      "movaps %%xmm1,%%xmm6\n" \
      "psrlq $31,%%xmm6\n" \
      "andps %%xmm5,%%xmm1\n" \
      "paddq %%xmm6,%%xmm1\n" \
      "movaps %%xmm1,16(%1)\n" \
      "movaps %%xmm7,16(%2)\n" \
      "\n" \
      "movaps 32(%1),%%xmm2\n" \
      "movaps %%xmm2,%%xmm7\n" \
      "movaps 32(%2),%%xmm6\n" \
      "pmuludq 16(%3),%%xmm2\n" \
      "pmuludq 32(%3),%%xmm6\n" \
      "paddq (%3),%%xmm2\n" \
      "psubq %%xmm6,%%xmm2\n" \
      "movaps %%xmm2,%%xmm6\n" \
      "psrlq $31,%%xmm6\n" \
      "andps %%xmm5,%%xmm2\n" \
      "paddq %%xmm6,%%xmm2\n" \
      "movaps %%xmm2,32(%1)\n" \
      "movaps %%xmm7,32(%2)\n" \
      "\n" \
      "movaps 48(%1),%%xmm3\n" \
      "movaps %%xmm3,%%xmm7\n" \
      "movaps 48(%2),%%xmm6\n" \
      "pmuludq 16(%3),%%xmm3\n" \
      "pmuludq 32(%3),%%xmm6\n" \
      "paddq (%3),%%xmm3\n" \
      "psubq %%xmm6,%%xmm3\n" \
      "movaps %%xmm3,%%xmm6\n" \
      "psrlq $31,%%xmm6\n" \
      "andps %%xmm5,%%xmm3\n" \
      "paddq %%xmm6,%%xmm3\n" \
      "movaps %%xmm3,48(%1)\n" \
      "movaps %%xmm7,48(%2)\n" \
      "\n" \
      "shufps $136,%%xmm1,%%xmm0\n" \
      "shufps $136,%%xmm3,%%xmm2\n" \
      "psrld  $30,%%xmm0\n" \
      "psrld  $30,%%xmm2\n" \
      "packssdw %%xmm2,%%xmm0\n" \
      "\n" \
      "movaps 64(%1),%%xmm1\n" \
      "movaps %%xmm1,%%xmm7\n" \
      "movaps 64(%2),%%xmm6\n" \
      "pmuludq 16(%3),%%xmm1\n" \
      "pmuludq 32(%3),%%xmm6\n" \
      "paddq (%3),%%xmm1\n" \
      "psubq %%xmm6,%%xmm1\n" \
      "movaps %%xmm1,%%xmm6\n" \
      "psrlq $31,%%xmm6\n" \
      "andps %%xmm5,%%xmm1\n" \
      "paddq %%xmm6,%%xmm1\n" \
      "movaps %%xmm1,64(%1)\n" \
      "movaps %%xmm7,64(%2)\n" \
      "\n" \
      "movaps 80(%1),%%xmm2\n" \
      "movaps %%xmm2,%%xmm7\n" \
      "movaps 80(%2),%%xmm6\n" \
      "pmuludq 16(%3),%%xmm2\n" \
      "pmuludq 32(%3),%%xmm6\n" \
      "paddq (%3),%%xmm2\n" \
      "psubq %%xmm6,%%xmm2\n" \
      "movaps %%xmm2,%%xmm6\n" \
      "psrlq $31,%%xmm6\n" \
      "andps %%xmm5,%%xmm2\n" \
      "paddq %%xmm6,%%xmm2\n" \
      "movaps %%xmm2,80(%1)\n" \
      "movaps %%xmm7,80(%2)\n" \
      "\n" \
      "movaps 96(%1),%%xmm3\n" \
      "movaps %%xmm3,%%xmm7\n" \
      "movaps 96(%2),%%xmm6\n" \
      "pmuludq 16(%3),%%xmm3\n" \
      "pmuludq 32(%3),%%xmm6\n" \
      "paddq (%3),%%xmm3\n" \
      "psubq %%xmm6,%%xmm3\n" \
      "movaps %%xmm3,%%xmm6\n" \
      "psrlq $31,%%xmm6\n" \
      "andps %%xmm5,%%xmm3\n" \
      "paddq %%xmm6,%%xmm3\n" \
      "movaps %%xmm3,96(%1)\n" \
      "movaps %%xmm7,96(%2)\n" \
      "\n" \
      "movaps 112(%1),%%xmm4\n" \
      "movaps %%xmm4,%%xmm7\n" \
      "movaps 112(%2),%%xmm6\n" \
      "pmuludq 16(%3),%%xmm4\n" \
      "pmuludq 32(%3),%%xmm6\n" \
      "paddq (%3),%%xmm4\n" \
      "psubq %%xmm6,%%xmm4\n" \
      "movaps %%xmm4,%%xmm6\n" \
      "psrlq $31,%%xmm6\n" \
      "andps %%xmm5,%%xmm4\n" \
      "paddq %%xmm6,%%xmm4\n" \
      "movaps %%xmm4,112(%1)\n" \
      "movaps %%xmm7,112(%2)\n" \
      "\n" \
      "shufps $136,%%xmm2,%%xmm1\n" \
      "shufps $136,%%xmm4,%%xmm3\n" \
      "psrld  $30,%%xmm1\n" \
      "psrld  $30,%%xmm3\n" \
      "packssdw %%xmm3,%%xmm1\n" \
      "\n" \
      "packsswb %%xmm1,%%xmm0\n" \
      "psllw $7,%%xmm0\n" \
      "pmovmskb %%xmm0,%0\n" \
      "":"=r"(output1):"r"(state->xN),"r"(state->xP),"r"(gm31_Consts));
  asm volatile("\n" \
      "movaps 128(%1),%%xmm0\n" \
      "movaps %%xmm0,%%xmm7\n" \
      "movaps 128(%2),%%xmm6\n" \
      "pmuludq 16(%3),%%xmm0\n" \
      "pmuludq 32(%3),%%xmm6\n" \
      "paddq (%3),%%xmm0\n" \
      "psubq %%xmm6,%%xmm0\n" \
      "movaps %%xmm0,%%xmm6\n" \
      "psrlq $31,%%xmm6\n" \
      "andps %%xmm5,%%xmm0\n" \
      "paddq %%xmm6,%%xmm0\n" \
      "movaps %%xmm0,128(%1)\n" \
      "movaps %%xmm7,128(%2)\n" \
      "\n" \
      "movaps 144(%1),%%xmm1\n" \
      "movaps %%xmm1,%%xmm7\n" \
      "movaps 144(%2),%%xmm6\n" \
      "pmuludq 16(%3),%%xmm1\n" \
      "pmuludq 32(%3),%%xmm6\n" \
      "paddq (%3),%%xmm1\n" \
      "psubq %%xmm6,%%xmm1\n" \
      "movaps %%xmm1,%%xmm6\n" \
      "psrlq $31,%%xmm6\n" \
      "andps %%xmm5,%%xmm1\n" \
      "paddq %%xmm6,%%xmm1\n" \
      "movaps %%xmm1,144(%1)\n" \
      "movaps %%xmm7,144(%2)\n" \
      "\n" \
      "movaps 160(%1),%%xmm2\n" \
      "movaps %%xmm2,%%xmm7\n" \
      "movaps 160(%2),%%xmm6\n" \
      "pmuludq 16(%3),%%xmm2\n" \
      "pmuludq 32(%3),%%xmm6\n" \
      "paddq (%3),%%xmm2\n" \
      "psubq %%xmm6,%%xmm2\n" \
      "movaps %%xmm2,%%xmm6\n" \
      "psrlq $31,%%xmm6\n" \
      "andps %%xmm5,%%xmm2\n" \
      "paddq %%xmm6,%%xmm2\n" \
      "movaps %%xmm2,160(%1)\n" \
      "movaps %%xmm7,160(%2)\n" \
      "\n" \
      "movaps 176(%1),%%xmm3\n" \
      "movaps %%xmm3,%%xmm7\n" \
      "movaps 176(%2),%%xmm6\n" \
      "pmuludq 16(%3),%%xmm3\n" \
      "pmuludq 32(%3),%%xmm6\n" \
      "paddq (%3),%%xmm3\n" \
      "psubq %%xmm6,%%xmm3\n" \
      "movaps %%xmm3,%%xmm6\n" \
      "psrlq $31,%%xmm6\n" \
      "andps %%xmm5,%%xmm3\n" \
      "paddq %%xmm6,%%xmm3\n" \
      "movaps %%xmm3,176(%1)\n" \
      "movaps %%xmm7,176(%2)\n" \
      "\n" \
      "shufps $136,%%xmm1,%%xmm0\n" \
      "shufps $136,%%xmm3,%%xmm2\n" \
      "psrld  $30,%%xmm0\n" \
      "psrld  $30,%%xmm2\n" \
      "packssdw %%xmm2,%%xmm0\n" \
      "\n" \
      "movaps 192(%1),%%xmm1\n" \
      "movaps %%xmm1,%%xmm7\n" \
      "movaps 192(%2),%%xmm6\n" \
      "pmuludq 16(%3),%%xmm1\n" \
      "pmuludq 32(%3),%%xmm6\n" \
      "paddq (%3),%%xmm1\n" \
      "psubq %%xmm6,%%xmm1\n" \
      "movaps %%xmm1,%%xmm6\n" \
      "psrlq $31,%%xmm6\n" \
      "andps %%xmm5,%%xmm1\n" \
      "paddq %%xmm6,%%xmm1\n" \
      "movaps %%xmm1,192(%1)\n" \
      "movaps %%xmm7,192(%2)\n" \
      "\n" \
      "movaps 208(%1),%%xmm2\n" \
      "movaps %%xmm2,%%xmm7\n" \
      "movaps 208(%2),%%xmm6\n" \
      "pmuludq 16(%3),%%xmm2\n" \
      "pmuludq 32(%3),%%xmm6\n" \
      "paddq (%3),%%xmm2\n" \
      "psubq %%xmm6,%%xmm2\n" \
      "movaps %%xmm2,%%xmm6\n" \
      "psrlq $31,%%xmm6\n" \
      "andps %%xmm5,%%xmm2\n" \
      "paddq %%xmm6,%%xmm2\n" \
      "movaps %%xmm2,208(%1)\n" \
      "movaps %%xmm7,208(%2)\n" \
      "\n" \
      "movaps 224(%1),%%xmm3\n" \
      "movaps %%xmm3,%%xmm7\n" \
      "movaps 224(%2),%%xmm6\n" \
      "pmuludq 16(%3),%%xmm3\n" \
      "pmuludq 32(%3),%%xmm6\n" \
      "paddq (%3),%%xmm3\n" \
      "psubq %%xmm6,%%xmm3\n" \
      "movaps %%xmm3,%%xmm6\n" \
      "psrlq $31,%%xmm6\n" \
      "andps %%xmm5,%%xmm3\n" \
      "paddq %%xmm6,%%xmm3\n" \
      "movaps %%xmm3,224(%1)\n" \
      "movaps %%xmm7,224(%2)\n" \
      "\n" \
      "movaps 240(%1),%%xmm4\n" \
      "movaps %%xmm4,%%xmm7\n" \
      "movaps 240(%2),%%xmm6\n" \
      "pmuludq 16(%3),%%xmm4\n" \
      "pmuludq 32(%3),%%xmm6\n" \
      "paddq (%3),%%xmm4\n" \
      "psubq %%xmm6,%%xmm4\n" \
      "movaps %%xmm4,%%xmm6\n" \
      "psrlq $31,%%xmm6\n" \
      "andps %%xmm5,%%xmm4\n" \
      "paddq %%xmm6,%%xmm4\n" \
      "movaps %%xmm4,240(%1)\n" \
      "movaps %%xmm7,240(%2)\n" \
      "\n" \
      "shufps $136,%%xmm2,%%xmm1\n" \
      "shufps $136,%%xmm4,%%xmm3\n" \
      "psrld  $30,%%xmm1\n" \
      "psrld  $30,%%xmm3\n" \
      "packssdw %%xmm3,%%xmm1\n" \
      "\n" \
      "packsswb %%xmm1,%%xmm0\n" \
      "psllw $7,%%xmm0\n" \
      "pmovmskb %%xmm0,%0\n" \
      "shll $16,%0\n" \
      "":"=r"(output2):"r"(state->xN),"r"(state->xP),"r"(gm31_Consts));
  asm volatile("\n" \
      "addl %1,%0\n" \
      "\n" \
      "":"=r"(output2):"r"(output1),"0"(output2));
  return output2;
}

extern "C" __device__ __host__ void gm31_get_sse_state_(gm31_state* state,gm31_sse_state* sse_state){
  int i;
  for(i=0;i<32;i++) {
    sse_state->xN[2*i]=state->xN[i]; sse_state->xP[2*i]=state->xP[i];
    sse_state->xN[2*i+1]=sse_state->xP[2*i+1]=0;
  }
}

extern "C" __device__ __host__ unsigned gm31_mod_g(unsigned long long x){ // returns x (mod g)
  unsigned long long F,G; G=x;
  do{ F=(G>>31); G = (G-(F<<31)+F); } while(G>gm31_g);
  return G;
}

extern "C" __device__ __host__ unsigned gm31_CNext(unsigned N,unsigned P){
  unsigned long long curr1,curr2,curr3;
  curr1=(unsigned long long)gm31_k*(unsigned long long)N; 
  curr2=(unsigned long long)gm31_q*(unsigned long long)P; 
  curr3=gm31_mod_g(gm31_qg+curr1-curr2); return curr3;
}

extern "C" __device__ __host__ unsigned gm31_CNext2(unsigned N,unsigned P,unsigned myk,unsigned myq){
  unsigned long long curr1,curr2,curr3;
  curr1=(unsigned long long)myk*(unsigned long long)N;
  curr2=(unsigned long long)myq*(unsigned long long)P;
  curr3=gm31_mod_g((unsigned long long)myq*(unsigned long long)gm31_g+curr1-curr2);
  return curr3;
}

extern "C" __device__ __host__ unsigned gm31_GetNextN(unsigned x0,unsigned x1,unsigned n){ // returns x_{2^n}
  unsigned myk=gm31_k,myq=gm31_q,i,x=x1;
  for(i=0;i<n;i++){
    x=gm31_CNext2(x,x0,myk,myq);
    myk=gm31_CNext2(myk,2,myk,myq);
    myq=gm31_CNext2(myq,0,myq,0);
  }
  return x;
}

extern "C" __device__ __host__ unsigned gm31_GetNextAny(unsigned x0,unsigned x1,unsigned long long N){ // returns x_N
  unsigned long long i; unsigned xp=x0,xn=x1,xpnew,xnnew,shift=0;
  i=N; while(i>0){
    if(i%2==1){                        // xp,xn ----> 2^shift
      xpnew=gm31_GetNextN(xp,xn,shift);
      xnnew=gm31_GetNextN(xn,gm31_CNext(xn,xp),shift);
      xp=xpnew; xn=xnnew;
    }
    i/=2; shift++;
  }
  return xp;
}

extern "C" __device__ __host__ void gm31_skipahead_(gm31_state* state, unsigned long long offset){
  unsigned xn,xp,j; 
  for(j=0;j<32;j++){
    xp=gm31_GetNextAny(state->xP[j],state->xN[j],offset);
    xn=gm31_GetNextAny(state->xP[j],state->xN[j],offset+1);
    state->xP[j]=xp; state->xN[j]=xn;
  }
}

extern "C" __device__ __host__ void gm31_init_(gm31_state* state){
  unsigned x0=554937932UL,x1=1253942293UL,xp,xn,j;
  for(j=0;j<32;j++){
    xp=gm31_GetNextAny(x0,x1,144115183781032008ULL);
    xn=gm31_GetNextAny(x0,x1,144115183781032009ULL);
    state->xP[j]=xp; state->xN[j]=xn; x0=xp; x1=xn;
  }
}

extern "C" __device__ __host__ void gm31_init_short_sequence_(gm31_state* state,unsigned SequenceNumber){
  gm31_init_(state);                     // 0 <= SequenceNumber < 10^9;   length of each sequence <= 8*10^7
  gm31_skipahead_(state,82927047ULL*(unsigned long long)SequenceNumber);
}

extern "C" __device__ __host__ void gm31_init_medium_sequence_(gm31_state* state,unsigned SequenceNumber){
  gm31_init_(state);                     // 0 <= SequenceNumber < 10^7;   length of each sequence <= 8*10^9
  gm31_skipahead_(state,8799201913ULL*(unsigned long long)SequenceNumber);
}

extern "C" __device__ __host__ void gm31_init_long_sequence_(gm31_state* state,unsigned SequenceNumber){
  gm31_init_(state);                     // 0 <= SequenceNumber < 10^5;   length of each sequence <= 8*10^11
  gm31_skipahead_(state,828317697521ULL*(unsigned long long)SequenceNumber);
}

extern "C" __device__ __host__ unsigned int gm31_generate_(gm31_state* state){
  int i; unsigned temp,sum=0,bit=1;
  for(i=0;i<32;i++){
    temp=gm31_CNext(state->xN[i],state->xP[i]);
    state->xP[i]=state->xN[i]; state->xN[i]=temp;
    sum += ((temp<gm31_halfg)?0:bit); bit*=2;
  }
  return sum;
}

extern "C" __device__ __host__ float gm31_generate_uniform_float_(gm31_state* state){
  int i; unsigned temp,sum=0,bit=1;
  for(i=0;i<32;i++){
    temp=gm31_CNext(state->xN[i],state->xP[i]);
    state->xP[i]=state->xN[i]; state->xN[i]=temp;
    sum += ((temp<gm31_halfg)?0:bit); bit*=2;
  }
  return ((float) sum) * 2.3283064365386963e-10;
}

extern "C" __host__ void gm31_print_state_(gm31_state* state){int i;
    printf("Generator State:\nxN={");
    for(i=0;i<32;i++) {printf("%u",state->xN[i]%gm31_g); printf((i<31)?",":"}\nxP={");}
    for(i=0;i<32;i++) {printf("%u",state->xP[i]%gm31_g); printf((i<31)?",":"}\n\n");}
}

extern "C" __host__ void gm31_print_sse_state_(gm31_sse_state* state){int i;
    printf("Generator State:\nxN={");
    for(i=0;i<32;i++) {printf("%u",state->xN[2*i]%gm31_g); printf((i<31)?",":"}\nxP={");}
    for(i=0;i<32;i++) {printf("%u",state->xP[2*i]%gm31_g); printf((i<31)?",":"}\n\n");}
}

__global__ void gm31_kernel_generate_array(gm31_state* state, unsigned int* out, long* length) {
    unsigned temp,i,orbit,seqNum; long offset;

    __shared__ unsigned xP[gm31_THREADS];  // one generator per s=32 threads, i.e. one orbit
    __shared__ unsigned xN[gm31_THREADS];  // per thread, i.e. blockDim.x orbits per block
    __shared__ unsigned  a[gm31_THREADS];  // array "a" contains corresponding parts of output

    orbit   = threadIdx.x % 32;
    seqNum  = (threadIdx.x + blockIdx.x * blockDim.x)>>5;  // RNG_sequence index
    offset  = seqNum*(*length);                            // start of the section in the output array

    xP[threadIdx.x]=gm31_GetNextAny(state->xP[orbit],state->xN[orbit],offset);
    xN[threadIdx.x]=gm31_GetNextAny(state->xP[orbit],state->xN[orbit],offset+1);

    for(i=0;i<(*length);i++){       // each s=32 threads result in "length" values in the output array

      temp = gm31_CNext( xN[threadIdx.x], xP[threadIdx.x] );
      xP[threadIdx.x] = xN[threadIdx.x]; xN[threadIdx.x] = temp; 
      a[threadIdx.x]  = (temp < gm31_halfg ? 0 : (1<<orbit) );

      __syncthreads();
      if((orbit&3)==0) a[threadIdx.x] = a[threadIdx.x]+a[threadIdx.x+1]+a[threadIdx.x+2]+a[threadIdx.x+3];
      __syncthreads();
      if((orbit&15)==0) a[threadIdx.x] = a[threadIdx.x]+a[threadIdx.x+4]+a[threadIdx.x+8]+a[threadIdx.x+12];
      __syncthreads();
      if(orbit==0){ out[offset+i]=a[threadIdx.x]+a[threadIdx.x+16]; }

    }
}

extern "C" __host__ void gm31_generate_gpu_array_(gm31_state* state, unsigned int* dev_out, unsigned int* length){

   long          mylength = (*length)/gm31_ARRAY_SECTIONS;
   gm31_state*   dev_state;
   long*         dev_length;

   if((mylength*gm31_ARRAY_SECTIONS)<(*length)) mylength++;

   gm31_CUDA_CALL(cudaMalloc((void**)&dev_state,sizeof(gm31_state)));
   gm31_CUDA_CALL(cudaMalloc((void**)&dev_length,sizeof(long)));
   gm31_CUDA_CALL(cudaMemcpy(dev_state,state,sizeof(gm31_state),cudaMemcpyHostToDevice));
   gm31_CUDA_CALL(cudaMemcpy(dev_length,&mylength,sizeof(long),cudaMemcpyHostToDevice));

   gm31_kernel_generate_array<<<gm31_BLOCKS,gm31_THREADS>>>(dev_state,dev_out,dev_length);
   gm31_CUDA_CALL(cudaGetLastError());
   
   gm31_CUDA_CALL(cudaFree(dev_state)); gm31_CUDA_CALL(cudaFree(dev_length));

}

__global__ void gm31_kernel_generate_array_float(gm31_state* state, float* out, long* length) {
    unsigned temp,sum,i,orbit,seqNum; long offset;

    __shared__ unsigned xP[gm31_THREADS];  // one generator per s=32 threads, i.e. one orbit
    __shared__ unsigned xN[gm31_THREADS];  // per thread, i.e. blockDim.x orbits per block
    __shared__ unsigned  a[gm31_THREADS];  // array "a" contains corresponding parts of output

    orbit   = threadIdx.x % 32;
    seqNum  = (threadIdx.x + blockIdx.x * blockDim.x)>>5;  // RNG_sequence index
    offset  = seqNum*(*length);                            // start of the section in the output array

    xP[threadIdx.x]=gm31_GetNextAny(state->xP[orbit],state->xN[orbit],offset);
    xN[threadIdx.x]=gm31_GetNextAny(state->xP[orbit],state->xN[orbit],offset+1);

    for(i=0;i<(*length);i++){       // each s=32 threads result in "length" values in the output array

      temp = gm31_CNext( xN[threadIdx.x], xP[threadIdx.x] );
      xP[threadIdx.x] = xN[threadIdx.x]; xN[threadIdx.x] = temp; 
      a[threadIdx.x]  = (temp < gm31_halfg ? 0 : (1<<orbit) );

      __syncthreads();
      if((orbit&3)==0) a[threadIdx.x] = a[threadIdx.x]+a[threadIdx.x+1]+a[threadIdx.x+2]+a[threadIdx.x+3];
      __syncthreads();
      if((orbit&15)==0) a[threadIdx.x] = a[threadIdx.x]+a[threadIdx.x+4]+a[threadIdx.x+8]+a[threadIdx.x+12];
      __syncthreads();
      if(orbit==0){ sum=a[threadIdx.x]+a[threadIdx.x+16];  out[offset+i]=((float)sum) * 2.3283064365386963e-10; }

    }
}

extern "C" __host__ void gm31_generate_gpu_array_float_(gm31_state* state, float* dev_out, unsigned int* length){

   long          mylength = (*length)/gm31_ARRAY_SECTIONS;
   gm31_state*   dev_state;
   long*         dev_length;

   if((mylength*gm31_ARRAY_SECTIONS)<(*length)) mylength++;

   gm31_CUDA_CALL(cudaMalloc((void**)&dev_state,sizeof(gm31_state)));
   gm31_CUDA_CALL(cudaMalloc((void**)&dev_length,sizeof(long)));
   gm31_CUDA_CALL(cudaMemcpy(dev_state,state,sizeof(gm31_state),cudaMemcpyHostToDevice));
   gm31_CUDA_CALL(cudaMemcpy(dev_length,&mylength,sizeof(long),cudaMemcpyHostToDevice));

   gm31_kernel_generate_array_float<<<gm31_BLOCKS,gm31_THREADS>>>(dev_state,dev_out,dev_length);
   gm31_CUDA_CALL(cudaGetLastError());
   
   gm31_CUDA_CALL(cudaFree(dev_state)); gm31_CUDA_CALL(cudaFree(dev_length));

}

__global__ void gm31_kernel_generate_array_double(gm31_state* state, double* out, long* length) {
    unsigned temp,sum,i,orbit,seqNum; long offset;

    __shared__ unsigned xP[gm31_THREADS];  // one generator per s=32 threads, i.e. one orbit
    __shared__ unsigned xN[gm31_THREADS];  // per thread, i.e. blockDim.x orbits per block
    __shared__ unsigned  a[gm31_THREADS];  // array "a" contains corresponding parts of output

    orbit   = threadIdx.x % 32;
    seqNum  = (threadIdx.x + blockIdx.x * blockDim.x)>>5;  // RNG_sequence index
    offset  = seqNum*(*length);                            // start of the section in the output array

    xP[threadIdx.x]=gm31_GetNextAny(state->xP[orbit],state->xN[orbit],offset);
    xN[threadIdx.x]=gm31_GetNextAny(state->xP[orbit],state->xN[orbit],offset+1);

    for(i=0;i<(*length);i++){       // each s=32 threads result in "length" values in the output array

      temp = gm31_CNext( xN[threadIdx.x], xP[threadIdx.x] );
      xP[threadIdx.x] = xN[threadIdx.x]; xN[threadIdx.x] = temp; 
      a[threadIdx.x]  = (temp < gm31_halfg ? 0 : (1<<orbit) );

      __syncthreads();
      if((orbit&3)==0) a[threadIdx.x] = a[threadIdx.x]+a[threadIdx.x+1]+a[threadIdx.x+2]+a[threadIdx.x+3];
      __syncthreads();
      if((orbit&15)==0) a[threadIdx.x] = a[threadIdx.x]+a[threadIdx.x+4]+a[threadIdx.x+8]+a[threadIdx.x+12];
      __syncthreads();
      if(orbit==0){ sum=a[threadIdx.x]+a[threadIdx.x+16];  out[offset+i]=((double)sum) * 2.3283064365386963e-10; }

    }
}

extern "C" __host__ void gm31_generate_gpu_array_double_(gm31_state* state, double* dev_out, unsigned int* length){

   long          mylength = (*length)/gm31_ARRAY_SECTIONS;
   gm31_state*   dev_state;
   long*         dev_length;

   if((mylength*gm31_ARRAY_SECTIONS)<(*length)) mylength++;

   gm31_CUDA_CALL(cudaMalloc((void**)&dev_state,sizeof(gm31_state)));
   gm31_CUDA_CALL(cudaMalloc((void**)&dev_length,sizeof(long)));
   gm31_CUDA_CALL(cudaMemcpy(dev_state,state,sizeof(gm31_state),cudaMemcpyHostToDevice));
   gm31_CUDA_CALL(cudaMemcpy(dev_length,&mylength,sizeof(long),cudaMemcpyHostToDevice));

   gm31_kernel_generate_array_double<<<gm31_BLOCKS,gm31_THREADS>>>(dev_state,dev_out,dev_length);
   gm31_CUDA_CALL(cudaGetLastError());
   
   gm31_CUDA_CALL(cudaFree(dev_state)); gm31_CUDA_CALL(cudaFree(dev_length));

}

extern "C" __host__ void gm31_generate_array_(gm31_state* state, unsigned int* out, unsigned int* length){

   long          mylength = (*length)/gm31_ARRAY_SECTIONS;
   gm31_state*   dev_state;
   unsigned int* dev_out;
   long*         dev_length;

   if((mylength*gm31_ARRAY_SECTIONS)<(*length)) mylength++;

   gm31_CUDA_CALL(cudaMalloc((void**)&dev_state,sizeof(gm31_state)));
   gm31_CUDA_CALL(cudaMalloc((void**)&dev_out,mylength*gm31_ARRAY_SECTIONS*sizeof(unsigned int)));
   gm31_CUDA_CALL(cudaMalloc((void**)&dev_length,sizeof(long)));
   gm31_CUDA_CALL(cudaMemcpy(dev_state,state,sizeof(gm31_state),cudaMemcpyHostToDevice));
   gm31_CUDA_CALL(cudaMemcpy(dev_length,&mylength,sizeof(long),cudaMemcpyHostToDevice));

   gm31_kernel_generate_array<<<gm31_BLOCKS,gm31_THREADS>>>(dev_state,dev_out,dev_length);
   gm31_CUDA_CALL(cudaGetLastError());
   
   gm31_CUDA_CALL(cudaMemcpy(out,dev_out,(*length)*sizeof(unsigned int),cudaMemcpyDeviceToHost));
   gm31_CUDA_CALL(cudaFree(dev_state)); gm31_CUDA_CALL(cudaFree(dev_out)); 
   gm31_CUDA_CALL(cudaFree(dev_length));

}
