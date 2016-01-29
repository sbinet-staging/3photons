// (c) Copyright 2013 Lev Barash, Landau Institute for Theoretical Physics, Russian Academy of Sciences
// This is supplement to the paper:
// L.Yu. Barash, L.N. Shchur, "PRAND: GPU accelerated parallel random number generation library: Using most reliable algorithms and applying parallelism of modern GPUs and CPUs".
// e-mail: barash @ itp.ac.ru (remove space)

#include<stdio.h>

#define gq58x4_CUDA_CALL(x) do { if((x) != cudaSuccess) { printf("Error: %s at %s:%d\n",cudaGetErrorString(cudaGetLastError()),__FILE__,__LINE__); exit(1);}} while(0)

#define gq58x4_BLOCKS  128
#define gq58x4_THREADS 128
#define gq58x4_ARRAY_SECTIONS (gq58x4_BLOCKS*gq58x4_THREADS/8)

#define gq58x4_k 8
#define gq58x4_q 48
#define gq58x4_g       288230374541099008ULL
#define gq58x4_gdiv16  18014398408818688ULL

typedef unsigned long long lt;

typedef struct{
  lt xN[8] __attribute__ ((aligned(16))),
     xP[8] __attribute__ ((aligned(16)));
} gq58x4_state;

typedef gq58x4_state gq58x4_sse_state;

lt gq58x4_sse_Consts[8] __attribute__ ((aligned(16))) =
   {13835057977972752384ULL,13835057977972752384ULL,1610612736ULL,1610612736ULL,
    288230371923853311ULL,288230371923853311ULL,288230374541099008ULL,288230374541099008ULL};

__host__ unsigned int gq58x4_sse_generate_(gq58x4_sse_state* state){
  unsigned output;
  asm volatile("movaps (%3),%%xmm0\n" \

      "movaps (%2),%%xmm1\n" \
      "movaps (%1),%%xmm4\n" \
      "movaps %%xmm4,(%2)\n" \
      "psllq  $3,%%xmm4\n" \
      "paddq  %%xmm0,%%xmm4\n" \
      "psllq  $4,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm4\n" \
      "psllq  $1,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm4\n" \
      "movaps %%xmm4,%%xmm1\n" \
      "psrlq  $58,%%xmm1\n" \
      "psllq  $29,%%xmm1\n" \
      "movaps %%xmm1,%%xmm3\n" \
      "psllq  $1,%%xmm3\n" \
      "paddq  %%xmm1,%%xmm3\n" \
      "psllq  $29,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm4\n" \
      "paddq  %%xmm3,%%xmm4\n" \
      "movaps %%xmm4,%%xmm1\n" \
      "paddq  16(%3),%%xmm1\n" \
      "pshufd $245,%%xmm1,%%xmm3\n" \
      "pcmpgtd 32(%3),%%xmm3\n" \
      "pand    48(%3),%%xmm3\n" \
      "psubq   %%xmm3,%%xmm4\n" \
      "movaps %%xmm4,(%1)\n" \
      "movaps %%xmm4,%%xmm1\n" \
      "paddq  %%xmm4,%%xmm1\n" \
      "paddq  %%xmm4,%%xmm1\n" \
      "psrlq  $29,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm4\n" \


      "movaps 16(%2),%%xmm1\n" \
      "movaps 16(%1),%%xmm5\n" \
      "movaps %%xmm5,16(%2)\n" \
      "psllq  $3,%%xmm5\n" \
      "paddq  %%xmm0,%%xmm5\n" \
      "psllq  $4,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm5\n" \
      "psllq  $1,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm5\n" \
      "movaps %%xmm5,%%xmm1\n" \
      "psrlq  $58,%%xmm1\n" \
      "psllq  $29,%%xmm1\n" \
      "movaps %%xmm1,%%xmm3\n" \
      "psllq  $1,%%xmm3\n" \
      "paddq  %%xmm1,%%xmm3\n" \
      "psllq  $29,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm5\n" \
      "paddq  %%xmm3,%%xmm5\n" \
      "movaps %%xmm5,%%xmm1\n" \
      "paddq  16(%3),%%xmm1\n" \
      "pshufd $245,%%xmm1,%%xmm3\n" \
      "pcmpgtd 32(%3),%%xmm3\n" \
      "pand    48(%3),%%xmm3\n" \
      "psubq   %%xmm3,%%xmm5\n" \
      "movaps %%xmm5,16(%1)\n" \
      "movaps %%xmm5,%%xmm1\n" \
      "paddq  %%xmm5,%%xmm1\n" \
      "paddq  %%xmm5,%%xmm1\n" \
      "psrlq  $29,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm5\n" \

      "movaps 32(%2),%%xmm1\n" \
      "movaps 32(%1),%%xmm6\n" \
      "movaps %%xmm6,32(%2)\n" \
      "psllq  $3,%%xmm6\n" \
      "paddq  %%xmm0,%%xmm6\n" \
      "psllq  $4,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm6\n" \
      "psllq  $1,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm6\n" \
      "movaps %%xmm6,%%xmm1\n" \
      "psrlq  $58,%%xmm1\n" \
      "psllq  $29,%%xmm1\n" \
      "movaps %%xmm1,%%xmm3\n" \
      "psllq  $1,%%xmm3\n" \
      "paddq  %%xmm1,%%xmm3\n" \
      "psllq  $29,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm6\n" \
      "paddq  %%xmm3,%%xmm6\n" \
      "movaps %%xmm6,%%xmm1\n" \
      "paddq  16(%3),%%xmm1\n" \
      "pshufd $245,%%xmm1,%%xmm3\n" \
      "pcmpgtd 32(%3),%%xmm3\n" \
      "pand    48(%3),%%xmm3\n" \
      "psubq   %%xmm3,%%xmm6\n" \
      "movaps %%xmm6,32(%1)\n" \
      "movaps %%xmm6,%%xmm1\n" \
      "paddq  %%xmm6,%%xmm1\n" \
      "paddq  %%xmm6,%%xmm1\n" \
      "psrlq  $29,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm6\n" \

      "movaps 48(%2),%%xmm1\n" \
      "movaps 48(%1),%%xmm7\n" \
      "movaps %%xmm7,48(%2)\n" \
      "psllq  $3,%%xmm7\n" \
      "paddq  %%xmm0,%%xmm7\n" \
      "psllq  $4,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm7\n" \
      "psllq  $1,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm7\n" \
      "movaps %%xmm7,%%xmm1\n" \
      "psrlq  $58,%%xmm1\n" \
      "psllq  $29,%%xmm1\n" \
      "movaps %%xmm1,%%xmm3\n" \
      "psllq  $1,%%xmm3\n" \
      "paddq  %%xmm1,%%xmm3\n" \
      "psllq  $29,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm7\n" \
      "paddq  %%xmm3,%%xmm7\n" \
      "movaps %%xmm7,%%xmm1\n" \
      "paddq  16(%3),%%xmm1\n" \
      "pshufd $245,%%xmm1,%%xmm3\n" \
      "pcmpgtd 32(%3),%%xmm3\n" \
      "pand    48(%3),%%xmm3\n" \
      "psubq   %%xmm3,%%xmm7\n" \
      "movaps %%xmm7,48(%1)\n" \
      "movaps %%xmm7,%%xmm1\n" \
      "paddq  %%xmm7,%%xmm1\n" \
      "paddq  %%xmm7,%%xmm1\n" \
      "psrlq  $29,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm7\n" \

      "psrlq  $54,%%xmm4\n" \
      "psrlq  $54,%%xmm5\n" \
      "psrlq  $54,%%xmm6\n" \
      "psrlq  $54,%%xmm7\n" \
      "packssdw  %%xmm5,%%xmm4\n" \
      "packssdw  %%xmm7,%%xmm6\n" \
      "packssdw  %%xmm6,%%xmm4\n" \
      "packsswb  %%xmm4,%%xmm4\n" \
      "movaps  %%xmm4,%%xmm0\n" \
      "psrldq   $4,%%xmm0\n" \
      "pslld    $4,%%xmm0\n" \
      "pxor    %%xmm0,%%xmm4\n"
      "movd    %%xmm4,%0\n" \
      "":"=&r"(output):"r"(state->xN),"r"(state->xP),"r"(gq58x4_sse_Consts));
      return output;
}

__device__ __host__ void gq58x4_get_sse_state_(gq58x4_state* state,gq58x4_sse_state* sse_state){
  int i; for(i=0;i<8;i++) {sse_state->xN[i]=state->xN[i]; sse_state->xP[i]=state->xP[i];}
}


__device__ __host__ lt gq58x4_mod_g(lt x){ // returns x (mod g)
  lt F,G; F = (x>>58); G = x-(F<<58)+(F<<29)+(F<<30);
  return ((G>=gq58x4_g) ? (G-gq58x4_g) : G);
}

__device__ __host__ lt gq58x4_MyMult(lt A,lt B){ // returns AB (mod gq58x4_g), where it is implied that A,B<gq58x4_g;
  lt A1,A0,B1,B0,curr,x,m;
  A1=A>>32; B1=B>>32; A0=A-(A1<<32)+(12*A1); B0=B-(B1<<32)+(12*B1);
  if(A0>>32) {A0-=4294967284ULL; A1++;}
  if(B0>>32) {B0-=4294967284ULL; B1++;}
  curr=A1*B0+B1*A0; m=curr>>26; x=curr-(m<<26);
  curr=((3*m+(x<<4))<<28)+(gq58x4_g-12*x)+(144*A1*B1)+(gq58x4_mod_g(A0*B0));
  return gq58x4_mod_g(curr);
}

__device__ __host__ lt gq58x4_CNext2(lt N,lt P,lt myk,lt myq){   // returns (myk*N-myq*P) (mod gq58x4_g)
  lt curr1,curr2;
  curr1=gq58x4_MyMult(myk,N); curr2=gq58x4_MyMult(myq,P);
  if(curr1>=curr2) return (curr1-curr2); else return (gq58x4_g+curr1-curr2);
}

__device__ __host__ lt gq58x4_CNext(lt N,lt P){ // returns (8N-48P) (mod gq58x4_g)
  return gq58x4_mod_g((N+6*(gq58x4_g-P))<<3);
}


__device__ __host__ lt gq58x4_GetNextN(lt x0,lt x1,unsigned int n){ //returns x_{2^n}
  lt myk=gq58x4_k,myq=gq58x4_q,i,x=x1;
  for(i=0;i<n;i++){
    x=gq58x4_CNext2(x,x0,myk,myq);
    myk=gq58x4_CNext2(myk,2,myk,myq);
    myq=gq58x4_CNext2(myq,0,myq,0);
  }
  return x;
}

__device__ __host__ lt gq58x4_GetNextAny(lt x0,lt x1,lt N64,lt N0){ //N=2^64*N64+N0+1
  lt i,xp=x0,xn=x1,xpnew,xnnew,shift=0;
  i=N0; while(i>0){
    if(i%2==1){                        // xp,xn ----> 2^shift
      xpnew=gq58x4_GetNextN(xp,xn,shift);
      xnnew=gq58x4_GetNextN(xn,gq58x4_CNext(xn,xp),shift);
      xp=xpnew; xn=xnnew;
    }
    i/=2; shift++;
  }
  i=N64; shift=64; while(i>0){
    if(i%2==1){                        // xp,xn ----> 2^shift
      xpnew=gq58x4_GetNextN(xp,xn,shift);
      xnnew=gq58x4_GetNextN(xn,gq58x4_CNext(xn,xp),shift);
      xp=xpnew; xn=xnnew;
    }
    i/=2; shift++;
  }
  return xp;                       // returns x_N, where N=2^64*N64+N0+1
}

__device__ __host__ void gq58x4_skipahead_(gq58x4_state* state, lt offset64, lt offset0){ // offset=offset64*2^64+offset0+1
  lt xn,xp,j; 
  for(j=0;j<8;j++){
    xp=gq58x4_GetNextAny(state->xP[j],state->xN[j],offset64,offset0);
    xn=gq58x4_GetNextAny(state->xP[j],state->xN[j],offset64,offset0+1);
    state->xP[j]=xp; state->xN[j]=xn;
  }
}

__device__ __host__ void gq58x4_init_(gq58x4_state* state){
  lt x0=100152853817629549ULL,x1=132388305121829306ULL,xp,xn,j;
  for(j=0;j<8;j++){
    xp=gq58x4_GetNextAny(x0,x1,0,35048736516210783ULL);
    xn=gq58x4_GetNextAny(x0,x1,0,35048736516210784ULL);
    state->xP[j]=xp; state->xN[j]=xn; x0=xp; x1=xn;
  }
}

__device__ __host__ void gq58x4_init_short_sequence_(gq58x4_state* state,unsigned SequenceNumber){
  gq58x4_init_(state);                     // 0 <= SequenceNumber < 3*10^8;   length of each sequence <= 8*10^7
  gq58x4_skipahead_(state,0,82927047ULL*(unsigned long long)SequenceNumber);
}

__device__ __host__ void gq58x4_init_medium_sequence_(gq58x4_state* state,unsigned SequenceNumber){
  gq58x4_init_(state);                     // 0 <= SequenceNumber < 3*10^6;   length of each sequence <= 8*10^9
  gq58x4_skipahead_(state,0,8799201913ULL*(unsigned long long)SequenceNumber);
}

__device__ __host__ void gq58x4_init_long_sequence_(gq58x4_state* state,unsigned SequenceNumber){
  gq58x4_init_(state);                     // 0 <= SequenceNumber < 3*10^4;   length of each sequence <= 8*10^11
  gq58x4_skipahead_(state,0,828317697521ULL*(unsigned long long)SequenceNumber);
}

__device__ __host__ unsigned int gq58x4_generate_(gq58x4_state* state){
  unsigned int sum=0; int i; lt temp;
  for(i=0;i<8;i++){ 
    temp=gq58x4_mod_g((state->xN[i]+6*(gq58x4_g-state->xP[i]))<<3);
    state->xP[i]=state->xN[i]; state->xN[i]=temp;
    sum+= ((temp/gq58x4_gdiv16)<<((i<4)?(8*i):(8*i-28)));
  }
  return sum;
}

__device__ __host__ float gq58x4_generate_uniform_float_(gq58x4_state* state){
  unsigned int sum=0; int i; lt temp;
  for(i=0;i<8;i++){ 
    temp=gq58x4_mod_g((state->xN[i]+6*(gq58x4_g-state->xP[i]))<<3);
    state->xP[i]=state->xN[i]; state->xN[i]=temp;
    sum+= ((temp/gq58x4_gdiv16)<<((i<4)?(8*i):(8*i-28)));
  }
  return ((float) sum) * 2.3283064365386963e-10;
}

__host__ void gq58x4_print_state_(gq58x4_state* state){int i;
    printf("Generator State:\nxN={");
    for(i=0;i<8;i++) {printf("%llu",state->xN[i]%gq58x4_g); printf((i<7)?",":"}\nxP={");}
    for(i=0;i<8;i++) {printf("%llu",state->xP[i]%gq58x4_g); printf((i<7)?",":"}\n\n");}
}

__host__ void gq58x4_print_sse_state_(gq58x4_sse_state* state){int i;
    printf("Generator State:\nxN={");
    for(i=0;i<8;i++) {printf("%llu",state->xN[i]%gq58x4_g); printf((i<7)?",":"}\nxP={");}
    for(i=0;i<8;i++) {printf("%llu",state->xP[i]%gq58x4_g); printf((i<7)?",":"}\n\n");}
}

__global__ void gq58x4_kernel_generate_array(gq58x4_state* state, unsigned int* out, long* length) {
    unsigned sum,i,j,orbit,seqNum,shift; long offset; lt temp;

    __shared__ lt xP[gq58x4_THREADS];        // one generator per s=8 threads, i.e. one orbit
    __shared__ lt xN[gq58x4_THREADS];        // per thread, i.e. blockDim.x orbits per block
    __shared__ unsigned  a[gq58x4_THREADS];  // array "a" contains corresponding parts of output

    orbit   = threadIdx.x % 8;
    seqNum  = (threadIdx.x + blockIdx.x * blockDim.x)>>3;  // RNG_sequence index
    offset  = seqNum*(*length);                            // start of the section in the output array

    xP[threadIdx.x]=gq58x4_GetNextAny(state->xP[orbit],state->xN[orbit],0,offset);
    xN[threadIdx.x]=gq58x4_GetNextAny(state->xP[orbit],state->xN[orbit],0,offset+1);

    shift=((orbit<4)?(8*orbit):(8*orbit-28));

    for(i=0;i<(*length);i++){

      temp = gq58x4_CNext( xN[threadIdx.x], xP[threadIdx.x] );
      xP[threadIdx.x] = xN[threadIdx.x]; xN[threadIdx.x] = temp; 
      a[threadIdx.x]  = ((temp/gq58x4_gdiv16)<<shift);

      __syncthreads();              // each s=8 threads result in "length" values in the output array

      if(orbit==0){ sum=0; for(j=0;j<8;j++) sum+=a[threadIdx.x+j]; out[offset+i]=sum; }

    }
}

__host__ void gq58x4_generate_gpu_array_(gq58x4_state* state, unsigned int* dev_out, long length){

   long          mylength = length/gq58x4_ARRAY_SECTIONS;
   gq58x4_state*   dev_state;
   long*         dev_length;

   if((mylength*gq58x4_ARRAY_SECTIONS)<length) mylength++;

   gq58x4_CUDA_CALL(cudaMalloc((void**)&dev_state,sizeof(gq58x4_state)));
   gq58x4_CUDA_CALL(cudaMalloc((void**)&dev_length,sizeof(long)));
   gq58x4_CUDA_CALL(cudaMemcpy(dev_state,state,sizeof(gq58x4_state),cudaMemcpyHostToDevice));
   gq58x4_CUDA_CALL(cudaMemcpy(dev_length,&mylength,sizeof(long),cudaMemcpyHostToDevice));

   gq58x4_kernel_generate_array<<<gq58x4_BLOCKS,gq58x4_THREADS>>>(dev_state,dev_out,dev_length);
   gq58x4_CUDA_CALL(cudaGetLastError());
   
   gq58x4_CUDA_CALL(cudaFree(dev_state)); gq58x4_CUDA_CALL(cudaFree(dev_length));

}

__global__ void gq58x4_kernel_generate_array_float(gq58x4_state* state, float* out, long* length) {
    unsigned sum,i,j,orbit,seqNum,shift; long offset; lt temp;

    __shared__ lt xP[gq58x4_THREADS];        // one generator per s=8 threads, i.e. one orbit
    __shared__ lt xN[gq58x4_THREADS];        // per thread, i.e. blockDim.x orbits per block
    __shared__ unsigned  a[gq58x4_THREADS];  // array "a" contains corresponding parts of output

    orbit   = threadIdx.x % 8;
    seqNum  = (threadIdx.x + blockIdx.x * blockDim.x)>>3;  // RNG_sequence index
    offset  = seqNum*(*length);                            // start of the section in the output array

    xP[threadIdx.x]=gq58x4_GetNextAny(state->xP[orbit],state->xN[orbit],0,offset);
    xN[threadIdx.x]=gq58x4_GetNextAny(state->xP[orbit],state->xN[orbit],0,offset+1);

    shift=((orbit<4)?(8*orbit):(8*orbit-28));

    for(i=0;i<(*length);i++){

      temp = gq58x4_CNext( xN[threadIdx.x], xP[threadIdx.x] );
      xP[threadIdx.x] = xN[threadIdx.x]; xN[threadIdx.x] = temp; 
      a[threadIdx.x]  = ((temp/gq58x4_gdiv16)<<shift);

      __syncthreads();              // each s=8 threads result in "length" values in the output array

      if(orbit==0){ sum=0; for(j=0;j<8;j++) sum+=a[threadIdx.x+j]; out[offset+i]=((float)sum)* 2.3283064365386963e-10; }

    }
}

__host__ void gq58x4_generate_gpu_array_float_(gq58x4_state* state, float* dev_out, long length){

   long          mylength = length/gq58x4_ARRAY_SECTIONS;
   gq58x4_state*   dev_state;
   long*         dev_length;

   if((mylength*gq58x4_ARRAY_SECTIONS)<length) mylength++;

   gq58x4_CUDA_CALL(cudaMalloc((void**)&dev_state,sizeof(gq58x4_state)));
   gq58x4_CUDA_CALL(cudaMalloc((void**)&dev_length,sizeof(long)));
   gq58x4_CUDA_CALL(cudaMemcpy(dev_state,state,sizeof(gq58x4_state),cudaMemcpyHostToDevice));
   gq58x4_CUDA_CALL(cudaMemcpy(dev_length,&mylength,sizeof(long),cudaMemcpyHostToDevice));

   gq58x4_kernel_generate_array_float<<<gq58x4_BLOCKS,gq58x4_THREADS>>>(dev_state,dev_out,dev_length);
   gq58x4_CUDA_CALL(cudaGetLastError());
   
   gq58x4_CUDA_CALL(cudaFree(dev_state)); gq58x4_CUDA_CALL(cudaFree(dev_length));

}

__global__ void gq58x4_kernel_generate_array_double(gq58x4_state* state, double* out, long* length) {
    unsigned sum,i,j,orbit,seqNum,shift; long offset; lt temp;

    __shared__ lt xP[gq58x4_THREADS];        // one generator per s=8 threads, i.e. one orbit
    __shared__ lt xN[gq58x4_THREADS];        // per thread, i.e. blockDim.x orbits per block
    __shared__ unsigned  a[gq58x4_THREADS];  // array "a" contains corresponding parts of output

    orbit   = threadIdx.x % 8;
    seqNum  = (threadIdx.x + blockIdx.x * blockDim.x)>>3;  // RNG_sequence index
    offset  = seqNum*(*length);                            // start of the section in the output array

    xP[threadIdx.x]=gq58x4_GetNextAny(state->xP[orbit],state->xN[orbit],0,offset);
    xN[threadIdx.x]=gq58x4_GetNextAny(state->xP[orbit],state->xN[orbit],0,offset+1);

    shift=((orbit<4)?(8*orbit):(8*orbit-28));

    for(i=0;i<(*length);i++){

      temp = gq58x4_CNext( xN[threadIdx.x], xP[threadIdx.x] );
      xP[threadIdx.x] = xN[threadIdx.x]; xN[threadIdx.x] = temp; 
      a[threadIdx.x]  = ((temp/gq58x4_gdiv16)<<shift);

      __syncthreads();              // each s=8 threads result in "length" values in the output array

      if(orbit==0){ sum=0; for(j=0;j<8;j++) sum+=a[threadIdx.x+j]; out[offset+i]=((double)sum)* 2.3283064365386963e-10; }

    }
}

__host__ void gq58x4_generate_gpu_array_double_(gq58x4_state* state, double* dev_out, long length){

   long          mylength = length/gq58x4_ARRAY_SECTIONS;
   gq58x4_state*   dev_state;
   long*         dev_length;

   if((mylength*gq58x4_ARRAY_SECTIONS)<length) mylength++;

   gq58x4_CUDA_CALL(cudaMalloc((void**)&dev_state,sizeof(gq58x4_state)));
   gq58x4_CUDA_CALL(cudaMalloc((void**)&dev_length,sizeof(long)));
   gq58x4_CUDA_CALL(cudaMemcpy(dev_state,state,sizeof(gq58x4_state),cudaMemcpyHostToDevice));
   gq58x4_CUDA_CALL(cudaMemcpy(dev_length,&mylength,sizeof(long),cudaMemcpyHostToDevice));

   gq58x4_kernel_generate_array_double<<<gq58x4_BLOCKS,gq58x4_THREADS>>>(dev_state,dev_out,dev_length);
   gq58x4_CUDA_CALL(cudaGetLastError());
   
   gq58x4_CUDA_CALL(cudaFree(dev_state)); gq58x4_CUDA_CALL(cudaFree(dev_length));

}

__host__ void gq58x4_generate_array_(gq58x4_state* state, unsigned int* out, long length){

   long          mylength = length/gq58x4_ARRAY_SECTIONS;
   gq58x4_state*   dev_state;
   unsigned int* dev_out;
   long*         dev_length;

   if((mylength*gq58x4_ARRAY_SECTIONS)<length) mylength++;

   gq58x4_CUDA_CALL(cudaMalloc((void**)&dev_state,sizeof(gq58x4_state)));
   gq58x4_CUDA_CALL(cudaMalloc((void**)&dev_out,mylength*gq58x4_ARRAY_SECTIONS*sizeof(unsigned int)));
   gq58x4_CUDA_CALL(cudaMalloc((void**)&dev_length,sizeof(long)));
   gq58x4_CUDA_CALL(cudaMemcpy(dev_state,state,sizeof(gq58x4_state),cudaMemcpyHostToDevice));
   gq58x4_CUDA_CALL(cudaMemcpy(dev_length,&mylength,sizeof(long),cudaMemcpyHostToDevice));

   gq58x4_kernel_generate_array<<<gq58x4_BLOCKS,gq58x4_THREADS>>>(dev_state,dev_out,dev_length);
   gq58x4_CUDA_CALL(cudaGetLastError());
   
   gq58x4_CUDA_CALL(cudaMemcpy(out,dev_out,length*sizeof(unsigned int),cudaMemcpyDeviceToHost));
   gq58x4_CUDA_CALL(cudaFree(dev_state)); gq58x4_CUDA_CALL(cudaFree(dev_out)); 
   gq58x4_CUDA_CALL(cudaFree(dev_length));

}