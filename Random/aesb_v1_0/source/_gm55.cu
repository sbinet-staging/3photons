// (c) Copyright 2013 Lev Barash, Landau Institute for Theoretical Physics, Russian Academy of Sciences
// This is supplement to the paper:
// L.Yu. Barash, L.N. Shchur, "PRAND: GPU accelerated parallel random number generation library: Using most reliable algorithms and applying parallelism of modern GPUs and CPUs".
// e-mail: barash @ itp.ac.ru (remove space)

#include<stdio.h>

#define gm55_CUDA_CALL(x) do { if((x) != cudaSuccess) { printf("Error: %s at %s:%d\n",cudaGetErrorString(cudaGetLastError()),__FILE__,__LINE__); exit(1);}} while(0)

#define gm55_BLOCKS  128
#define gm55_THREADS 128
#define gm55_ARRAY_SECTIONS (gm55_BLOCKS*gm55_THREADS/8)

#define gm55_k 256
#define gm55_q 176
#define gm55_g       36028797018961904ULL
#define gm55_gdiv16  2251799813685119ULL

typedef unsigned long long lt;

typedef struct{
  lt xN[8] __attribute__ ((aligned(16))), 
     xP[8] __attribute__ ((aligned(16)));
} gm55_state;

typedef gm55_state gm55_sse_state;

lt gm55_sse_Consts[8] __attribute__ ((aligned(16))) = {396316767208580944ULL,396316767208580944ULL,
   2064ULL,2064ULL,36028792732385279ULL,36028792732385279ULL,36028797018961904ULL,36028797018961904ULL};

__host__ unsigned int gm55_sse_generate_(gm55_sse_state* state){
  unsigned output;
  asm volatile("movaps (%3),%%xmm0\n" \

      "movaps (%2),%%xmm1\n" \
      "movaps (%1),%%xmm4\n" \
      "movaps %%xmm4,(%2)\n" \
      "psllq  $4,%%xmm4\n" \
      "paddq  %%xmm0,%%xmm4\n" \
      "movaps %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "psllq  $3,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm4\n" \
      "movaps %%xmm4,%%xmm2\n" \
      "psrlq  $51,%%xmm2\n" \
      "movaps %%xmm2,%%xmm3\n" \
      "psllq  $7,%%xmm3\n" \
      "paddq  %%xmm2,%%xmm3\n" \
      "psllq  $51,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm4\n" \
      "paddq  %%xmm3,%%xmm4\n" \
      "psllq  $4,%%xmm4\n" \
      "movaps %%xmm4,%%xmm1\n" \
      "paddq  16(%3),%%xmm1\n" \
      "pshufd $245,%%xmm1,%%xmm3\n" \
      "pcmpgtd 32(%3),%%xmm3\n" \
      "pand    48(%3),%%xmm3\n" \
      "psubq   %%xmm3,%%xmm4\n" \
      "movaps %%xmm4,(%1)\n" \

      "movaps 16(%2),%%xmm1\n" \
      "movaps 16(%1),%%xmm5\n" \
      "movaps %%xmm5,16(%2)\n" \
      "psllq  $4,%%xmm5\n" \
      "paddq  %%xmm0,%%xmm5\n" \
      "movaps %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "psllq  $3,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm5\n" \
      "movaps %%xmm5,%%xmm2\n" \
      "psrlq  $51,%%xmm2\n" \
      "movaps %%xmm2,%%xmm3\n" \
      "psllq  $7,%%xmm3\n" \
      "paddq  %%xmm2,%%xmm3\n" \
      "psllq  $51,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm5\n" \
      "paddq  %%xmm3,%%xmm5\n" \
      "psllq  $4,%%xmm5\n" \
      "movaps %%xmm5,%%xmm1\n" \
      "paddq  16(%3),%%xmm1\n" \
      "pshufd $245,%%xmm1,%%xmm3\n" \
      "pcmpgtd 32(%3),%%xmm3\n" \
      "pand    48(%3),%%xmm3\n" \
      "psubq   %%xmm3,%%xmm5\n" \
      "movaps %%xmm5,16(%1)\n" \

      "movaps 32(%2),%%xmm1\n" \
      "movaps 32(%1),%%xmm6\n" \
      "movaps %%xmm6,32(%2)\n" \
      "psllq  $4,%%xmm6\n" \
      "paddq  %%xmm0,%%xmm6\n" \
      "movaps %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "psllq  $3,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm6\n" \
      "movaps %%xmm6,%%xmm2\n" \
      "psrlq  $51,%%xmm2\n" \
      "movaps %%xmm2,%%xmm3\n" \
      "psllq  $7,%%xmm3\n" \
      "paddq  %%xmm2,%%xmm3\n" \
      "psllq  $51,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm6\n" \
      "paddq  %%xmm3,%%xmm6\n" \
      "psllq  $4,%%xmm6\n" \
      "movaps %%xmm6,%%xmm1\n" \
      "paddq  16(%3),%%xmm1\n" \
      "pshufd $245,%%xmm1,%%xmm3\n" \
      "pcmpgtd 32(%3),%%xmm3\n" \
      "pand    48(%3),%%xmm3\n" \
      "psubq   %%xmm3,%%xmm6\n" \
      "movaps %%xmm6,32(%1)\n" \

      "movaps 48(%2),%%xmm1\n" \
      "movaps 48(%1),%%xmm7\n" \
      "movaps %%xmm7,48(%2)\n" \
      "psllq  $4,%%xmm7\n" \
      "paddq  %%xmm0,%%xmm7\n" \
      "movaps %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "psllq  $3,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm7\n" \
      "movaps %%xmm7,%%xmm2\n" \
      "psrlq  $51,%%xmm2\n" \
      "movaps %%xmm2,%%xmm3\n" \
      "psllq  $7,%%xmm3\n" \
      "paddq  %%xmm2,%%xmm3\n" \
      "psllq  $51,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm7\n" \
      "paddq  %%xmm3,%%xmm7\n" \
      "psllq  $4,%%xmm7\n" \
      "movaps %%xmm7,%%xmm1\n" \
      "paddq  16(%3),%%xmm1\n" \
      "pshufd $245,%%xmm1,%%xmm3\n" \
      "pcmpgtd 32(%3),%%xmm3\n" \
      "pand    48(%3),%%xmm3\n" \
      "psubq   %%xmm3,%%xmm7\n" \
      "movaps %%xmm7,48(%1)\n" \

      "psrlq  $51,%%xmm4\n" \
      "psrlq  $51,%%xmm5\n" \
      "psrlq  $51,%%xmm6\n" \
      "psrlq  $51,%%xmm7\n" \
      "packssdw  %%xmm5,%%xmm4\n" \
      "packssdw  %%xmm7,%%xmm6\n" \
      "packssdw  %%xmm6,%%xmm4\n" \
      "packsswb  %%xmm4,%%xmm4\n" \
      "movaps  %%xmm4,%%xmm0\n" \
      "psrldq   $4,%%xmm0\n" \
      "pslld    $4,%%xmm0\n" \
      "pxor    %%xmm0,%%xmm4\n"
      "movd    %%xmm4,%0\n" \
      "":"=&r"(output):"r"(state->xN),"r"(state->xP),"r"(gm55_sse_Consts));
      return output;
}

__device__ __host__ void gm55_get_sse_state_(gm55_state* state,gm55_sse_state* sse_state){
  int i; for(i=0;i<8;i++) {sse_state->xN[i]=state->xN[i]; sse_state->xP[i]=state->xP[i];}
}

__device__ __host__ lt gm55_mod_g(lt x){ // returns x (mod g)
  lt F,G; F = (x>>55); G = x-(F<<55)+(2064*F);
  return ((G>=gm55_g) ? (G-gm55_g) : G);
}

__device__ __host__ lt gm55_MyMult(lt A,lt B){ // returns AB (mod gm55_g), where it is implied that A,B<gm55_g;
  lt A1,A0,B1,B0,curr,x,m;
  A1=A>>28; B1=B>>27; A0=A-(A1<<28); B0=B-(B1<<27);
  curr=2*A1*B0+B1*A0; m=curr>>28; x=curr-(m<<28);
  curr=(x<<27)+2064*m+(gm55_mod_g(129*A1*B1)<<4)+A0*B0;
  return gm55_mod_g(curr);
}

__device__ __host__ lt gm55_CNext2(lt N,lt P,lt myk,lt myq){   // returns (myk*N-myq*P) (mod gm55_g)
  lt curr1,curr2;
  curr1=gm55_MyMult(myk,N); curr2=gm55_MyMult(myq,P);
  if(curr1>=curr2) return (curr1-curr2); else return (gm55_g+curr1-curr2);
}

__device__ __host__ lt gm55_CNext(lt N,lt P){ // returns (256*N-176*P) (mod gm55_g)
  return gm55_mod_g((N<<8)+176*(gm55_g-P));
}

__device__ __host__ lt gm55_GetNextN(lt x0,lt x1,unsigned int n){ //returns x_{2^n}
  lt myk=gm55_k,myq=gm55_q,i,x=x1;
  for(i=0;i<n;i++){
    x=gm55_CNext2(x,x0,myk,myq);
    myk=gm55_CNext2(myk,2,myk,myq);
    myq=gm55_CNext2(myq,0,myq,0);
  }
  return x;
}

__device__ __host__ lt gm55_GetNextAny(lt x0,lt x1,lt N64,lt N0){ //N=2^64*N64+N0+1
  lt i,xp=x0,xn=x1,xpnew,xnnew,shift=0;
  i=N0; while(i>0){
    if(i%2==1){                        // xp,xn ----> 2^shift
      xpnew=gm55_GetNextN(xp,xn,shift);
      xnnew=gm55_GetNextN(xn,gm55_CNext(xn,xp),shift);
      xp=xpnew; xn=xnnew;
    }
    i/=2; shift++;
  }
  i=N64; shift=64; while(i>0){
    if(i%2==1){                        // xp,xn ----> 2^shift
      xpnew=gm55_GetNextN(xp,xn,shift);
      xnnew=gm55_GetNextN(xn,gm55_CNext(xn,xp),shift);
      xp=xpnew; xn=xnnew;
    }
    i/=2; shift++;
  }
  return xp;                       // returns x_N, where N=2^64*N64+N0+1
}

__device__ __host__ void gm55_skipahead_(gm55_state* state, lt offset64, lt offset0){ // offset=offset64*2^64+offset0+1
  lt xn,xp,j; 
  for(j=0;j<8;j++){
    xp=gm55_GetNextAny(state->xP[j],state->xN[j],offset64,offset0);
    xn=gm55_GetNextAny(state->xP[j],state->xN[j],offset64,offset0+1);
    state->xP[j]=xp; state->xN[j]=xn;
  }
}

__device__ __host__ void gm55_init_(gm55_state* state){
  lt x0=gm55_mod_g(100152853817629549ULL),x1=gm55_mod_g(132388305121829306ULL),xp,xn,j;
  for(j=0;j<8;j++){
    xp=gm55_GetNextAny(x0,x1,7730941120ULL,2741045636569588180ULL);
    xn=gm55_GetNextAny(x0,x1,7730941120ULL,2741045636569588181ULL);
    state->xP[j]=xp; state->xN[j]=xn; x0=xp; x1=xn;
  }
}

__device__ __host__ void gm55_init_short_sequence_(gm55_state* state,lt SequenceNumber){ // 0 <= SequenceNumber < 10^18
  lt n1,n2;                                                                // length of each sequence < 10^10
  gm55_init_(state);
  n1=SequenceNumber/892447987; n2=SequenceNumber%892447987;
  gm55_skipahead_(state,n1,n1*4193950067);
  gm55_skipahead_(state,0,n2*20669825409); // thus we are skipping ahead (SequenceNumber*20669825409) numbers
}

__device__ __host__ void gm55_init_long_sequence_(gm55_state* state,lt SequenceNumber){ // 0 <= SequenceNumber < 4*10^9
  gm55_init_(state);                                                      // length of each sequence  <   10^20
  gm55_skipahead_(state,8*SequenceNumber,2699204111*SequenceNumber);
}

__device__ __host__ unsigned int gm55_generate_(gm55_state* state){
  unsigned int sum=0; int i; lt temp;
  for(i=0;i<8;i++){ 
    temp=gm55_mod_g(((state->xN[i])<<8)+gm55_q*(gm55_g-state->xP[i]));
    state->xP[i]=state->xN[i]; state->xN[i]=temp;
    sum+= ((temp/gm55_gdiv16)<<((i<4)?(8*i):(8*i-28)));
  }
  return sum;
}

__device__ __host__ float gm55_generate_uniform_float_(gm55_state* state){
  unsigned int sum=0; int i; lt temp;
  for(i=0;i<8;i++){ 
    temp=gm55_mod_g(((state->xN[i])<<8)+gm55_q*(gm55_g-state->xP[i]));
    state->xP[i]=state->xN[i]; state->xN[i]=temp;
    sum+= ((temp/gm55_gdiv16)<<((i<4)?(8*i):(8*i-28)));
  }
  return ((float) sum) * 2.3283064365386963e-10;
}

__host__ void gm55_print_state_(gm55_state* state){int i;
    printf("Generator State:\nxN={");
    for(i=0;i<8;i++) {printf("%llu",state->xN[i]%gm55_g); printf((i<7)?",":"}\nxP={");}
    for(i=0;i<8;i++) {printf("%llu",state->xP[i]%gm55_g); printf((i<7)?",":"}\n\n");}
}

__host__ void gm55_print_sse_state_(gm55_sse_state* state){int i;
    printf("Generator State:\nxN={");
    for(i=0;i<8;i++) {printf("%llu",state->xN[i]%gm55_g); printf((i<7)?",":"}\nxP={");}
    for(i=0;i<8;i++) {printf("%llu",state->xP[i]%gm55_g); printf((i<7)?",":"}\n\n");}
}

__global__ void gm55_kernel_generate_array(gm55_state* state, unsigned int* out, long* length) {
    unsigned sum,i,j,orbit,seqNum,shift; long offset; lt temp;

    __shared__ lt xP[gm55_THREADS];        // one generator per s=8 threads, i.e. one orbit
    __shared__ lt xN[gm55_THREADS];        // per thread, i.e. blockDim.x orbits per block
    __shared__ unsigned  a[gm55_THREADS];  // array "a" contains corresponding parts of output

    orbit   = threadIdx.x % 8;
    seqNum  = (threadIdx.x + blockIdx.x * blockDim.x)>>3;  // RNG_sequence index
    offset  = seqNum*(*length);                            // start of the section in the output array

    xP[threadIdx.x]=gm55_GetNextAny(state->xP[orbit],state->xN[orbit],0,offset);
    xN[threadIdx.x]=gm55_GetNextAny(state->xP[orbit],state->xN[orbit],0,offset+1);

    shift=((orbit<4)?(8*orbit):(8*orbit-28));

    for(i=0;i<(*length);i++){

      temp = gm55_CNext( xN[threadIdx.x], xP[threadIdx.x] );
      xP[threadIdx.x] = xN[threadIdx.x]; xN[threadIdx.x] = temp; 
      a[threadIdx.x]  = ((temp/gm55_gdiv16)<<shift);

      __syncthreads();              // each s=8 threads result in "length" values in the output array

      if(orbit==0){ sum=0; for(j=0;j<8;j++) sum+=a[threadIdx.x+j]; out[offset+i]=sum; }

    }
}

__host__ void gm55_generate_gpu_array_(gm55_state* state, unsigned int* dev_out, long length){

   long          mylength = length/gm55_ARRAY_SECTIONS;
   gm55_state*   dev_state;
   long*         dev_length;

   if((mylength*gm55_ARRAY_SECTIONS)<length) mylength++;

   gm55_CUDA_CALL(cudaMalloc((void**)&dev_state,sizeof(gm55_state)));
   gm55_CUDA_CALL(cudaMalloc((void**)&dev_length,sizeof(long)));
   gm55_CUDA_CALL(cudaMemcpy(dev_state,state,sizeof(gm55_state),cudaMemcpyHostToDevice));
   gm55_CUDA_CALL(cudaMemcpy(dev_length,&mylength,sizeof(long),cudaMemcpyHostToDevice));

   gm55_kernel_generate_array<<<gm55_BLOCKS,gm55_THREADS>>>(dev_state,dev_out,dev_length);
   gm55_CUDA_CALL(cudaGetLastError());
   
   gm55_CUDA_CALL(cudaFree(dev_state)); gm55_CUDA_CALL(cudaFree(dev_length));

}

__global__ void gm55_kernel_generate_array_float(gm55_state* state, float* out, long* length) {
    unsigned sum,i,j,orbit,seqNum,shift; long offset; lt temp;

    __shared__ lt xP[gm55_THREADS];        // one generator per s=8 threads, i.e. one orbit
    __shared__ lt xN[gm55_THREADS];        // per thread, i.e. blockDim.x orbits per block
    __shared__ unsigned  a[gm55_THREADS];  // array "a" contains corresponding parts of output

    orbit   = threadIdx.x % 8;
    seqNum  = (threadIdx.x + blockIdx.x * blockDim.x)>>3;  // RNG_sequence index
    offset  = seqNum*(*length);                            // start of the section in the output array

    xP[threadIdx.x]=gm55_GetNextAny(state->xP[orbit],state->xN[orbit],0,offset);
    xN[threadIdx.x]=gm55_GetNextAny(state->xP[orbit],state->xN[orbit],0,offset+1);

    shift=((orbit<4)?(8*orbit):(8*orbit-28));

    for(i=0;i<(*length);i++){

      temp = gm55_CNext( xN[threadIdx.x], xP[threadIdx.x] );
      xP[threadIdx.x] = xN[threadIdx.x]; xN[threadIdx.x] = temp; 
      a[threadIdx.x]  = ((temp/gm55_gdiv16)<<shift);

      __syncthreads();              // each s=8 threads result in "length" values in the output array

      if(orbit==0){ sum=0; for(j=0;j<8;j++) sum+=a[threadIdx.x+j]; out[offset+i]=((float)sum)* 2.3283064365386963e-10; }

    }
}

__host__ void gm55_generate_gpu_array_float_(gm55_state* state, float* dev_out, long length){

   long          mylength = length/gm55_ARRAY_SECTIONS;
   gm55_state*   dev_state;
   long*         dev_length;

   if((mylength*gm55_ARRAY_SECTIONS)<length) mylength++;

   gm55_CUDA_CALL(cudaMalloc((void**)&dev_state,sizeof(gm55_state)));
   gm55_CUDA_CALL(cudaMalloc((void**)&dev_length,sizeof(long)));
   gm55_CUDA_CALL(cudaMemcpy(dev_state,state,sizeof(gm55_state),cudaMemcpyHostToDevice));
   gm55_CUDA_CALL(cudaMemcpy(dev_length,&mylength,sizeof(long),cudaMemcpyHostToDevice));

   gm55_kernel_generate_array_float<<<gm55_BLOCKS,gm55_THREADS>>>(dev_state,dev_out,dev_length);
   gm55_CUDA_CALL(cudaGetLastError());
   
   gm55_CUDA_CALL(cudaFree(dev_state)); gm55_CUDA_CALL(cudaFree(dev_length));

}

__global__ void gm55_kernel_generate_array_double(gm55_state* state, double* out, long* length) {
    unsigned sum,i,j,orbit,seqNum,shift; long offset; lt temp;

    __shared__ lt xP[gm55_THREADS];        // one generator per s=8 threads, i.e. one orbit
    __shared__ lt xN[gm55_THREADS];        // per thread, i.e. blockDim.x orbits per block
    __shared__ unsigned  a[gm55_THREADS];  // array "a" contains corresponding parts of output

    orbit   = threadIdx.x % 8;
    seqNum  = (threadIdx.x + blockIdx.x * blockDim.x)>>3;  // RNG_sequence index
    offset  = seqNum*(*length);                            // start of the section in the output array

    xP[threadIdx.x]=gm55_GetNextAny(state->xP[orbit],state->xN[orbit],0,offset);
    xN[threadIdx.x]=gm55_GetNextAny(state->xP[orbit],state->xN[orbit],0,offset+1);

    shift=((orbit<4)?(8*orbit):(8*orbit-28));

    for(i=0;i<(*length);i++){

      temp = gm55_CNext( xN[threadIdx.x], xP[threadIdx.x] );
      xP[threadIdx.x] = xN[threadIdx.x]; xN[threadIdx.x] = temp; 
      a[threadIdx.x]  = ((temp/gm55_gdiv16)<<shift);

      __syncthreads();              // each s=8 threads result in "length" values in the output array

      if(orbit==0){ sum=0; for(j=0;j<8;j++) sum+=a[threadIdx.x+j]; out[offset+i]=((double)sum)* 2.3283064365386963e-10; }

    }
}

__host__ void gm55_generate_gpu_array_double_(gm55_state* state, double* dev_out, long length){

   long          mylength = length/gm55_ARRAY_SECTIONS;
   gm55_state*   dev_state;
   long*         dev_length;

   if((mylength*gm55_ARRAY_SECTIONS)<length) mylength++;

   gm55_CUDA_CALL(cudaMalloc((void**)&dev_state,sizeof(gm55_state)));
   gm55_CUDA_CALL(cudaMalloc((void**)&dev_length,sizeof(long)));
   gm55_CUDA_CALL(cudaMemcpy(dev_state,state,sizeof(gm55_state),cudaMemcpyHostToDevice));
   gm55_CUDA_CALL(cudaMemcpy(dev_length,&mylength,sizeof(long),cudaMemcpyHostToDevice));

   gm55_kernel_generate_array_double<<<gm55_BLOCKS,gm55_THREADS>>>(dev_state,dev_out,dev_length);
   gm55_CUDA_CALL(cudaGetLastError());
   
   gm55_CUDA_CALL(cudaFree(dev_state)); gm55_CUDA_CALL(cudaFree(dev_length));

}

__host__ void gm55_generate_array_(gm55_state* state, unsigned int* out, long length){

   long          mylength = length/gm55_ARRAY_SECTIONS;
   gm55_state*   dev_state;
   unsigned int* dev_out;
   long*         dev_length;

   if((mylength*gm55_ARRAY_SECTIONS)<length) mylength++;

   gm55_CUDA_CALL(cudaMalloc((void**)&dev_state,sizeof(gm55_state)));
   gm55_CUDA_CALL(cudaMalloc((void**)&dev_out,mylength*gm55_ARRAY_SECTIONS*sizeof(unsigned int)));
   gm55_CUDA_CALL(cudaMalloc((void**)&dev_length,sizeof(long)));
   gm55_CUDA_CALL(cudaMemcpy(dev_state,state,sizeof(gm55_state),cudaMemcpyHostToDevice));
   gm55_CUDA_CALL(cudaMemcpy(dev_length,&mylength,sizeof(long),cudaMemcpyHostToDevice));

   gm55_kernel_generate_array<<<gm55_BLOCKS,gm55_THREADS>>>(dev_state,dev_out,dev_length);
   gm55_CUDA_CALL(cudaGetLastError());
   
   gm55_CUDA_CALL(cudaMemcpy(out,dev_out,length*sizeof(unsigned int),cudaMemcpyDeviceToHost));
   gm55_CUDA_CALL(cudaFree(dev_state)); gm55_CUDA_CALL(cudaFree(dev_out)); 
   gm55_CUDA_CALL(cudaFree(dev_length));

}
