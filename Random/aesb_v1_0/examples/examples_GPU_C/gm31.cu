
#include<stdio.h>
#include<gm31.h>

#define NN       100000001UL

void PrintOutSum(unsigned int * out, long length){
   long i; unsigned int  sum=0;
   for(i=0;i<length;i++) sum+=out[i];
   printf("Fractional part of the total sum of first %ld generated numbers: %f\n",length,sum/4294967296.);
}

int main(void){
   unsigned int* out;
   out=(unsigned int*)malloc(NN*sizeof(unsigned int));
   if(out==NULL){
     printf("Not enough memory"); return 1;
   }
   gm31_state state; 
   gm31_init_(&state); 
   gm31_print_state_(&state);
   gm31_generate_array_(&state,out,NN);
   printf("%ld GM31 pseudorandom numbers generated using GPGPU.\n",NN);
   PrintOutSum(out,NN-1);
   printf("Last output value: %f\n",out[NN-1]/4294967296.);
   free(out);
   return 0;
}

