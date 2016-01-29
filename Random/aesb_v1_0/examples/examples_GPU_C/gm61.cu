
#include<stdio.h>
#include<gm61.h>

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
   gm61_state state; 
   gm61_init_(&state); 
   gm61_print_state_(&state);
   gm61_generate_array_(&state,out,NN);
   printf("%ld GM61 pseudorandom numbers generated using GPGPU.\n",NN);
   PrintOutSum(out,NN-1);
   printf("Last output value: %f\n",out[NN-1]/4294967296.);
   free(out);
   return 0;
}

