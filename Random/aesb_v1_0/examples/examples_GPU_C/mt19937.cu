
#include<stdio.h>
#include<mt19937.h>

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
   mt19937_state state; 
   mt19937_init_device_consts_();
   mt19937_init_(&state); 
   mt19937_print_state_(&state);
   mt19937_generate_array_(&state,out,NN);
   printf("%ld MT19937 pseudorandom numbers generated using GPGPU.\n",NN);
   PrintOutSum(out,NN-1);
   printf("Last output value: %f\n",out[NN-1]/4294967296.);
   mt19937_free_device_consts_();
   free(out);
   return 0;
}

