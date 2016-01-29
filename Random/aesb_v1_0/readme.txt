

To compile the library:


  cd source
  make


To compile, for example, "mrg32k3a.cu" in the current directory, with use of the library:


  nvcc -O2 -arch=sm_20 -maxrregcount 32 -o mrg32k3a mrg32k3a.cu -I../../include -L../../lib -lprand


To compile all the examples in C in the current directory, with use of the library:


  bash -c 'for i in *.cu; do nvcc -O2 -arch=sm_20 -maxrregcount 32 -o ${i%???} $i -I../../include -L../../lib -lprand; done'


To compile all the examples in FORTRAN in the current directory, with use of the library:


  bash -c 'for i in *.f90; do nvcc -c -arch=sm_20 -maxrregcount 32 -O2 ${i%????}.cu; gfortran -o ${i%????} $i ${i%????}.o -L /usr/local/cuda/lib64 -lcudart; rm ${i%????}.o; done'
