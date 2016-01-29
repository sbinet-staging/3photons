
PROGRAM GQ58X1_Example

   PARAMETER (N=100000001)
   INTEGER i,sum,out(N)
   REAL val,rsum

   TYPE GQ58X1_state
     INTEGER z(132)
   END TYPE GQ58X1_state
   TYPE(GQ58X1_state) state

   CALL GQ58X1_Init ( state )
   CALL GQ58X1_Print_State (state)
   CALL GQ58X1_Generate_Array (state,out,N)

   sum = 0

   DO i=1,N-1
     sum = sum + out(i)
   END DO

   WRITE(*,1) N
   rsum = sum/4294967296.
   val = out(N)/4294967296.
   IF (rsum<0) THEN
     rsum=rsum+1
   ENDIF
   IF (val<0) THEN
     val=val+1
   ENDIF
   WRITE(*,2) N-1,rsum
   WRITE(*,3) val

1  FORMAT(I9," GQ58.1 pseudorandom numbers generated using GPGPU.")
2  FORMAT("Fractional part of the total sum of first ",I9," generated numbers: ",F8.6)
3  FORMAT("Last output value: ",F8.6)

END PROGRAM GQ58X1_Example

