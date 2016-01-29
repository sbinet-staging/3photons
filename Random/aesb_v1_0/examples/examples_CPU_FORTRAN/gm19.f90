
PROGRAM GM19_Example

   PARAMETER (N=100000000)
   INTEGER i,GM19_generate,sum
   REAL val,rsum

   TYPE GM19_state
     INTEGER z(68)
   END TYPE GM19_state
   TYPE(GM19_state) state

   CALL GM19_Init ( state )
   CALL GM19_Print_State (state)

   sum = 0

   DO i=1,N
     sum = sum + GM19_generate (state)
   END DO

   WRITE(*,1) N
   rsum = sum/4294967296.
   val = GM19_generate (state)/4294967296.
   IF (rsum<0) THEN
     rsum=rsum+1
   ENDIF
   IF (val<0) THEN
     val=val+1
   ENDIF
   WRITE(*,2) rsum
   WRITE(*,3) val

1  FORMAT(I9," GM19 pseudorandom numbers generated using general instructions of CPU.")
2  FORMAT("Fractional part of the total sum of generated numbers: ",F8.6)
3  FORMAT("Next output value: ",F8.6)

END PROGRAM GM19_Example

