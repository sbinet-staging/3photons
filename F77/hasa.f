C--------------------------------------------------------------------------
      REAL*8 FUNCTION RN(IDMY)
C--------------------------------------------------------------------------
C    ***************************************************************
C    * RANDOM NUMBER FUNCTION TAKEN FROM KNUTH RANF                *
C    * (SEMINUMERICAL ALGORITHMS).                                 *
C    * METHOD IS X(N)=MOD(X(N-55)-X(N-24),1/FMODUL)                *
C    * NO PROVISION YET FOR CONTROL OVER THE SEED NUMBER.          *
C    *                                                             *
C    * RANF GIVES ONE RANDOM NUMBER BETWEEN 0 AND 1.               *
C    * IRN55 GENERATES 55 RANDOM NUMBERS BETWEEN 0 AND 1/FMODUL.   *
C    * IN55  INITIALIZES THE 55 NUMBERS AND WARMS UP THE SEQUENCE. *
C    ***************************************************************
      IMPLICIT NONE
      INTEGER IA(55),NCALL,MCALL,IDMY
      REAL*8 FMODUL
      PARAMETER (FMODUL=1.D-09)
      DATA NCALL/0/
      DATA MCALL/55/
      SAVE IA, NCALL, MCALL

      IF (NCALL.EQ.0) THEN
        CALL IN55(IA,234612947)
        NCALL = 1
      ENDIF
      IF (MCALL.EQ.0) THEN
        CALL IRN55(IA)
        MCALL=55
      ENDIF
      RN=IA(MCALL)*FMODUL
      MCALL=MCALL-1
      END

C--------------------------------------------------------------------------
      SUBROUTINE IN55(IA,IX)
C--------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER IA(55),IX,J,K,I,II,MODULO
      PARAMETER (MODULO=1000000000)
      IA(55)=IX
      J=IX
      K=1
      DO I=1,54
        II=MOD(21*I,55)
        IA(II)=K
        K=J-K
        IF(K.LT.0)K=K+MODULO
        J=IA(II)
      ENDDO
      DO I=1,10
        CALL IRN55(IA)
      ENDDO
      RETURN
      END

C--------------------------------------------------------------------------
      SUBROUTINE IRN55(IA)
C--------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER IA(55),MODULO,I,J
      PARAMETER (MODULO=1000000000)
      DO I=1,24
        J=IA(I)-IA(I+31)
        IF(J.LT.0)J=J+MODULO
        IA(I)=J
      ENDDO
      DO I=25,55
        J=IA(I)-IA(I-24)
        IF(J.LT.0)J=J+MODULO
        IA(I)=J
      ENDDO
      RETURN
      END
