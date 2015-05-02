C----------------------------------------------------------------------
      LOGICAL FUNCTION CUT(XMUET)
C----------------------------------------------------------------------
C
C Determine si un evenement passe ou non les coupures
C
C----------------------------------------------------------------------
!     $ use OMP_LIB
      IMPLICIT NONE
      double precision XMUET
      INTEGER I
      INCLUDE 'ppp.for'
      INCLUDE 'angle.for'
      INCLUDE 'cutpar.for'

C   COUPURE S'IL Y A LIEU
      CUT=.FALSE.
      DO I=1,INP
        CUT=CUT.OR.(POUT(4,I).LT.EMIN)
      ENDDO
      CUT=CUT.OR.(ABS(COS1P1).GT.ACUT).OR.(ABS(COS2P1).GT.ACUT).OR.
     &  (ABS(COS3P1).GT.ACUT).OR.
     &  (COS12.GT.BCUT).OR.(COS13.GT.BCUT).OR.(COS23.GT.BCUT)
     &  .OR.(ABS(COSN).LT.SINCUT)


      RETURN
      END
C----------------------------------------------------------------------
