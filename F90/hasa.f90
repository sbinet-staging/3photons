!--------------------------------------------------------------------------
FUNCTION RN (IDMY)
!--------------------------------------------------------------------------
!    ***************************************************************
!    * RANDOM NUMBER FUNCTION TAKEN FROM KNUTH RANF                *
!    * (SEMINUMERICAL ALGORITHMS).                                 *
!    * METHOD IS X (N)=MOD (X (N-55)-X (N-24), 1/FMODUL)           *
!    * NO PROVISION YET FOR CONTROL OVER THE SEED NUMBER.          *
!    *                                                             *
!    * RANF GIVES ONE RANDOM NUMBER BETWEEN 0 AND 1.               *
!    * IRN55 GENERATES 55 RANDOM NUMBERS BETWEEN 0 AND 1/FMODUL.   *
!    * IN55  INITIALIZES THE 55 NUMBERS AND WARMS UP THE SEQUENCE. *
!    ***************************************************************
  use precision
  implicit none
  REAL (pr) :: RN
  INTEGER   :: IDMY   ! Dummy argument parameter
  INTEGER, SAVE :: IA (55), NCALL, MCALL
  REAL (pr), PARAMETER :: FMODUL = 1.E-09_pr
  DATA NCALL/0/
  DATA MCALL/55/

  IF (NCALL.EQ.0) THEN
     CALL IN55 (IA, 234612947)
     NCALL = 1
  ENDIF
  IF (MCALL.EQ.0) THEN
     CALL IRN55 (IA)
     MCALL=55
  ENDIF
  RN=IA (MCALL)*FMODUL
  MCALL=MCALL-1
END FUNCTION RN

!--------------------------------------------------------------------------
SUBROUTINE IN55 (IA, IX)
!--------------------------------------------------------------------------
  implicit none
  INTEGER :: IA (55), IX
  INTEGER, PARAMETER :: MODULO = 1000000000
  INTEGER :: J, K, I, II ! I is a dummy loop variable

  IA (55)=IX
  J=IX
  K=1
  DO I=1, 54
     II=MOD (21*I, 55)
     IA (II)=K
     K=J-K
     IF (K.LT.0)K=K+MODULO
     J=IA (II)
  ENDDO
  DO I=1, 10
     CALL IRN55 (IA)
  ENDDO
  RETURN
END SUBROUTINE IN55

!--------------------------------------------------------------------------
SUBROUTINE IRN55 (IA)
!--------------------------------------------------------------------------
  implicit none
  INTEGER IA (55)
  INTEGER, PARAMETER :: MODULO = 1000000000
  INTEGER I, J ! I is a dummy loop variable
  DO I=1, 24
     J=IA (I)-IA (I+31)
     IF (J.LT.0)J=J+MODULO
     IA (I)=J
  ENDDO
  DO I=25, 55
     J=IA (I)-IA (I-24)
     IF (J.LT.0)J=J+MODULO
     IA (I)=J
  ENDDO
  RETURN
END SUBROUTINE IRN55
!--------------------------------------------------------------------------
