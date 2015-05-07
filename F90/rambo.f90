!--------------------------------------------------------------------------
SUBROUTINE RAMBO (N, ET, XM, P, WT)
!--------------------------------------------------------------------------
!
!                       RAMBO
!
!    RA(NDOM)  M(OMENTA)  B(EAUTIFULLY)  O(RGANIZED)
!
! ****** VERSION LIGHT & LIGHTSPEED ******
!
!    A DEMOCRATIC MULTI-PARTICLE PHASE SPACE GENERATOR
!    AUTHORS@D  S.D. ELLIS,  R. KLEISS,  W.J. STIRLING
!    THIS IS VERSION 1.0 -  WRITTEN BY R. KLEISS
!
!    N  = NUMBER OF PARTICLES (>1, IN THIS VERSION <101)
!    ET = TOTAL CENTRE-OF-MASS ENERGY
!    XM = PARTICLE MASSES ( DIM=100 )
!    P  = PARTICLE MOMENTA ( DIM=(4, 100) )
!    WT = WEIGHT OF THE EVENT
!
!--------------------------------------------------------------------------
  use precision
  implicit none
  INTEGER N
  INTEGER :: I, K, ITMAX = 6
  INTEGER, PARAMETER :: NP = 100
  INTEGER, SAVE      :: IBEGIN = 0, NM
  REAL (pr) :: ET
  REAL (pr),SAVE :: TWOPI, PO2LOG, Z(NP), XMT
  REAL (pr) :: C, S, F, RMAS, G, A, X, BQ
  REAL (pr) :: ACC = 1.0e-14_pr
  REAL (pr) :: WT
  REAL (pr) :: XM(NP), P(4, NP), Q(4, NP), R(4), B(3)
  INTEGER IWARN(5)
  DATA IWARN/5*0/
  REAL (pr) :: RN

! INITIALIZATION STEP@D FACTORIALS FOR THE PHASE SPACE WEIGHT
  IF(IBEGIN.EQ.0) THEN
     IBEGIN = 1
     TWOPI  = 8.0_pr*ATAN(1.0_pr)
!        NORM=(TWOPI**(4-3*INP))/DBLE(ITOT)
     PO2LOG = LOG (TWOPI / 4.0_pr)
     Z (2) = PO2LOG
     DO K = 3, NP
        Z (K) = Z (K-1) + PO2LOG -2.0_pr * LOG (DBLE (K-2))
!          Z(K)=Z(K-1)+PO2LOG-2_pr*LOG(DFLOAT(K-2))
     ENDDO
     DO K = 3, NP
        Z (K) = (Z (K) - LOG (DBLE (K-1)))
     ENDDO

! CHECK ON THE NUMBER OF PARTICLES
     IF(N.GT.1.AND.N.LT.101) GOTO 104
     PRINT 1001, N
     STOP

! CHECK WHETHER TOTAL ENERGY IS SUFFICIENT; COUNT NONZERO MASSES
104  XMT=0._pr
     NM=0
     DO I = 1, N
        IF (XM(I).NE.0._pr) NM = NM + 1
        XMT = XMT + ABS (XM (I))
     ENDDO
     IF (XMT.LE.ET) GOTO 201
     PRINT 1002, XMT, ET
     STOP
  ENDIF

! THE PARAMETER VALUES ARE NOW ACCEPTED

! GENERATE N MASSLESS MOMENTA IN INFINITE PHASE SPACE
201 DO I=1, N
     C=2._pr*RN(1)-1._pr
     S=SQRT(1._pr-C*C)
     F=TWOPI*RN(2)
     Q(4, I)=-LOG(RN(3)*RN(4))
     Q(3, I)=Q(4, I)*C
     Q(2, I)=Q(4, I)*S*COS(F)
     Q(1, I)=Q(4, I)*S*SIN(F)
     !write (*, *) 'C, F, Q (4, I): ', C, F, Q (4, I)
  ENDDO

! CALCULATE THE PARAMETERS OF THE CONFORMAL TRANSFORMATION
  DO I = 1, 4
     R(I)=0._pr
  ENDDO
  DO I = 1, N
     DO K = 1, 4
        R (K) = R (K) + Q (K, I)
     ENDDO
  ENDDO
  RMAS = SQRT (R (4)**2 - R (3)**2 - R (2)**2 - R (1)**2)
  DO K = 1, 3
     B (K) = -R (K) / RMAS
  ENDDO
  G=R(4)/RMAS
  A=1._pr/(1._pr+G)
  X=ET/RMAS

! TRANSFORM THE Q'S CONFORMALLY INTO THE P'S
  DO I=1, N
     BQ=B(1)*Q(1, I)+B(2)*Q(2, I)+B(3)*Q(3, I)
     DO K=1, 3
        P(K, I)=X*(Q(K, I)+B(K)*(Q(4, I)+A*BQ))
     ENDDO
     P(4, I)=X*(G*Q(4, I)+BQ)
  ENDDO

! CALCULE LE POIDS ET LES AVERTISSEMENTS EVENTUELS
  WT = PO2LOG
  IF (N.NE.2) WT = (2._pr * N - 4._pr) * LOG (ET) + Z (N)
  IF (WT.GE.-180._pr) GOTO 208
  IF (IWARN(1).LE.5) PRINT 1004, WT
  IWARN (1) = IWARN (1) + 1
208 IF (WT.LE. 174._pr) GOTO 209
  IF (IWARN (2).LE.5) PRINT 1005, WT
  IWARN (2) = IWARN (2)+1

! RENVOIE LES IMPULSIONS SANS MASSES PONDEREES
209 WT=EXP(WT)
  RETURN

1001 FORMAT(' RAMBO FAILS@D # OF PARTICLES =', I5, ' IS NOT ALLOWED')
1002 FORMAT(' RAMBO FAILS@D TOTAL MASS =', D15.6, ' IS NOT', &
            ' SMALLER THAN TOTAL ENERGY =', D15.6)
1004 FORMAT(' RAMBO WARNS@D WEIGHT = EXP(', F20.9, ') MAY UNDERFLOW')
1005 FORMAT(' RAMBO WARNS@D WEIGHT = EXP(', F20.9, ') MAY  OVERFLOW')

END SUBROUTINE RAMBO
