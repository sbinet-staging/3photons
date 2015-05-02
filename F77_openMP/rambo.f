C--------------------------------------------------------------------------
      SUBROUTINE RAMBO(N,ET,XM,P,WT)
C--------------------------------------------------------------------------
C
C                       RAMBO
C
C    RA(NDOM)  M(OMENTA)  B(EAUTIFULLY)  O(RGANIZED)
C
C ****** VERSION LIGHT & LIGHTSPEED ******
C
C    A DEMOCRATIC MULTI-PARTICLE PHASE SPACE GENERATOR
C    AUTHORS@D  S.D. ELLIS,  R. KLEISS,  W.J. STIRLING
C    THIS IS VERSION 1.0 -  WRITTEN BY R. KLEISS
C
C    N  = NUMBER OF PARTICLES (>1, IN THIS VERSION <101)
C    ET = TOTAL CENTRE-OF-MASS ENERGY
C    XM = PARTICLE MASSES ( DIM=100 )
C    P  = PARTICLE MOMENTA ( DIM=(4,100) )
C    WT = WEIGHT OF THE EVENT
C
C--------------------------------------------------------------------------
!     $ use OMP_LIB
      IMPLICIT NONE
      INTEGER N,I,K,NM,IBEGIN,NP
!      INTEGER ITMAX
      PARAMETER(NP=100)
      DOUBLE PRECISION ET,XMT,C,S,F,RMAS,G,A,X,BQ!,ACC
      DOUBLE PRECISION TWOPI,PO2LOG
      parameter (TWOPI=8.D0*ATAN (1.D0))
      parameter (PO2LOG=LOG (TWOPI/4.D0))
      DOUBLE PRECISION WT
      DOUBLE PRECISION XM(NP),P(4,NP),Q(4,NP),Z(NP),R(4),B(3)
      INTEGER IWARN(5)
      INTEGER iii

      DOUBLE PRECISION RN
      DOUBLE PRECISION grnd
      DOUBLE PRECISION RANDOM_NUMBER (4)
!      SAVE TWOPI, PO2LOG
      SAVE Z, XMT
      SAVE IBEGIN, NM
      integer rang
      data rang /0/
!      DATA ACC/1.D-14/
!      DATA ITMAX/6/
      DATA IBEGIN/0/,IWARN/5*0/
!$    integer OMP_GET_THREAD_NUM

C INITIALIZATION STEP@D FACTORIALS FOR THE PHASE SPACE WEIGHT
      IF(IBEGIN.EQ.0) THEN
        IBEGIN=1
!        TWOPI=8.D0*ATAN(1.D0)
!$     rang = OMP_GET_THREAD_NUM ()
        call sgrnd (rang+1) ! sgenrand ()
c        NORM=(TWOPI**(4-3*INP))/DBLE(ITOT)
!        PO2LOG=LOG(TWOPI/4.D0)
        Z(2)=PO2LOG
        DO K=3,NP
          Z(K)=Z(K-1)+PO2LOG-2.D0*LOG(DBLE(K-2))
!          Z(K)=Z(K-1)+PO2LOG-2.D0*LOG(DFLOAT(K-2))
        ENDDO
        DO K=3,NP
!         Z(K)=(Z(K)-LOG(DFLOAT(K-1)))
         Z(K)=(Z(K)-LOG(DBLE(K-1)))
        ENDDO

C CHECK ON THE NUMBER OF PARTICLES
        IF(N.GT.1.AND.N.LT.101) GOTO 104
        PRINT 1001,N
        STOP

C CHECK WHETHER TOTAL ENERGY IS SUFFICIENT; COUNT NONZERO MASSES
  104   XMT=0.D0
        NM=0
        DO I=1,N
          IF(XM(I).NE.0.D0) NM=NM+1
          XMT=XMT+ABS(XM(I))
        ENDDO
        IF(XMT.LE.ET) GOTO 201
        PRINT 1002,XMT,ET
        STOP
      ENDIF

C THE PARAMETER VALUES ARE NOW ACCEPTED

C GENERATE N MASSLESS MOMENTA IN INFINITE PHASE SPACE
  201 DO I=1,N
         do iii=1, 4
!            RANDOM_NUMBER (iii) = RN (iii)
            RANDOM_NUMBER (iii) = grnd()
         end do
        C=2.D0*RANDOM_NUMBER (1) - 1.D0
        S=SQRT (1.D0-C*C)
        F=TWOPI*RANDOM_NUMBER (2)
        Q (4,I) = -LOG(RANDOM_NUMBER (3) * RANDOM_NUMBER (4))
        Q (3,I) = Q (4,I) * C
        Q (2,I) = Q (4,I) * S *COS (F)
        Q (1,I) = Q (4,I) * S *SIN (F)
      ENDDO

C CALCULATE THE PARAMETERS OF THE CONFORMAL TRANSFORMATION
      DO I=1,4
        R(I)=0.D0
      ENDDO
      DO I=1,N
        DO K=1,4
          R(K)=R(K)+Q(K,I)
        ENDDO
      ENDDO
      RMAS=SQRT(R(4)**2-R(3)**2-R(2)**2-R(1)**2)
      DO K=1,3
        B(K)=-R(K)/RMAS
      ENDDO
      G=R(4)/RMAS
      A=1.D0/(1.D0+G)
      X=ET/RMAS

C TRANSFORM THE Q'S CONFORMALLY INTO THE P'S
      DO I=1,N
        BQ=B(1)*Q(1,I)+B(2)*Q(2,I)+B(3)*Q(3,I)
        DO K=1,3
          P(K,I)=X*(Q(K,I)+B(K)*(Q(4,I)+A*BQ))
        ENDDO
        P(4,I)=X*(G*Q(4,I)+BQ)
      ENDDO

C CALCULE LE POIDS ET LES AVERTISSEMENTS EVENTUELS
      WT=PO2LOG
      IF(N.NE.2) WT=(2.D0*N-4.D0)*DLOG(ET)+Z(N)
      IF(WT.GE.-180.D0) GOTO 208
      IF(IWARN(1).LE.5) PRINT 1004,WT
      IWARN(1)=IWARN(1)+1
  208 IF(WT.LE. 174.D0) GOTO 209
      IF(IWARN(2).LE.5) PRINT 1005,WT
      IWARN(2)=IWARN(2)+1

C RENVOIE LES IMPULSIONS SANS MASSES PONDEREES
 209  WT=EXP(WT)
      RETURN

 1001 FORMAT(' RAMBO FAILS@D # OF PARTICLES =',I5,' IS NOT ALLOWED')
 1002 FORMAT(' RAMBO FAILS@D TOTAL MASS =',D15.6,' IS NOT',
     . ' SMALLER THAN TOTAL ENERGY =',D15.6)
 1004 FORMAT(' RAMBO WARNS@D WEIGHT = EXP(',F20.9,') MAY UNDERFLOW')
 1005 FORMAT(' RAMBO WARNS@D WEIGHT = EXP(',F20.9,') MAY  OVERFLOW')

      END
