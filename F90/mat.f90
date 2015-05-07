!----------------------------------------------------------------------
!
! Calcule l'élément de matrice e+e- en 3 photons
!
!----------------------------------------------------------------------
SUBROUTINE MATRX (ETOT)
  use precision
  use scalar
  use result
  use param
  implicit none
  REAL (pr) :: ETOT
  INTEGER l1, l2, l3
  INTEGER LL1, LL2, LL3
  INTEGER K
  COMPLEX (pr) :: A, BP, BM
  COMPLEX (pr) :: BPPP, BPPM, BPMM, BMMM, APPM, APMM

! LES - VIENNENT DES Photons SORTANTS
  DO LL1 = 1, 2
     l1 = - (2*LL1-3)
     DO LL2 = 1, 2
        l2 = - (2*LL2-3)
        DO LL3 = 1, 2
           l3 = - (2*LL3-3)
! CALCULS DES AMPLITUDES D'HELICITE
           IF  ((l1.EQ.1).AND. (l2.EQ.1).AND. (l3.EQ.1)) THEN
              A  = 0._pr
              BP = 0._pr
              BM = BPPP (3, 4, 5)
           ENDIF
           IF  ((l1.EQ.-1).AND. (l2.EQ.-1).AND. (l3.EQ.-1)) THEN
              A  = 0._pr
              BP = 0._pr
              BM = BMMM (3, 4, 5)
           ENDIF
           IF  ((l1.EQ.1).AND. (l2.EQ.1).AND. (l3.EQ.-1)) THEN
              A  = APPM (3, 4, 5)
              BP = BPPM (3, 4, 5)
              BM = 0._pr
           ENDIF
           IF  ((l1.EQ.1).AND. (l2.EQ.-1).AND. (l3.EQ.1)) THEN
              A  = APPM (5, 3, 4)
              BP = BPPM (5, 3, 4)
              BM = 0._pr
           ENDIF
           IF  ((l1.EQ.-1).AND. (l2.EQ.1).AND. (l3.EQ.1)) THEN
              A  = APPM (4, 5, 3)
              BP = BPPM (4, 5, 3)
              BM = 0._pr
           ENDIF
           IF  ((l1.EQ.1).AND. (l2.EQ.-1).AND. (l3.EQ.-1)) THEN
              A  = APMM (3, 4, 5)
              BP = BPMM (3, 4, 5)
              BM = 0._pr
           ENDIF
           IF  ((l1.EQ.-1).AND. (l2.EQ.-1).AND. (l3.EQ.1)) THEN
              A  = APMM (5, 3, 4)
              BP = BPMM (5, 3, 4)
              BM = 0._pr
           ENDIF
           IF  ((l1.EQ.-1).AND. (l2.EQ.1).AND. (l3.EQ.-1)) THEN
              A  = APMM (4, 3, 5)
              BP = BPMM (4, 3, 5)
              BM = 0._pr
           ENDIF
! COUPLAGES
           A  = A *GA
           BP = BP*GBP
           BM = BM*GBM
! LES DIFFERENTS TERMES DE L'ELEMENT DE MATRICE AU CARRE
           M2 (LL1, LL2, LL3, 1) = ABS (A )**2
           M2 (LL1, LL2, LL3, 2) = ABS (BP)**2
           M2 (LL1, LL2, LL3, 3) = ABS (BM)**2
           M2 (LL1, LL2, LL3, 4) = 2.*REAL (A *CONJG (BP))
           M2 (LL1, LL2, LL3, 5) = 2.*IMAG (A *CONJG (BP))

!           M2 (LL1, LL2, LL3, 1) = CDABS (A )**2
!           M2 (LL1, LL2, LL3, 2) = CDABS (BP)**2
!           M2 (LL1, LL2, LL3, 3) = CDABS (BM)**2
!           M2 (LL1, LL2, LL3, 4) = 2.*DREAL (A *DCONJG (BP))
!           M2 (LL1, LL2, LL3, 5) = 2.*DIMAG (A *DCONJG (BP))

        ENDDO
     ENDDO
  ENDDO

  IF  (IMPR) THEN
     DO K = 1, NRESUL
        WRITE (*, *) K
        WRITE (*, *)'  ---     --+     -+-     +--     '// &
        '-++     +-+     ++-     +++'
          WRITE (*, 900) &
               M2 (1, 1, 1, K), M2 (1, 1, 2, K),  &
               M2 (1, 2, 1, K), M2 (2, 1, 1, K),  &
               M2 (1, 2, 2, K), M2 (2, 1, 2, K),  &
               M2 (2, 2, 1, K), M2 (2, 2, 2, K)
     ENDDO
  ENDIF

900 FORMAT (8E8.2)

  RETURN
END SUBROUTINE MATRX
!----------------------------------------------------------------------

!----------------------------------------------------------------------
SUBROUTINE PRECALC
!----------------------------------------------------------------------
  use precision
  use ppp
  use spinor
  use scalar
  use angle
  implicit none
  INTEGER J, K
  REAL (pr) :: P (5, 4)
  REAL (pr) :: NX, NY, NZ, NN
  REAL (pr) :: TX, XX (5)
  COMPLEX (pr) :: FX (5), CX

  DO J = 1, 4
     P (1, J) = P1 (J)
     P (2, J) = P2 (J)
     P (3, J) = POUT (J, 1)
     P (4, J) = POUT (J, 2)
     P (5, J) = POUT (J, 3)
  ENDDO
  DO K = 1, 5
     S (K, K) = 0._pr
     PS (K, K) = 0._pr
     TX = SQRT (P (K, 4)+P (K, 3))
     XX (K) = TX
     FX (K) = CMPLX (P (K, 1), P (K, 2))/TX
  ENDDO

  DO J = 1, 4
     DO K = J+1, 5
! PRODUIT SPINORIEL DE M.MANGANO,  S.PARKE
        CX = FX (J)*XX (K)-FX (K)*XX (J)
        S (J, K) = CX
        S (K, J) = -CX
        T (K, J) = CONJG (CX)
        T (J, K) = -CONJG (CX)
! PRODUIT SCALAIRE DE LORENTZ
        PS (J, K) = CX*CONJG (CX)/2.0_pr
        PS (K, J) = PS (J, K)
     ENDDO
  ENDDO

  COS1P1 = 1._pr-PS (1, 3)/ (P (1, 4)*P (3, 4))
  COS2P1 = 1._pr-PS (1, 4)/ (P (1, 4)*P (4, 4))
  COS3P1 = 1._pr-PS (1, 5)/ (P (1, 4)*P (5, 4))

  COS12 = 1._pr-PS (3, 4)/ (POUT (4, 1)*POUT (4, 2))
  COS13 = 1._pr-PS (3, 5)/ (POUT (4, 1)*POUT (4, 3))
  COS23 = 1._pr-PS (4, 5)/ (POUT (4, 2)*POUT (4, 3))

  NX = POUT (2, 1)*POUT (3, 2)-POUT (2, 2)*POUT (3, 1)
  NY = POUT (3, 1)*POUT (1, 2)-POUT (3, 2)*POUT (1, 1)
  NZ = POUT (1, 1)*POUT (2, 2)-POUT (1, 2)*POUT (2, 1)
  NN = SQRT (NX**2+NY**2+NZ**2)
  COSN =  (NX*P (1, 1)+NY*P (1, 2)+NZ*P (1, 3))/P (1, 4)/NN

  COSAC =  (P (3, 2)*P (4, 2)+P (3, 3)*P (4, 3))/SQRT ( &
       (P (3, 2)*P (3, 2)+P (3, 3)*P (3, 3))* &
       (P (4, 2)*P (4, 2)+P (4, 3)*P (4, 3)))

  RETURN
END SUBROUTINE PRECALC
!----------------------------------------------------------------------
