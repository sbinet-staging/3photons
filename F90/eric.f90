!----------------------------------------------------------------------
subroutine ERIC (UNE, BREPEM, CONVERS, PI)
!----------------------------------------------------------------------
!
! Affiche les resultats dans le parametrage d'Eric
!
!----------------------------------------------------------------------
  use precision
  use param
  use ppp
  use resfin
  implicit none
  INTEGER UNE
  REAL (pr) :: BREPEM, CONVERS, PI
  REAL (pr) :: MUTH

  MUTH=BREPEM/(8*9*5*PI**2*MZ0*GZ0)*CONVERS
  WRITE (UNE, *)
  WRITE (UNE, *)   '       :        -          +'
  WRITE (UNE, 950) 'sigma0  : ', SPM2(1, 1)/2._pr, SPM2(2, 1)/2._pr
  WRITE (UNE, 950) 'alpha0  : ', SPM2(1, 5)/2._pr, SPM2(2, 5)/2._pr
  WRITE (UNE, 950) 'beta0   : ', -SPM2(1, 4)/2._pr, -SPM2(2, 4)/2._pr
  WRITE (UNE, 950) 'lambda0 : ', (-SPM2(1, 2)+SPM2(1, 3))/2._pr, &
                                 (-SPM2(2, 2)+SPM2(2, 3))/2._pr
  WRITE (UNE, 950) 'mu0     : ', (SPM2(1, 2)+SPM2(1, 3))/2._pr, &
                                 (SPM2(2, 2)+SPM2(2, 3))/2._pr
  WRITE (UNE, 950) 'mu/lamb : ',  &
                     (SPM2(1, 2)+SPM2(1, 3))/(-SPM2(1, 2)+SPM2(1, 3)), &
                     (SPM2(2, 2)+SPM2(2, 3))/(-SPM2(2, 2)+SPM2(2, 3))
  WRITE (UNE, 950) 'mu (th) : ', MUTH
  WRITE (UNE, 950) 'mu (num): ',  &
           (SPM2(1, 2)+SPM2(1, 3)+SPM2(2, 2)+SPM2(2, 3))/4._pr
  WRITE (UNE, 950) 'rapport : ',  &
          (SPM2(1, 2)+SPM2(1, 3)+SPM2(2, 2)+SPM2(2, 3))/4._pr/MUTH

  RETURN
950 FORMAT(A, 2E12.4)
END SUBROUTINE ERIC
!----------------------------------------------------------------------
