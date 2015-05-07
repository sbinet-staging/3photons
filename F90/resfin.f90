!
!   Common pour passer les résultats finaux
!
module resfin
  use precision
  implicit none
  integer, parameter :: NRES = 8
  REAL (pr) :: SPM2DIF(NRES),SPM2(2,NRES),VAR(2,NRES)
  type resfin_type
     real (pr) :: spm2dif (NRES), spm2 (2, NRES), var (2, NRES)
  end type resfin_type
end module resfin
