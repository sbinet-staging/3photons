!
!   Common pour passer les résultats finaux
!
!> \brief array of final results: differential cross-section, sum and variance
!> \author Vincent C. LAFAGE
!> \date 2007-06-13 ISO
module resfin
  use precision
  implicit none
  integer, parameter :: NRES = 8 !< number of type of results
  REAL (pr) :: &
       SPM2DIF (NRES), & !< element of cross-section
       SPM2 (2, NRES), & !< total cross-section
       VAR (2, NRES)     !< variance of the sum
  !> \brief type replacing the previous common
  type resfin_type
     real (pr) :: &
          spm2dif (NRES), & !< element of cross-section
          spm2 (2, NRES), & !< total cross-section
          var (2, NRES)     !< variance of the sum
  end type resfin_type
end module resfin
