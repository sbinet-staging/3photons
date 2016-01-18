!
!   Common pour passer les resultats engendrees par mat(rice)
!
!> \brief array of squared matrix elements contribution with detail of helicities
!> \author Vincent C. LAFAGE
!> \date 2007-06-13 ISO
module result
  use precision
  implicit none
  integer, parameter :: NRESUL = 8 !< number of type of results
  REAL (pr) :: M2 (2, 2, 2, NRESUL)     !< array of squared matrix elements NRESUL contribution with detail of outgoing helicities configuration
  !> \brief type replacing the previous common
  type result_typ
     real (pr) :: m2 (2, 2, 2, NRESUL)  !< array of squared matrix elements NRESUL contribution with detail of outgoing helicities configuration
  end type result_typ
end module result
