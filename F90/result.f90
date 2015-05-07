!
!   Common pour passer les resultats engendrees par mat(rice)
!
module result
  use precision
  implicit none
  integer, parameter :: NRESUL = 8
  REAL (pr) :: M2(2,2,2,NRESUL)
  type result_typ
     real (pr) :: m2 (2, 2, 2, NRESUL)
  end type result_typ
end module result
