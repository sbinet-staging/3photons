!
!   Common pour passer les coupures
!
module cutpar
  use precision
  implicit none
  REAL (pr) :: ACUT,BCUT,EMIN,SINCUT
  type cutpar_typ
     real (pr) :: acut, bcut, emin, sincut
  end type cutpar_typ
end module cutpar
