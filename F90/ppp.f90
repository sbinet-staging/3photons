!
!   Common pour passer les impulsions engendrees par Rambo
!
module ppp
  use precision
  implicit none
  INTEGER, PARAMETER :: INP = 3    ! Number of impulsions generated
  REAL (pr) :: POUT(4,100),P1(4),P2(4)
  type ppp_typ
     real (pr) :: pout (4, 100), p1 (4), p2 (4)
  end type ppp_typ
end module ppp
