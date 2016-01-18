!
!   Common pour passer les impulsions engendrees par Rambo
!
!> \brief array of generated momenta
!> \author Vincent C. LAFAGE
!> \date 2007-06-13 ISO
module ppp
  use precision
  implicit none
  INTEGER, PARAMETER :: INP = 3    !< Number of generated impulsions
  REAL (pr) :: &
       POUT (4, 100), & !< array of outgoing 4-momenta
       P1 (4),        & !< incoming (electron) 4-momentum
       P2 (4)           !< incoming (positron) 4-momentum
  !> \brief type replacing the previous common
  type ppp_typ
     real (pr) :: &
          pout (4, 100), & !< array of outgoing 4-momenta   
          p1 (4),        & !< incoming (electron) 4-momentum
          p2 (4)           !< incoming (positron) 4-momentum
  end type ppp_typ
end module ppp
