!
!   Common pour passer les probabilites supplementaires
!
!> \brief set of extra probabilities (for histogramming)
!> \author Vincent C. LAFAGE
!> \date 2007-06-13 ISO
module xtrpro
  use precision
  implicit none
  REAL*4 :: PRPLUS, PRMOINS     !<  extra probabilities (for histogramming)
  REAL (pr) :: EE1,EE2
  !> \brief type replacing the previous common
  type xtrpro_typ
     real*4 :: prplus, prmoins  !<  extra probabilities (for histogramming)
     real (pr) :: ee1, ee2
  end type xtrpro_typ
end module xtrpro
! present dans mc.f (naturellement)
! present dans book.f
