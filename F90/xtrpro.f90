!
!   Common pour passer les probabilites supplementaires
!
module xtrpro
  use precision
  implicit none
  REAL*4 :: PRPLUS,PRMOINS
  REAL (pr) :: EE1,EE2
  type xtrpro_typ
     real*4 :: prplus, prmoins
     real (pr) :: ee1, ee2
  end type xtrpro_typ
end module xtrpro
! present dans mc.f (naturellement)
! present dans book.f
