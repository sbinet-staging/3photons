!
!   Common pour passer les parametres physiques
!
module param
  use precision
  implicit none
  REAL (pr) :: MZ0,GZ0,GA,GBP,GBM,POLP,POLM,POLP2,POLM2
  LOGICAL IMPR
  type param_typ
     real (pr) :: mz0, gz0, ga, gbp, gbm, polp, polm, polp2, polm2
     logical impr
  end type param_typ
end module param
