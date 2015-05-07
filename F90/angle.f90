!
!   Common pour passer les angles
!
module angle
  use precision
  implicit none
  REAL (pr) :: COS12,COS1P1,COS2P1,COS3P1,COS13,COS23,COSN,COSAC
  type angle_typ
     real (pr) :: cos12, cos1p1, cos2p1, cos3p1, cos13, cos23, cosn, cosac
  end type angle_typ
end module angle
