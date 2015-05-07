!
!   Common pour passer les produits spinoriels engendres par mat3
!
module spinor
  use precision
  implicit none
  COMPLEX (pr) :: S (5, 5), T (5, 5)
  REAL (pr) :: RAC8
  !REAL (pr), parameter :: RAC8 = 2_pr * sqrt (2._pr)
  type spinor_typ
     complex (pr) :: s (5, 5), t (5, 5)
  end type spinor_typ
end module spinor
