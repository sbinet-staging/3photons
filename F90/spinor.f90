!
!   Common pour passer les produits spinoriels engendres par mat3
!
!> \brief array of massless momenta spinor products
!> \author Vincent C. LAFAGE
!> \date 2007-06-13 ISO
module spinor
  use precision
  implicit none
  COMPLEX (pr), dimension (5, 5)  :: &
       S, & !< massless momenta spinor inner products Gram matrix
       T    !< massless momenta conjugate spinor inner products Gram matrix

  REAL (pr) :: RAC8 !< √8
  !REAL (pr), parameter :: RAC8 = 2_pr * sqrt (2._pr)
  !> \brief type replacing the previous common
  type spinor_typ
     complex (pr), dimension (5, 5) :: &
          s, & !< massless momenta spinor products Gram matrix
          t    !< massless momenta conjugate spinor products Gram matrix
  end type spinor_typ
end module spinor
