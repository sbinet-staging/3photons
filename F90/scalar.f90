!
!   Common pour passer les produits scalaires engendres par mat3
!
!> \brief array of momenta Lorentz scalar products
!> \author Vincent C. LAFAGE
!> \date 2007-06-13 ISO
module scalar
  use precision
  implicit none
  REAL (pr), dimension (5, 5) :: PS !< symmetric array with all possible Lorentz 4-scalar products for the problem (Gram matrix)
  !> \brief type replacing the previous common
  type scalar_typ
     real (pr), dimension (5,5) :: ps !< symmetric array with all possible Lorentz 4-scalar product for the problem (Gram matrix)
  end type scalar_typ
end module scalar
