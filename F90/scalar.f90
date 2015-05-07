!
!   Common pour passer les produits scalaires engendres par mat3
!
module scalar
  use precision
  implicit none
  REAL (pr) :: PS(5,5)
  type scalar_typ
     real (pr) :: ps (5, 5)
  end type scalar_typ
end module scalar
